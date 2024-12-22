#include <string.h>
#include <stdio.h>
#include "rs_debug.h"

// Simple helper function for Galois Field modulo operations
static inline void MOD_NN(int& x) {
    while (x < 0) x += nn;
    x = x % nn;
}

void RS_ENCODER_REF::Init(void) {
    if (!bInitialized) {
        bInitialized = true;
    }
}

RS_ENCODER_REF::RS_ENCODER_REF(int CorrectableErrors) {
    // Initialize key parameters
    tt = CorrectableErrors;  // Number of correctable errors
    kk = nn - 2*tt;         // Number of data symbols (message length)

    // Choose generator polynomial starting power
    // For RS codes, we can choose b0 to make the generator polynomial symmetric
    // This doesn't affect correctness but can simplify some implementations
    b0 = (kk+1)/2;

    // Generate the generator polynomial
    RSGenPoly();
    RsVerify::verify_generator(gg, tt);
}

RS_ENCODER_REF::~RS_ENCODER_REF(void) {
}

void RS_ENCODER_REF::RSGenPoly(void) {
    // This function generates the generator polynomial for the RS code
    // The generator polynomial g(x) is the product of (x + α^i) terms
    // where i ranges from b0 to b0+2t-1

    int i, j, iMod;

    // Initialize first term (x + α^b0)
    gg[0] = Pow2Poly[b0];  // Constant term
    gg[1] = 1;             // Coefficient of x

    // Multiply by subsequent terms (x + α^(b0+i))
    for (i = 2; i <= 2*tt; i++) {
        gg[i] = 1;

        // Multiply current polynomial by (x + α^(b0+i-1))
        for (j = i-1; j > 0; j--) {
            if (gg[j] != 0) {
                iMod = Poly2Pow[gg[j]] + b0 + i-1;
                MOD_NN(iMod);
                gg[j] = gg[j-1] ^ Pow2Poly[iMod];
            } else {
                gg[j] = gg[j-1];
            }
        }

        // Handle constant term
        iMod = Poly2Pow[gg[0]] + b0 + i-1;
        MOD_NN(iMod);
        gg[0] = Pow2Poly[iMod];
    }

    // Convert generator polynomial coefficients to power form for encoding
    for (i = 0; i <= 2*tt; i++) {
        gg[i] = Poly2Pow[gg[i]];
    }
}

void RS_ENCODER_REF::RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]) {
    // Systematic encoder for RS codes
    // Input: data[0..kk-1] contains the message
    // Output: bb[0..2*tt-1] contains the parity bytes
    //
    // The codeword is formed by:
    // c(x) = m(x)*x^(n-k) + b(x)
    // where m(x) is the message polynomial and b(x) is the remainder when
    // m(x)*x^(n-k) is divided by g(x)

    // Clear parity buffer
    memset(bb, 0, 2*tt);

    // For each message byte
    for (int i = 0; i < kk; i++) {
        // Calculate feedback term
        GF feedback = bb[0] ^ data[i];

        // If feedback is non-zero, multiply feedback by each generator coefficient
        // and add to the shift register
        if (feedback != 0) {
            for (int j = 0; j < 2*tt-1; j++) {
                if (gg[j] != GF_INFINITY) {
                    int iMod = gg[j] + Poly2Pow[feedback];
                    MOD_NN(iMod);
                    bb[j] = bb[j+1] ^ Pow2Poly[iMod];
                } else {
                    bb[j] = bb[j+1];
                }
            }
            bb[2*tt-1] = feedback;
        } else {
            // Just shift if feedback is zero
            for (int j = 0; j < 2*tt-1; j++) {
                bb[j] = bb[j+1];
            }
            bb[2*tt-1] = 0;
        }
    }
}

int RS_ENCODER_REF::RSDecode(GF recd[nn]) {
    RsVerify::verify_received_word(recd, nn);

    int i, j, q, iMod;
    int u = 0;  // We already had this

    // Allocate working buffers
    GF syndromes[2*MAX_TT+1];    // Syndromes
    GF z[MAX_TT+2];              // must be at least tt+1
    GF root[MAX_TT];             // Error locations in power form
    GF loc[MAX_TT];              // Error position indices
    GF err[n];                   // Error values
    GF reg[2*MAX_TT+1];          // Scratch space for Chien search
    int syn_error = 0;           // Flag for non-zero syndrome
    int count = 0;               // Number of errors found

    // Step 1: Calculate Syndromes
    // For an error-free codeword, evaluating recd(x) at αⁱ (i=b0...b0+2t-1) should give zero
    memset(syndromes, 0, sizeof(syndromes));

    // For each position in received word
    int iPowInit = 0;
    for (int j = 0; j < nn; j++) {
        GF symbol = recd[j];
        if (symbol != 0) {
            // For each syndrome, calculate recd(αⁱ)
            int iPow0 = iPowInit + Poly2Pow[symbol];
            MOD_NN(iPow0);

            for (int i = 1; i <= 2*tt; i+=2) {
                MOD_NN(iPow0);
                syndromes[i] ^= Pow2Poly[iPow0];
                iPow0 += j;
                MOD_NN(iPow0);
                syndromes[i+1] ^= Pow2Poly[iPow0];
                iPow0 += j;
            }
        }
        iPowInit += b0;
        MOD_NN(iPowInit);
    }

    // Convert syndromes to power form and check for errors
    for (int i = 1; i <= 2*tt; i++) {
        if (syndromes[i] != 0) {
            syn_error = 1;
        }
        syndromes[i] = Poly2Pow[syndromes[i]];
    }

    RsVerify::verify_syndromes(syndromes, tt);
    RsDebug::print_syndromes("Initial", syndromes, tt);

    if (syn_error) {  // If errors, attempt to correct
        // Initialize arrays for Berlekamp-Massey
        GF d[2*MAX_TT+2];        // Like rs.cpp
        int l[2*MAX_TT+2];       // Degree of each elp
        int u_lu[2*MAX_TT+2];    // Difference between step number and degree
        GF elp[2*MAX_TT+2][2*MAX_TT];  // Error locator polynomials

        // Initial conditions
        d[0] = 0;           // power form
        d[1] = syndromes[1];  // power form
        elp[0][0] = 0;      // power form
        elp[1][0] = 1;      // polynomial form
        for (i = 1; i < nn-kk; i++) {
            elp[0][i] = GF_INFINITY;   // power form
            elp[1][i] = 0;   // polynomial form
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;

        do {
            u++;
            if (d[u] == GF_INFINITY) {
                l[u+1] = l[u];
                for (i = 0; i <= l[u]; i++) {
                    elp[u+1][i] = elp[u][i];
                    elp[u][i] = Poly2Pow[elp[u][i]];
                }
            }
            else {
                // Search for words with greatest u_lu[q] for which d[q]!=0
                q = u-1;
                while ((d[q] == GF_INFINITY) && (q > 0))
                    q--;

                /* have found first non-zero d[q]  */
                if (q>0)
                {
                    j=q;
                    do
                    {
                        j--;
                        if ((d[j] != GF_INFINITY) && (u_lu[q]<u_lu[j]))
                            q = j;
                    } while (j>0);
                }

                /* have now found q such that d[u]!=0 and u_lu[q] is maximum */
                /* store degree of new elp polynomial */
                if (l[u] > l[q]+u-q)
                    l[u+1] = l[u];
                else
                    l[u+1] = l[q]+u-q;

                /* form new elp(x) */
                for (i = 0; i < nn-kk; i++)
                    elp[u+1][i] = 0;

                for (i=0; i<=l[q]; i++)
                    if (elp[q][i] != GF_INFINITY)
                    {
                        iMod = d[u]-d[q]+elp[q][i];
                        MOD_NN(iMod);
                        elp[u+1][i+u-q] = Pow2Poly[iMod];
                    }
                for (i=0; i<=l[u]; i++)
                {
                    elp[u+1][i] ^= elp[u][i];
                    elp[u][i] = Poly2Pow[elp[u][i]];  /*convert old elp value to power*/
                }
            }
            u_lu[u+1] = u-l[u+1];

            /* form (u+1)th discrepancy */
            if (u < nn-kk)    /* no discrepancy computed on last iteration */
            {
                d[u+1] = Pow2Poly[syndromes[u+1]];
                for (i = 1; i <= l[u+1]; i++)
                {
                    if ((syndromes[u+1-i] != GF_INFINITY) && (elp[u+1][i] != 0))
                    {
                        iMod = syndromes[u+1-i]+Poly2Pow[elp[u+1][i]];
                        MOD_NN(iMod);
                        d[u+1] ^= Pow2Poly[iMod];
                    }
                }
                d[u+1] = Poly2Pow[d[u+1]];    /* put d[u+1] into power form */
            }

            RsDebug::print_berlekamp_step(u, d[u], l[u], elp[u], l[u]);  // Updated debug
        } while ((u < nn-kk) && (l[u+1] <= tt));

        u++;
        RsVerify::verify_lambda(elp[u], l[u]);

        if (l[u] <= tt)         /* can correct error */
        {
            /* put elp into power form */
            for (i=0; i<=l[u]; i++)
                elp[u][i] = Poly2Pow[elp[u][i]];

            /* find roots of the error location polynomial */
            for (i=1; i <= l[u]; i++)
                reg[i] = elp[u][i];
            count = 0;
            for (i=1; i <= nn; i++)
            {
                q = 1;
                for (j = 1; j <= l[u]; j++)
                {
                    if (reg[j] != GF_INFINITY)
                    {
                        iMod = reg[j]+j;
                        MOD_NN(iMod);
                        q ^= Pow2Poly[iMod];
                        reg[j] = iMod;
                    }
                }
                if (!q)  /* store root and error location number indices */
                {
                    root[count] = i;
                    loc[count] = nn-i;
                    count++;
                }
            }
            printf("(count, l[u]) = (%d, %d)\n", count, l[u]);
            if (count == l[u])    /* no. roots = degree of elp hence <= tt errors */
            {
                /* form polynomial z(x) */
                for (i = 1; i <= l[u]; i++) /* Z[0] = 1 always - do not need */
                {
                    GF zi = Pow2Poly[syndromes[i]] ^ Pow2Poly[elp[u][i]];
                    for (j=1; j<i; j++)
                    {
                        if ((syndromes[j] != GF_INFINITY) && (elp[u][i-j] != GF_INFINITY))
                        {
                            iMod = elp[u][i-j] + syndromes[j];
                            MOD_NN(iMod);
                            zi ^= Pow2Poly[iMod];
                        }
                    }
                    z[i] = Poly2Pow[zi];         /* put into power form */
                }

                // Evaluate errors at locations given by error location
                // numbers loc[i]
                //
                for (i = 0; i < nn; i++)
                    err[i] = 0;

                // Compute numerator of error term first
                //
                for (i=0; i < l[u]; i++)
                {
                    int jPow = root[i];
                    err[loc[i]] = 1;       /* accounts for z[0] */
                    for (j=1; j<=l[u]; j++)
                    {
                        if (z[j] != GF_INFINITY)
                        {
                            iMod = z[j]+jPow;
                            MOD_NN(iMod);
                            err[loc[i]] ^= Pow2Poly[iMod];
                        }
                        jPow += root[i]; // jPow = root[i]*j;
                        MOD_NN(jPow);
                    } /* z(x) evaluated at X(l)^(-1) */

                    // term X(l)^(1-b0)
                    //
                    if (err[loc[i]] != 0)
                    {
                        err[loc[i]] = Pow2Poly[(Poly2Pow[err[loc[i]]]+root[i]*(b0+nn-1))%nn];
                    }
                    if (err[loc[i]] != 0)
                    {
                        err[loc[i]] = Poly2Pow[err[loc[i]]];
                        q = 0;     /* form denominator of error term */
                        for (j = 0; j < l[u]; j++)
                            if (j != i)
                            {
                                iMod = loc[j]+root[i];
                                MOD_NN(iMod);
                                q += Poly2Pow[1^Pow2Poly[iMod]];
                                MOD_NN(q);
                            }
                        // q = q % nn;
                        iMod = err[loc[i]]-q;
                        MOD_NN(iMod);
                        err[loc[i]] = Pow2Poly[iMod];
                        recd[loc[i]] ^= err[loc[i]];  /*recd[i] must be in polynomial form */
                        RsDebug::print_error_location(loc[i], err[loc[i]], recd[loc[i]] ^ err[loc[i]], recd[loc[i]]);
                    }
                }
                return count; // Number of corrections.
            }
            else
            {    /* no. roots != degree of elp => >tt errors and cannot solve */
                return -2; // NOT ENOUGH ROOTS BY CHIEN SEARCH
            }
        }
        else
        {         /* elp has degree has degree >tt hence cannot solve */
            return -3; // LAMDA DEGREE TOO HIGH
        }
    }
    else
    {
        /* no non-zero syndromes => no errors: output received codeword */
        return 0;
    }
}

bool RS_ENCODER_REF::bInitialized = false;
