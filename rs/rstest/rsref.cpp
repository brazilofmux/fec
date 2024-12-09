#include "rs.h"
#include <string.h>
#include <stdio.h>
#include "rs_debug.h"

// Tables for converting between power form and polynomial form
GF Pow2Poly[n], Poly2Pow[n];
#define PrimitivePolynomial  0x1D    // [x^8] + x^4 + x^3 + x^2 + 1

// Simple modulo function for GF operations
static inline void MOD_NN(int& x) {
    while (x < 0) x += nn;
    x = x % nn;
}

// Generate GF(2^8) from the primitive polynomial
static void RSGenField(void) {
    int iPoly = 1;

    // Handle special infinite case
    Pow2Poly[n-1] = 0;
    Poly2Pow[0] = GF_INFINITY;

    // Generate power/polynomial conversion tables
    for (int i = 0; i < n-1; i++) {
        Pow2Poly[i] = iPoly;
        Poly2Pow[iPoly] = i;
        iPoly <<= 1;
        if (iPoly > n-1) {
            iPoly ^= PrimitivePolynomial;
            iPoly &= 0xFF;
        }
    }
}

void RS_Init(void)
{
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        RSGenField();
	RsVerify::verify_tables();
    }
}

RS_ENCODER::RS_ENCODER(int CorrectableErrors) {
    tt = CorrectableErrors;
    kk = nn-2*tt;

    // Choose generator polynomial starting power
    b0 = (kk+1)/2;
    RSGenPoly();
    RsVerify::verify_generator(gg, tt);
}

RS_ENCODER::~RS_ENCODER(void) {
}

void RS_ENCODER::RSGenPoly(void) {
    // Generate the generator polynomial
    int i, j, iMod;

    gg[0] = Pow2Poly[b0];  // First root
    gg[1] = 1;

    // Multiply by subsequent factors
    for (i = 2; i <= 2*tt; i++) {
        gg[i] = 1;
        for (j = i-1; j > 0; j--) {
            if (gg[j] != 0) {
                iMod = Poly2Pow[gg[j]]+b0+i-1;
                MOD_NN(iMod);
                gg[j] = gg[j-1]^Pow2Poly[iMod];
            } else {
                gg[j] = gg[j-1];
            }
        }
        iMod = Poly2Pow[gg[0]]+b0+i-1;
        MOD_NN(iMod);
        gg[0] = Pow2Poly[iMod];
    }

    // Convert to power form for encoding
    for (i = 0; i <= 2*tt; i++) {
        gg[i] = Poly2Pow[gg[i]];
    }
}

void RS_ENCODER::RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]) {
    // Simple systematic encoder without optimization tricks
    memset(bb, 0, 2*tt);

    for (int i = 0; i < kk; i++) {
        GF feedback = bb[0]^data[i];
        for (int j = 0; j < 2*tt-1; j++) {
            if (gg[j] != GF_INFINITY && feedback != GF_INFINITY) {
                int iMod = gg[j]+Poly2Pow[feedback];
                MOD_NN(iMod);
                bb[j] = bb[j+1]^Pow2Poly[iMod];
            } else {
                bb[j] = bb[j+1];
            }
        }
        bb[2*tt-1] = feedback;
    }
}
int RS_ENCODER::RSDecode(GF recd[nn]) {
    RsDebug::init(true);  // Enable debugging
    RsVerify::verify_received_word(recd, nn);

    // Declare all variables up front
    int i, j, q, iMod;
    int u = 0;  // We already had this

    GF syndromes[2*MAX_TT+1];    // Syndromes
    GF lambda[2*MAX_TT+1];       // Error locator polynomial
    GF root[MAX_TT];             // Error locations
    GF loc[MAX_TT];              // Error position indices
    GF reg[2*MAX_TT+1];          // Scratch space for Chien search
    int deg_lambda = 0;          // Degree of error locator polynomial
    int syn_error = 0;
    int count = 0;               // Number of errors found

    // First form the syndromes: evaluate recd(x) at roots of g(x)
    memset(syndromes, 0, sizeof(syndromes));

    // For each position in received word
    int iPowInit = 0, iPow0;
    for (int j = 0; j < nn; j++) {
        GF RECD = recd[j];
        if (RECD != 0) {
            iPow0 = iPowInit + Poly2Pow[RECD];
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
        // Initialize for Berlekamp-Massey algorithm
        lambda[0] = 1;
        lambda[1] = Pow2Poly[syndromes[1]];
        for (int i = 2; i <= 2*tt; i++) {
            lambda[i] = 0;
        }

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

        // Convert lambda to power form
        for (int i = 0; i <= 2*tt; i++) {
            lambda[i] = Poly2Pow[lambda[i]];
        }

        // Find degree of lambda
        deg_lambda = 2*tt;
        while (deg_lambda > 0 && lambda[deg_lambda] == GF_INFINITY) {
            deg_lambda--;
        }
	RsVerify::verify_lambda(lambda, deg_lambda);

        if (deg_lambda <= 2*tt) {
            // Find roots of error locator polynomial using Chien search
            memcpy(reg, lambda, sizeof(reg));
            count = 0;

            for (int i = 1; i <= nn; i++) {
                GF q = 1;
                for (int j = 1; j <= deg_lambda; j++) {
                    if (reg[j] != GF_INFINITY) {
                        int iMod = reg[j] + j;
                        MOD_NN(iMod);
                        q ^= Pow2Poly[iMod];
                        reg[j] = iMod;
                    }
                }
                if (q == 0) {
                    root[count] = i;
                    loc[count] = nn - i;
                    count++;
                }
            }

            if (count == deg_lambda) {
                // Calculate error values
                for (int i = 0; i < count; i++) {
                    GF err = 0;
                    int pos = loc[i];
                    // Compute error value at position
                    for (int j = 0; j <= deg_lambda; j++) {
                        if (lambda[j] != GF_INFINITY) {
                            int iMod = j * root[i];
                            MOD_NN(iMod);
                            err ^= Pow2Poly[iMod];
                        }
                    }
                    if (err != 0) {
                        GF original = recd[pos];
                        recd[pos] ^= err;  // Correct the error
                        RsDebug::print_error_location(pos, err, original, recd[pos]);
                    }
                }
                return count;  // Return number of corrections
            }
            else {
                return -2;  // Number of roots != degree of lambda
            }
        }
        return -3;  // deg_lambda > tt
    }
    return 0;  // No errors found
}
