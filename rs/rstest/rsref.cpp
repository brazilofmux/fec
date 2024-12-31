#include <string.h>
#include <stdio.h>
#include <iostream>
#include "rs.h"
#include "rsref.h"
#include "RsVerification.h"

// Simple helper function for Galois Field modulo operations
static inline void MOD_NN(int& x) {
    while (x >= nn) {
        x -= nn;
        x = (x >> 8) + (x & nn);
    }
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
    RsVerification::verify_generator(gg, tt);
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
    // Clear parity buffer
    memset(bb, 0, 2 * tt);

    // Systematic encoder for RS codes
    // Input: data[0..kk-1] contains the message
    // Output: bb[0..2*tt-1] contains the parity bytes
    //
    // The codeword is formed by:
    // c(x) = m(x)*x^(n-k) + b(x)
    // where m(x) is the message polynomial and b(x) is the remainder when
    // m(x)*x^(n-k) is divided by g(x)
    for (int i = kk - 1; i >= 0; i--) {
        GF feedback = Poly2Pow[data[i] ^ bb[2 * tt - 1]];  // Get feedback using current data and parity

        if (feedback != GF_INFINITY) {  // Non-zero feedback term
            for (int j = 2 * tt - 1; j > 0; j--) {
                if (gg[j] != GF_INFINITY) {
                    int modVal = gg[j] + feedback;
                    MOD_NN(modVal);
                    bb[j] = bb[j - 1] ^ Pow2Poly[modVal];
                } else {
                    bb[j] = bb[j - 1];
                }
            }
            // Update the first parity byte with gg[0] term
            int modVal = gg[0] + feedback;
            MOD_NN(modVal);
            bb[0] = Pow2Poly[modVal];
        } else {  // Zero feedback: Shift the parity bytes without XOR
            for (int j = 2 * tt - 1; j > 0; j--) {
                bb[j] = bb[j - 1];
            }
            bb[0] = 0;  // Clear the first parity byte
        }
    }
}

int RS_ENCODER_REF::RSDecode(GF recd[nn]) {
    /* Reed-Solomon Decoding Process
     *
     * Mathematical Background:
     * 1. A Reed-Solomon code is defined by evaluating codeword polynomials at powers
     *    of a primitive element α of the Galois Field GF(2^8).
     *
     * 2. The generator polynomial g(x) has roots α^b0, α^(b0+1), ..., α^(b0+2t-1)
     *    where t is the error-correction capability.
     *
     * 3. For a received word r(x) = c(x) + e(x) where:
     *    - c(x) is the transmitted codeword
     *    - e(x) is the error polynomial
     *
     * 4. The decoding process:
     *    a) Calculate syndromes Si = r(α^(b0+i)) for i = 1..2t
     *    b) Use Berlekamp-Massey to find error locator polynomial Λ(x)
     *    c) Find roots of Λ(x) to locate errors (Chien search)
     *    d) Calculate error values using Forney algorithm
     */
    RsVerification::verify_received_word(recd, nn);

    // Allocate working buffers
    GF syndromes[2*MAX_TT+1];    // Syndromes
    GF root[MAX_TT];             // Error locations in power form
    GF loc[MAX_TT];              // Error position indices
    GF err[n];                   // Error values
    GF reg[2*MAX_TT+1];          // Scratch space for Chien search
    int syn_error = 0;           // Flag for non-zero syndrome
    int count = 0;               // Number of errors found

    /* Step 1: Syndrome Calculation
     *
     * Theory: For a valid codeword c(x), c(α^i) = 0 for i = b0..b0+2t-1
     * Therefore, for received word r(x):
     *   Si = r(α^(b0+i)) = e(α^(b0+i))
     *
     * These syndrome values depend only on the errors, not the transmitted codeword.
     * If all syndromes are zero, r(x) is a valid codeword.
     */
#if 1
    // Initialize syndromes to zero
    for (int i = 0; i <= 2*tt; i++) {
        syndromes[i] = 0;
    }

    // For each position in received word
    int iPowInit = 0;
    for (int j = 0; j < nn; j++) {
        GF RECD = recd[j];
        if (RECD != 0) {
            int iPow0 = iPowInit + Poly2Pow[RECD];
            MOD_NN(iPow0);

            for (int i = 1; i <= 2*tt; i++) {
                MOD_NN(iPow0);
                syndromes[i] ^= Pow2Poly[iPow0];
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
#else
    memset(syndromes, 0, sizeof(syndromes));

    // For each position in received word
    int iPowInit = 0;
    for (int j = 0; j < nn; j++) {
        GF symbol = recd[j];
        if (symbol != 0) {
            // For each syndrome, calculate recd(αⁱ)
            int iPow0 = iPowInit + Poly2Pow[symbol];
            MOD_NN(iPow0);

            for (int i = 1; i <= 2*tt; i++) {
                MOD_NN(iPow0);
                syndromes[i] ^= Pow2Poly[iPow0];
                iPow0 += j;
            }
        }
        iPowInit += b0;
        MOD_NN(iPowInit);
    }

    /* Convert syndromes to power form
     * In GF(2^8), addition is XOR in polynomial form
     * but multiplication is addition of exponents in power form
     */
    for (int i = 1; i <= 2*tt; i++) {
        if (syndromes[i] != 0) {
            syn_error = 1;
        }
        syndromes[i] = Poly2Pow[syndromes[i]];
    }
#endif

    RsVerification::verify_syndromes(syndromes, tt);
    RsVerification::print_syndromes("Initial", syndromes, tt);

    if (syn_error) {
        // Since we have non-zero syndromes, attempt error correction

        /* Step 2: Berlekamp-Massey Algorithm
         *
         * Theory: The error locator polynomial Λ(x) has the property that:
         * Λ(x) = Π(1 - xXi) where Xi are the error locations (Λ(Xi) = 0)
         *
         * The algorithm iteratively constructs Λ(x) by:
         * 1. Using current polynomial to predict next syndrome
         * 2. Computing discrepancy between prediction and actual syndrome
         * 3. Updating polynomial based on discrepancy
         *
         * Key properties:
         * - deg(Λ) equals number of errors
         * - Polynomial coefficients are in GF(2^8)
         * - Updates maintain the syndrome equations
         */
        GF elp[2*MAX_TT+2][2*MAX_TT];  // Error locator polynomials
        GF d[2*MAX_TT+2];              // Discrepancy values
        int l[2*MAX_TT+2];             // Degree of each ELP
        int u_lu[2*MAX_TT+2];          // u-l(u) values
        int u = 0;                      // Step number

        // Initialize for first two steps
        d[0] = 0;           // power form
        d[1] = syndromes[1];// power form
        elp[0][0] = 0;      // power form
        elp[1][0] = 1;      // polynomial form
        for (int i = 1; i < nn-kk; i++) {
            elp[0][i] = GF_INFINITY;   // power form
            elp[1][i] = 0;             // polynomial form
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;

        do {
            u++;
            if (d[u] == GF_INFINITY) {
                // If discrepancy is zero, copy previous polynomial
                l[u+1] = l[u];
                for (int i = 0; i <= l[u]; i++) {
                    elp[u+1][i] = elp[u][i];
                    elp[u][i] = Poly2Pow[elp[u][i]];  // Convert old ELP to power form
                }
            } else {
                // Find earlier error polynomial with max u-l(u)
                int q = u-1;
                while ((d[q] == GF_INFINITY) && (q > 0)) {
                    q--;
                }
                if (q > 0) {
                    int j = q;
                    do {
                        j--;
                        if ((d[j] != GF_INFINITY) && (u_lu[q] < u_lu[j])) {
                            q = j;
                        }
                    } while (j > 0);
                }

                // Update polynomial degree
                if (l[u] > l[q] + u - q) {
                    l[u+1] = l[u];
                } else {
                    l[u+1] = l[q] + u - q;
                }

                // Form new error locator polynomial
                for (int i = 0; i < nn-kk; i++) {
                    elp[u+1][i] = 0;
                }

                // Add scaled version of selected earlier polynomial
                for (int i = 0; i <= l[q]; i++) {
                    if (elp[q][i] != GF_INFINITY) {
                        int iMod = d[u] - d[q] + elp[q][i];
                        MOD_NN(iMod);
                        elp[u+1][i+u-q] = Pow2Poly[iMod];
                    }
                }

                // Add current polynomial
                for (int i = 0; i <= l[u]; i++) {
                    elp[u+1][i] ^= elp[u][i];
                    elp[u][i] = Poly2Pow[elp[u][i]];  // Convert old ELP to power form
                }
            }

            u_lu[u+1] = u - l[u+1];

            // Calculate next discrepancy if not at end
            if (u < nn-kk) {
                // d[u+1] = S[u+1] + Σ(j=1 to l[u+1]) elp[u+1][j] * S[u+1-j]
                d[u+1] = Pow2Poly[syndromes[u+1]];
                for (int i = 1; i <= l[u+1]; i++) {
                    if ((syndromes[u+1-i] != GF_INFINITY) && (elp[u+1][i] != 0)) {
                        int iMod = syndromes[u+1-i] + Poly2Pow[elp[u+1][i]];
                        MOD_NN(iMod);
                        d[u+1] ^= Pow2Poly[iMod];
                    }
                }
                d[u+1] = Poly2Pow[d[u+1]];  // Convert to power form
            }

            RsVerification::print_berlekamp_step(u, d[u], l[u], elp[u], l[u]);

        } while ((u < nn-kk) && (l[u+1] <= tt));

        u++;
        RsVerification::verify_lambda(elp[u], l[u]);

        // We have a valid error locator polynomial if degree <= tt
        if (l[u] <= tt) {
            // Convert final error locator polynomial to power form
            for (int i = 0; i <= l[u]; i++) {
                elp[u][i] = Poly2Pow[elp[u][i]];
            }

            /* Step 3: Chien Search
             *
             * Theory: To find error locations, we need roots of Λ(x)
             * For each j, evaluate Λ(α^-j):
             * - If Λ(α^-j) = 0, then j is an error location
             * - We expect deg(Λ) roots if Λ(x) is correct
             *
             * Optimization: Instead of recomputing powers, use a register
             * that's updated by multiplying each coefficient by appropriate power
             */

            // Initialize evaluation register with coefficients
            for (int i = 1; i <= l[u]; i++) {
                reg[i] = elp[u][i];
            }

            count = 0;  // Count of roots found

            // Try each possible position
            for (int i = 1; i <= nn; i++) {
                GF q = 1;  // Accumulator for polynomial evaluation

                // For each term in the polynomial
                for (int j = 1; j <= l[u]; j++) {
                    if (reg[j] != GF_INFINITY) {
                        // Multiply coefficient by next power of x
                        int iMod = reg[j] + j;
                        MOD_NN(iMod);
                        q ^= Pow2Poly[iMod];

                        // Update register for next position
                        reg[j] = iMod;
                    }
                }

                if (!q) {  // We found a root
                    root[count] = i;         // Save root in power form
                    loc[count] = nn - i;     // Save error position
                    count++;
                }
            }

            // Verify we found the expected number of roots
            if (count == l[u]) {
                // Found correct number of roots, continue with error value computation
                /* Step 4: Forney Algorithm
                 *
                 * Theory: Error evaluator polynomial Ω(x) is computed from:
                 * Ω(x) = S(x)Λ(x) mod x^(2t)
                 * where S(x) is the syndrome polynomial
                 *
                 * For each error location Xk = α^-k:
                 * The error value ek is computed as:
                 * ek = -(Xk^(1-b0) * Ω(Xk^-1)) / Λ'(Xk^-1)
                 * where Λ'(x) is the formal derivative of Λ(x)
                 */

                // First compute error evaluator polynomial Ω(x)
                GF z[MAX_TT+2];  // Error evaluator polynomial

                // For each position (skipping z[0] which is always 1)
                for (int i = 1; i <= l[u]; i++) {
                    // Start with the syndrome term
                    GF zi = Pow2Poly[syndromes[i]] ^ Pow2Poly[elp[u][i]];

                    // Add correction terms
                    for (int j = 1; j < i; j++) {
                        if ((syndromes[j] != GF_INFINITY) && (elp[u][i-j] != GF_INFINITY)) {
                            int iMod = elp[u][i-j] + syndromes[j];
                            MOD_NN(iMod);
                            zi ^= Pow2Poly[iMod];
                        }
                    }
                    z[i] = Poly2Pow[zi];  // Store in power form
                }

                // Clear error array
                for (int i = 0; i < nn; i++) {
                    err[i] = 0;
                }

                // For each error location found
                for (int i = 0; i < l[u]; i++) {
                    // Position in codeword
                    int position = loc[i];

                    // Evaluate Ω(x) at α^(-j)
                    int jPow = root[i];
                    err[position] = 1;  // Start with 1 for z[0]

                    for (int j = 1; j <= l[u]; j++) {
                        if (z[j] != GF_INFINITY) {
                            int iMod = z[j] + jPow;
                            MOD_NN(iMod);
                            err[position] ^= Pow2Poly[iMod];
                        }
                        jPow += root[i];  // Next power
                        MOD_NN(jPow);
                    }

                    // Apply correction formula
                    if (err[position] != 0) {
                        // Multiply by α^(1-b0)
                        err[position] = Pow2Poly[(Poly2Pow[err[position]] +
                                                root[i]*(b0+nn-1)) % nn];

                        if (err[position] != 0) {
                            err[position] = Poly2Pow[err[position]];

                            // Compute formal derivative denominator
                            int q = 0;
                            for (int j = 0; j < l[u]; j++) {
                                if (j != i) {
                                    int iMod = loc[j] + root[i];
                                    MOD_NN(iMod);
                                    q += Poly2Pow[1 ^ Pow2Poly[iMod]];
                                    MOD_NN(q);
                                }
                            }

                            // Final error value computation
                            int iMod = err[position] - q;
                            MOD_NN(iMod);
                            err[position] = Pow2Poly[iMod];

                            // Apply correction to received word
                            GF original = recd[position];
                            recd[position] ^= err[position];
                            RsVerification::print_error_location(position, err[position],
                                                        original, recd[position]);
                        }
                    }
                }

                return count;  // Return number of errors corrected
            } else {
                // Number of roots doesn't match polynomial degree
                return RS_ERROR_CHIEN_SEARCH;
            }
        } else {
            // Error locator polynomial degree too high
            return RS_ERROR_LAMDA_ERROR;
        }
    }
    else
    {
        /* no non-zero syndromes => no errors: output received codeword */
        return 0;  // No errors found
    }
}

int RS_ENCODER_REF::RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras) {
    RsVerification::verify_received_word(recd, nn);

    /* Reed-Solomon Decoding with Erasures
     *
     * Mathematical Background:
     * 1. An erasure is an error whose position is known in advance
     * 2. Each erasure counts as one condition on the error pattern
     * 3. Total corrections = errors + erasures ≤ 2t
     * 4. Known erasure positions let us modify initial conditions for B-M algorithm
     *
     * Key differences from standard decoding:
     * 1. First compute erasure locator polynomial Γ(x) from known positions
     * 2. Initialize error locator polynomial Λ(x) = Γ(x)
     * 3. Continue with modified Berlekamp-Massey to find remaining errors
     * 4. Final polynomial will locate both errors and erasures
     */

    // Working buffers
    GF syndromes[2*MAX_TT+1];    // Si = r(α^(b0+i))
    GF root[n];                  // Roots of final locator polynomial
    GF loc[n];                   // Error/erasure positions
    GF err[n];                   // Error/erasure values
    GF reg[2*MAX_TT+1];         // Working register for Chien search
    int syn_error = 0;          // Non-zero syndrome flag

    // Additional buffers for erasure handling
    GF phi[2*MAX_TT+1];         // Erasure locator polynomial
    GF tmp_pol[2*MAX_TT+1];     // Temporary polynomial for computations
    GF lambda[2*MAX_TT+1];      // Error+erasure locator polynomial
    GF lambda_pr[2*MAX_TT+1];   // Formal derivative of lambda
    GF omega[2*MAX_TT+1];       // Error evaluator polynomial
    GF b[2*MAX_TT+1];           // Temporary for B-M algorithm
    GF T[2*MAX_TT+1];           // Another temporary for B-M

    /* Step 1: Compute erasure locator polynomial Γ(x)
     * Γ(x) = Π(1 - xX[i]) for each erasure position X[i]
     * where X[i] = α^pos[i]
     */

    memset(tmp_pol, 0, sizeof(tmp_pol));
    memset(phi, 0, sizeof(phi));

    if (no_eras > 0) {
        phi[0] = 1;  // constant term
        phi[1] = Pow2Poly[eras_pos[0]];  // First erasure factor

        // Multiply by (1 - xX[i]) for each erasure
        for (int i = 1; i < no_eras; i++) {
            GF U = eras_pos[i];
            for (int j = 1; j < i+2; j++) {
                if (phi[j-1] == 0) {
                    tmp_pol[j] = 0;
                } else {
                    int iMod = U + Poly2Pow[phi[j-1]];
                    MOD_NN(iMod);
                    tmp_pol[j] = Pow2Poly[iMod];
                }
            }
            for (int j = 1; j < i+2; j++) {
                phi[j] = phi[j] ^ tmp_pol[j];
            }
        }

        // Convert phi to power form
        for (int i = 0; i < nn-kk+1; i++) {
            phi[i] = Poly2Pow[phi[i]];
        }
    }

    /* Step 2: Calculate Syndromes
     * Same as standard decoding, but will be used with
     * the initialized Λ(x) = Γ(x)
     */

    // Initialize syndromes to zero
    for (int i = 0; i <= 2*tt; i++) {
        syndromes[i] = 0;
    }

    // For each position in received word
    int iPowInit = 0;
    for (int j = 0; j < nn; j++) {
        GF RECD = recd[j];
        if (RECD != 0) {
            int iPow0 = iPowInit + Poly2Pow[RECD];
            MOD_NN(iPow0);

            for (int i = 1; i <= 2*tt; i++) {
                MOD_NN(iPow0);
                syndromes[i] ^= Pow2Poly[iPow0];
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

    RsVerification::verify_syndromes(syndromes, tt);
    RsVerification::print_syndromes("Initial", syndromes, tt);

    if (syn_error) {
        /* Step 3: Modified Berlekamp-Massey Algorithm
         * Initialize with Λ(x) = Γ(x) if there were erasures
         * Otherwise start with Λ(x) = 1
         */
        int r = no_eras;
        int deg_phi = no_eras;
        int L = no_eras;

        if (no_eras > 0) {
            // Initialize lambda(x) and b(x) to phi(x)
            for (int i = 0; i < deg_phi+1; i++) {
                lambda[i] = Pow2Poly[phi[i]];
            }
            for (int i = deg_phi+1; i < 2*tt+1; i++) {
                lambda[i] = 0;
            }
            for (int i = 0; i < 2*tt+1; i++) {
                b[i] = lambda[i];
            }
        } else {
            lambda[0] = 1;
            for (int i = 1; i < 2*tt+1; i++) {
            lambda[i] = 0;
            }
            for (int i = 0; i < 2*tt+1; i++) {
                b[i] = lambda[i];
            }
        }

        while (++r <= 2*tt) {  // r is the step number
            // Compute discrepancy at the r-th step in poly-form
            GF discr_r = 0;
            for (int i = 0; i < r; i++) {
                if ((lambda[i] != 0) && (syndromes[r-i] != GF_INFINITY)) {
                    int iMod = Poly2Pow[lambda[i]] + syndromes[r-i];
                    MOD_NN(iMod);
                    discr_r ^= Pow2Poly[iMod];
                }
            }

            if (discr_r == 0) {
                // Shift b(x) by x
                for (int i = 2*tt; i > 0; i--) {
                    b[i] = b[i-1];
                }
                b[0] = 0;
            } else {
                // T(x) = lambda(x) - discr_r*x*b(x)
                T[0] = lambda[0];
                for (int i = 1; i < 2*tt+1; i++) {
                    if (b[i-1] == 0) {
                        T[i] = lambda[i];
                    } else {
                        int iMod = Poly2Pow[discr_r] + Poly2Pow[b[i-1]];
                        MOD_NN(iMod);
                        T[i] = lambda[i] ^ Pow2Poly[iMod];
                    }
                }

                if (2*L <= r+no_eras-1) {
                    L = r+no_eras-L;
                    // b(x) = inv(discr_r) * lambda(x)
                    for (int i = 0; i < 2*tt+1; i++) {
                        if (lambda[i] == 0) {
                            b[i] = 0;
                        } else {
                            int iMod = Poly2Pow[lambda[i]] - Poly2Pow[discr_r];
                            MOD_NN(iMod);
                            b[i] = Pow2Poly[iMod];
                        }
                    }
                    for (int i = 0; i < 2*tt+1; i++) {
                        lambda[i] = T[i];
                    }
                } else {
                    for (int i = 0; i < 2*tt+1; i++) {
                        lambda[i] = T[i];
                    }
                    // b(x) = x*b(x)
                    for (int i = 2*tt; i > 0; i--) {
                        b[i] = b[i-1];
                    }
                    b[0] = 0;
                }
            }
        }

        /* Step 4: Find roots (Chien search) and compute error values (Forney)
         * Same as standard decoding, but handles both errors and erasures
         */
        // Put lambda(x) into power form
        for (int i = 0; i < 2*tt+1; i++) {
            lambda[i] = Poly2Pow[lambda[i]];
        }

        // Compute deg(lambda(x))
        int deg_lambda = 2*tt;
        while ((lambda[deg_lambda] == GF_INFINITY) && (deg_lambda > 0)) {
            --deg_lambda;
        }

        RsVerification::verify_lambda(lambda, deg_lambda);

        if (deg_lambda <= 2*tt) {
            // Find roots of error+erasure locator polynomial using Chien Search
            for (int i = 1; i < 2*tt+1; i++) {
                reg[i] = lambda[i];
            }

            int count = 0;  // Number of roots found
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
                if (!q) {  // Store root and error location number
                    root[count] = i;
                    loc[count] = nn-i;
                    count++;
                }
            }

            if (deg_lambda == count) {  // Number of roots equals degree of lambda
                // Compute error+erasure evaluator polynomial
                for (int i = 0; i < 2*tt; i++) {
                    omega[i] = 0;
                    for (int j = 0; j < deg_lambda+1 && j < i+1; j++) {
                        if ((syndromes[i+1-j] != GF_INFINITY) && (lambda[j] != GF_INFINITY)) {
                            int iMod = syndromes[i+1-j] + lambda[j];
                            MOD_NN(iMod);
                            omega[i] ^= Pow2Poly[iMod];
                        }
                    }
                }
                omega[2*tt] = 0;

                // Compute formal derivative of lambda
                for (int i = 0; i < 2*tt; i += 2) {
                    lambda_pr[i+1] = 0;
                    lambda_pr[i] = Pow2Poly[lambda[i+1]];
                }
                lambda_pr[2*tt] = 0;

                // Find degree of omega
                int deg_omega = 2*tt;
                while ((omega[deg_omega] == 0) && (deg_omega > 0)) {
                    --deg_omega;
                }

                // Compute error values
                for (int j = 0; j < count; j++) {
                    GF pres_root = root[j];
                    GF pres_loc = loc[j];

                    // Calculate numerator
                    GF num1 = 0;
                    int iPow0 = 0;
                    for (int i = 0; i < deg_omega+1; i++) {
                        if (omega[i] != 0) {
                            int iMod = Poly2Pow[omega[i]] + iPow0;
                            MOD_NN(iMod);
                            num1 ^= Pow2Poly[iMod];
                        }
                        iPow0 += pres_root;
                        MOD_NN(iPow0);
                    }

                    // Additional term for erasure decoding
                    int iMod = b0-1;
                    MOD_NN(iMod);
                    iMod = (iMod*pres_root)%nn;
                    GF num2 = Pow2Poly[iMod];

                    // Calculate denominator (formal derivative evaluation)
                    GF den = 0;
                    iPow0 = 0;
                    for (int i = 0; i < deg_lambda+1; i++) {
                        if (lambda_pr[i] != 0) {
                            iMod = Poly2Pow[lambda_pr[i]] + iPow0;
                            MOD_NN(iMod);
                            den ^= Pow2Poly[iMod];
                        }
                        iPow0 += pres_root;
                        MOD_NN(iPow0);
                    }

                    err[pres_loc] = 0;
                    if (num1 != 0 && den != 0) {
                        iMod = Poly2Pow[num1] + Poly2Pow[num2] - Poly2Pow[den];
                        MOD_NN(iMod);
                        err[pres_loc] = Pow2Poly[iMod];
                    }

                    // Apply correction
                    recd[pres_loc] ^= err[pres_loc];
                }

                return count;  // Return number of corrections
            } else {
                return RS_ERROR_CHIEN_SEARCH;  // Number of roots != degree of lambda
            }
        } else {
            return RS_ERROR_LAMDA_ERROR;  // Degree of lambda too high
        }
    } else {
        return 0;  // No errors found
    }
}
bool RS_ENCODER_REF::bInitialized = false;
