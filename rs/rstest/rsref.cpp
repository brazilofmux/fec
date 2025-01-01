#include <string.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "rs.h"
#include "rsref.h"
#include "RsVerification.h"

static inline void MOD_NN(int& x) {
    while (x >= nn) {
        x -= nn;
        x = (x >> 8) + (x & nn);
    }
}

static inline int mod_nn(int x) {
    while (x >= nn) {
        x -= nn;
        x = (x >> 8) + (x & nn);
    }
    return x;
}

RS_ENCODER_REF::RS_ENCODER_REF(int CorrectableErrors) : tt(CorrectableErrors), kk(nn - 2 * CorrectableErrors) {
    b0 = (kk + 1) / 2;  // Symmetric generator polynomial
    RSGenPoly();
    RsVerification::verify_generator(generator_poly, tt);
}

void RS_ENCODER_REF::RSGenPoly() {
    generator_poly.assign(2 * tt + 1, 0);
    generator_poly[0] = Pow2Poly[b0];
    generator_poly[1] = 1;

    for (int i = 2; i <= 2 * tt; ++i) {
        generator_poly[i] = 1;
        for (int j = i - 1; j > 0; --j) {
            if (generator_poly[j] != 0) {
                int index = Poly2Pow[generator_poly[j]] + b0 + i - 1;
                generator_poly[j] = generator_poly[j - 1] ^ Pow2Poly[mod_nn(index)];
            }
            else {
                generator_poly[j] = generator_poly[j - 1];
            }
        }
        generator_poly[0] = Pow2Poly[mod_nn(Poly2Pow[generator_poly[0]] + b0 + i - 1)];
    }
    for (auto& coeff : generator_poly) {
        coeff = Poly2Pow[coeff];
    }
}

void RS_ENCODER_REF::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    std::fill(bb, bb + 2 * tt, 0);

    for (int i = 0; i < kk; ++i) {
        GF feedback = Poly2Pow[data[i] ^ bb[0]];  // Adjust for descending order
        if (feedback != GF_INFINITY) {
            for (int j = 0; j < 2 * tt - 1; ++j) {
                if (generator_poly[j + 1] != GF_INFINITY) {
                    bb[j] = bb[j + 1] ^ Pow2Poly[mod_nn(generator_poly[j + 1] + feedback)];
                } else {
                    bb[j] = bb[j + 1];
                }
            }
            bb[2 * tt - 1] = Pow2Poly[mod_nn(generator_poly[0] + feedback)];
        } else {
            for (int j = 0; j < 2 * tt - 1; ++j) {
                bb[j] = bb[j + 1];
            }
            bb[2 * tt - 1] = 0;
        }
    }
}

void RS_ENCODER_REF::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt + 1, 0);
    int iPowInit = 0;

    for (int j = 0; j < nn; j++) {
        if (recd[j] != 0) {
            int iPow0 = iPowInit + Poly2Pow[recd[j]];
            iPow0 = mod_nn(iPow0);

            for (int i = 1; i <= 2 * tt; i++) {
                iPow0 = mod_nn(iPow0);
                syndromes[i] ^= Pow2Poly[iPow0];
                iPow0 += j;
            }
        }
        iPowInit += b0;
        iPowInit = mod_nn(iPowInit);
    }

    for (int i = 1; i <= 2 * tt; i++) {
        syndromes[i] = Poly2Pow[syndromes[i]];
    }
}

void RS_ENCODER_REF::berlekamp_massey(const std::vector<GF>& syndromes, int tt,
                                      std::vector<GF>& lambda, int& lambda_deg) {
    const int max_deg = 2 * tt;

    // Initialize tables
    std::vector<GF> b(max_deg + 1, 0); // Previous lambda(x)
    lambda.assign(max_deg + 1, 0);     // Current lambda(x)
    lambda[0] = 1;                     // Polynomial form

    int L = 0;  // Degree of lambda(x)
    int m = 1;  // Distance since last update
    GF discr = 0; // Discrepancy (polynomial form)

    for (int r = 1; r <= max_deg; ++r) {
        // Compute discrepancy: discr = S[r] + sum(lambda[i] * S[r-i])
        discr = syndromes[r];
        for (int i = 1; i <= L; ++i) {
            if (lambda[i] != 0 && syndromes[r - i] != GF_INFINITY) {
                discr ^= Pow2Poly[mod_nn(Poly2Pow[lambda[i]] + syndromes[r - i])];
            }
        }

        RsVerification::print_berlekamp_step(r, discr, L, lambda, L);

        if (discr == 0) {
            ++m; // No update needed
        } else {
            std::vector<GF> temp = lambda; // Save current lambda(x)

            // Update lambda(x): lambda = lambda - discr * x^m * b
            for (int i = 0; i <= max_deg; ++i) {
                if (b[i] != 0) {
                    int index = mod_nn(Poly2Pow[discr] + Poly2Pow[b[i]]);
                    if (i + m <= max_deg) {
                        lambda[i + m] ^= Pow2Poly[index];
                    }
                }
            }

            if (2 * L < r) {
                L = r - L;
                b = temp; // Update b to old lambda(x)
                m = 1;
            } else {
                ++m;
            }
        }
    }

    lambda_deg = L;
}

int RS_ENCODER_REF::RSDecode(GF recd[nn]) {
    RsVerification::verify_received_word(recd, nn);

    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    bool has_error = false;
    for (int i = 1; i <= 2 * tt; i++) {
        if (syndromes[i] != 0) {
            has_error = true;
        }
    }

    RsVerification::verify_syndromes(syndromes, tt);
    RsVerification::print_syndromes("Initial", syndromes, tt);

    if (!has_error) {
        return 0; // No errors detected
    }

    std::vector<GF> lambda;
    int lambda_deg = 0;
    berlekamp_massey(syndromes, tt, lambda, lambda_deg);
    RsVerification::verify_lambda(lambda, lambda_deg);

    std::vector<GF> roots(lambda_deg);
    int num_errors = 0;
    for (int i = 0; i < nn; ++i) {
        GF eval = 1;
        for (int j = 1; j <= lambda_deg; ++j) {
            if (lambda[j] != GF_INFINITY) {
                eval ^= Pow2Poly[mod_nn(Poly2Pow[lambda[j]] + i * j)];
            }
        }
        if (eval == 0) {
            roots[num_errors++] = i;
        }
    }

    if (num_errors != lambda_deg) {
        return -2; // Error count mismatch
    }

    for (int i = 0; i < num_errors; ++i) {
        int position = nn - 1 - roots[i];
        GF error_val = 0;

        for (int j = 1; j <= lambda_deg; ++j) {
            if (syndromes[j] != GF_INFINITY) {
                error_val ^= Pow2Poly[mod_nn(Poly2Pow[syndromes[j]] + j * roots[i])];
            }
        }

        recd[position] ^= error_val;
    }

    return num_errors;
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
