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

#if 0
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
#endif

int RS_ENCODER_REF::RSDecode(GF recd[nn]) {
    return RSDecodeErasures(recd, nullptr, 0);
}

#define A0	(nn)

// RSDecodeErasures: Port of Phil Karn's Berlekamp-Massey implementation
// Brain-dead port for erasure decoding with no erasures initially
int RS_ENCODER_REF::RSDecodeErasures(GF data[nn], int eras_pos[], int no_eras) {
    RsVerification::verify_received_word(data, nn);

    int kk = nn - 2*tt;
    int b0 = (kk+1)/2;
    int deg_lambda, el, deg_omega;
    int i, j, r;
    GF u,q,tmp,num1,num2,den,discr_r;
    GF recd[nn];
    GF lambda[2 * MAX_TT + 1], s[2 * MAX_TT + 1];	/* Err+Eras Locator poly
                         * and syndrome poly */
    GF b[2 * MAX_TT + 1], t[2 * MAX_TT + 1], omega[2 * MAX_TT + 1];
    GF root[2 * MAX_TT], reg[2 * MAX_TT + 1], loc[2 * MAX_TT];
    int syn_error, count;

    /* data[] is in polynomial form, copy and convert to index form */
    for (i = nn-1; i >= 0; i--){
        recd[i] = Poly2Pow[data[i]];
    }
    /* first form the syndromes; i.e., evaluate recd(x) at roots of g(x)
     * namely @**(B0+i), i = 0, ... ,(NN-KK-1)
     */
    syn_error = 0;
    s[0] = 0;
    for (i = 1; i <= 2*tt; i++) {
        tmp = 0;
        for (j = 0; j < nn; j++)
            if (recd[j] != A0)	/* recd[j] in index form */
                tmp ^= Pow2Poly[mod_nn(recd[j] + (b0+i-1)*j)];
        syn_error |= tmp;	/* set flag if non-zero syndrome =>
                     * error */
        /* store syndrome in index form  */
        s[i] = Poly2Pow[tmp];
    }

    RsVerification::verify_syndromes(s, tt);
    RsVerification::print_syndromes("Initial", s, tt);

    if (!syn_error) {
        /*
         * if syndrome is zero, data[] is a codeword and there are no
         * errors to correct. So return data[] unmodified
         */
        return 0;
    }
    memset(&lambda[1], 0, 2*tt*sizeof(lambda[0]));
    lambda[0] = 1;
    if (no_eras > 0) {
        /* Init lambda to be the erasure locator polynomial */
        lambda[1] = Pow2Poly[eras_pos[0]];
        for (i = 1; i < no_eras; i++) {
            u = eras_pos[i];
            for (j = i+1; j > 0; j--) {
                tmp = Poly2Pow[lambda[j - 1]];
                if(tmp != A0)
                    lambda[j] ^= Pow2Poly[mod_nn(u + tmp)];
            }
        }
    }
    for(i=0;i<2*tt+1;i++)
        b[i] = Poly2Pow[lambda[i]];

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    r = no_eras;
    el = no_eras;
    while (++r <= 2*tt) {	/* r is the step number */
        /* Compute discrepancy at the r-th step in poly-form */
        discr_r = 0;
        for (i = 0; i < r; i++){
            if ((lambda[i] != 0) && (s[r - i] != A0)) {
                discr_r ^= Pow2Poly[mod_nn(Poly2Pow[lambda[i]] + s[r - i])];
            }
        }
        RsVerification::print_berlekamp_step(r, discr_r, r, std::vector<GF>(lambda, lambda + r + 1), r);

        discr_r = Poly2Pow[discr_r];	/* Index form */
        if (discr_r == A0) {
            /* 2 lines below: B(x) <-- x*B(x) */
            memmove(&b[1], b, 2*tt*sizeof(b[0]));
            b[0] = A0;
        } else {
            /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
            t[0] = lambda[0];
            for (i = 0 ; i < 2*tt; i++) {
                if(b[i] != A0)
                    t[i+1] = lambda[i+1] ^ Pow2Poly[mod_nn(discr_r + b[i])];
                else
                    t[i+1] = lambda[i+1];
            }
            if (2 * el <= r + no_eras - 1) {
                el = r + no_eras - el;
                /*
                 * 2 lines below: B(x) <-- inv(discr_r) *
                 * lambda(x)
                 */
                for (i = 0; i <= 2*tt; i++)
                    b[i] = (lambda[i] == 0) ? A0 : mod_nn(Poly2Pow[lambda[i]] - discr_r + nn);
            } else {
                /* 2 lines below: B(x) <-- x*B(x) */
                memmove(&b[1], b, 2*tt*sizeof(b[0]));
                b[0] = A0;
            }
            memcpy(lambda, t, (2*tt+1) * sizeof(lambda[0]));
        }
    }

    /* Convert lambda to index form and compute deg(lambda(x)) */
    deg_lambda = 0;
    for(i=0;i<2*tt+1;i++){
        lambda[i] = Poly2Pow[lambda[i]];
        if(lambda[i] != A0)
            deg_lambda = i;
    }

    RsVerification::verify_lambda(lambda, deg_lambda);

    /*
     * Find roots of the error+erasure locator polynomial. By Chien
     * Search
     */
    memcpy(&reg[1], &lambda[1], 2*tt*sizeof(lambda[0]));
    count = 0;		/* Number of roots of lambda(x) */
    for (i = 1; i <= nn; i++) {
        q = 1;
        for (j = deg_lambda; j > 0; j--)
            if (reg[j] != A0) {
                reg[j] = mod_nn(reg[j] + j);
                q ^= Pow2Poly[reg[j]];
            }
        if (!q) {
            /* store root (index-form) and error location number */
            root[count] = i;
            loc[count] = nn - i;
            count++;
        }
    }

    if (deg_lambda != count) {
        /*
         * deg(lambda) unequal to number of roots => uncorrectable
         * error detected
         */
        return -1;
    }
    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**(NN-KK)). in index form. Also find deg(omega).
     */
    deg_omega = 0;
    for (i = 0; i < 2*tt;i++){
        tmp = 0;
        j = (deg_lambda < i) ? deg_lambda : i;
        for(;j >= 0; j--){
            if ((s[i + 1 - j] != A0) && (lambda[j] != A0))
                tmp ^= Pow2Poly[mod_nn(s[i + 1 - j] + lambda[j])];
        }
        if(tmp != 0)
            deg_omega = i;
        omega[i] = Poly2Pow[tmp];
    }
    omega[2*tt] = A0;

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(B0-1) and den = lambda_pr(inv(X(l))) all in poly-form
     */
    for (j = count-1; j >=0; j--) {
        num1 = 0;
        for (i = deg_omega; i >= 0; i--) {
            if (omega[i] != A0)
                num1  ^= Pow2Poly[mod_nn(omega[i] + i * root[j])];
        }
        num2 = Pow2Poly[mod_nn(root[j] * (b0 - 1) + nn)];
        den = 0;

        /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
        for (i = std::min(deg_lambda,2*tt-1) & ~1; i >= 0; i -=2) {
            if(lambda[i+1] != A0)
                den ^= Pow2Poly[mod_nn(lambda[i+1] + i * root[j])];
        }
        if (den == 0) {
            return -1;
        }
        /* Apply error to data */
        if (num1 != 0) {
            data[loc[j]] ^= Pow2Poly[mod_nn(Poly2Pow[num1] + Poly2Pow[num2] + nn - Poly2Pow[den])];
        }
    }
    return count;
}
