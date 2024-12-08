#include "rs.h"
#include <string.h>

// Tables for converting between power form and polynomial form
static GF Pow2Poly[n], Poly2Pow[n];
#define PrimitivePolynomial  0x1D    // [x^8] + x^4 + x^3 + x^2 + 1

// Simple modulo function for GF operations
static inline void MOD_NN(int& x) {
    x = x % nn;
    if (x < 0) x += nn;  // Handle negative case just to be thorough
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
    }
}

RS_ENCODER::RS_ENCODER(int CorrectableErrors) {
    tt = CorrectableErrors;
    kk = nn-2*tt;

    // Choose generator polynomial starting power
    b0 = (kk+1)/2;
    RSGenPoly();
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
            if (gg[j-1] != 0) {
                iMod = Poly2Pow[gg[j-1]]+b0+i-1;
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
    GF syndromes[2*MAX_TT+1];    // Syndromes
    GF lambda[2*MAX_TT+1];       // Error locator polynomial
    GF omega[2*MAX_TT+1];        // Error evaluator polynomial
    GF root[MAX_TT];             // Error locations
    GF loc[MAX_TT];              // Error position indices
    GF reg[2*MAX_TT+1];          // Scratch space for Chien search
    int deg_lambda;              // Degree of error locator polynomial
    int syn_error = 0;
    int count = 0;               // Number of errors found

    // First form the syndromes: evaluate recd(x) at roots of g(x)
    memset(syndromes, 0, sizeof(syndromes));

    // For each position in received word
    for (int j = 0; j < nn; j++) {
        if (recd[j] != 0) {
            // For each syndrome
            int iPow0 = Poly2Pow[recd[j]];
            for (int i = 1; i <= 2*tt; i++) {
                int iMod = iPow0 + j*i;
                MOD_NN(iMod);
                syndromes[i] ^= Pow2Poly[iMod];
            }
        }
    }

    // Convert syndromes to power form and check for errors
    for (int i = 1; i <= 2*tt; i++) {
        if (syndromes[i] != 0) {
            syn_error = 1;
        }
        syndromes[i] = Poly2Pow[syndromes[i]];
    }

    if (syn_error) {  // If errors, attempt to correct
        // Initialize for Berlekamp-Massey algorithm
        lambda[0] = 1;
        lambda[1] = Pow2Poly[syndromes[1]];
        for (int i = 2; i <= 2*tt; i++) {
            lambda[i] = 0;
        }

        int L = 0;
        int k = 0;

        // Berlekamp-Massey algorithm
        for (int r = 0; r < 2*tt; r++) {
            // Compute discrepancy
            GF d = syndromes[r+1];
            for (int i = 1; i <= L; i++) {
                if (lambda[i] != 0 && syndromes[r-i+1] != GF_INFINITY) {
                    int iMod = Poly2Pow[lambda[i]] + syndromes[r-i+1];
                    MOD_NN(iMod);
                    d ^= Pow2Poly[iMod];
                }
            }

            if (d != 0) {
                GF T[2*MAX_TT+1];
                memcpy(T, lambda, sizeof(T));

                for (int i = 0; i + k < 2*tt; i++) {
                    if (lambda[i] != 0) {
                        int iMod = Poly2Pow[lambda[i]] + Poly2Pow[d];
                        MOD_NN(iMod);
                        lambda[i+k] ^= Pow2Poly[iMod];
                    }
                }

                if (2*L <= r) {
                    L = r + 1 - L;
                    k = 1;
                    memcpy(lambda, T, sizeof(lambda));
                } else {
                    k++;
                }
            } else {
                k++;
            }
        }

        // Convert lambda to power form
        for (int i = 0; i <= 2*tt; i++) {
            lambda[i] = Poly2Pow[lambda[i]];
        }

        // Find degree of lambda
        deg_lambda = 2*tt;
        while (deg_lambda > 0 && lambda[deg_lambda] == GF_INFINITY) {
            deg_lambda--;
        }

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
                        recd[pos] ^= err;  // Correct the error
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

int RS_ENCODER::RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras)
{
    GF syndromes[2*MAX_TT+1];    // Syndromes
    GF lambda[2*MAX_TT+1];       // Error locator polynomial
    GF reg[2*MAX_TT+1];          // Scratch for Chien search
    GF root[MAX_TT];             // Error locations
    GF loc[MAX_TT];              // Error position indices
    GF omega[2*MAX_TT+1];        // Error evaluator polynomial
    
    // Initialize erasure locator polynomial
    GF phi[2*MAX_TT+1], tmp[2*MAX_TT+1];
    memset(phi, 0, sizeof(phi));
    memset(tmp, 0, sizeof(tmp));
    
    phi[0] = 1;  // Initialize phi(x)
    
    // Build erasure locator polynomial
    if (no_eras > 0) {
        phi[1] = Pow2Poly[eras_pos[0]];
        for (int i = 1; i < no_eras; i++) {
            // Multiply phi(x) by (1 + x*alpha^eras_pos[i])
            int U = eras_pos[i];
            for (int j = i + 1; j > 0; j--) {
                if (phi[j-1] != 0) {
                    int iMod = Poly2Pow[phi[j-1]] + U;
                    iMod %= nn;
                    if (iMod < 0) iMod += nn;
                    tmp[j] = Pow2Poly[iMod];
                } else {
                    tmp[j] = 0;
                }
            }
            for (int j = 1; j < i + 2; j++) {
                phi[j] ^= tmp[j];
            }
        }
        
        // Convert phi to power form
        for (int i = 0; i <= no_eras; i++) {
            phi[i] = Poly2Pow[phi[i]];
        }
    }
    
    // Form syndromes
    memset(syndromes, 0, sizeof(syndromes));
    
    for (int j = 0; j < nn; j++) {
        if (recd[j] != 0) {
            int iPow0 = Poly2Pow[recd[j]];
            for (int i = 1; i <= 2*tt; i++) {
                int iMod = iPow0 + j*i;
                iMod %= nn;
                if (iMod < 0) iMod += nn;
                syndromes[i] ^= Pow2Poly[iMod];
            }
        }
    }
    
    // Convert syndromes to power form and check for errors
    int syn_error = 0;
    for (int i = 1; i <= 2*tt; i++) {
        if (syndromes[i] != 0) {
            syn_error = 1;
        }
        syndromes[i] = Poly2Pow[syndromes[i]];
    }
    
    if (!syn_error) {
        return 0;  // No errors
    }
    
    // Initialize Berlekamp-Massey
    // Start with erasure locator polynomial if we have erasures
    if (no_eras > 0) {
        for (int i = 0; i <= no_eras; i++) {
            lambda[i] = Pow2Poly[phi[i]];
        }
        for (int i = no_eras + 1; i <= 2*tt; i++) {
            lambda[i] = 0;
        }
    } else {
        lambda[0] = 1;
        for (int i = 1; i <= 2*tt; i++) {
            lambda[i] = 0;
        }
    }
    
    // Berlekamp-Massey algorithm
    int L = no_eras;
    int k = 0;
    
    for (int r = no_eras + 1; r <= 2*tt; r++) {
        // Compute discrepancy
        GF d = syndromes[r];
        for (int i = 1; i <= L; i++) {
            if (lambda[i] != 0 && syndromes[r-i] != GF_INFINITY) {
                int iMod = Poly2Pow[lambda[i]] + syndromes[r-i];
                iMod %= nn;
                if (iMod < 0) iMod += nn;
                d ^= Pow2Poly[iMod];
            }
        }
        d = Poly2Pow[d];
        
        if (d != GF_INFINITY) {
            // Save old lambda and compute new one
            GF T[2*MAX_TT+1];
            memcpy(T, lambda, sizeof(T));
            
            for (int i = 0; i + k < 2*tt; i++) {
                if (lambda[i] != 0) {
                    int iMod = Poly2Pow[lambda[i]] + Poly2Pow[d];
                    iMod %= nn;
                    if (iMod < 0) iMod += nn;
                    lambda[i+k] ^= Pow2Poly[iMod];
                }
            }
            
            if (2*L <= r + no_eras - 1) {
                L = r + no_eras - L;
                k = 1;
                memcpy(lambda, T, sizeof(lambda));
            } else {
                k++;
            }
        } else {
            k++;
        }
    }
    
    // Convert lambda to power form
    for (int i = 0; i <= 2*tt; i++) {
        lambda[i] = Poly2Pow[lambda[i]];
    }
    
    // Find degree of error locator polynomial
    int deg_lambda = 2*tt;
    while (deg_lambda > 0 && lambda[deg_lambda] == GF_INFINITY) {
        deg_lambda--;
    }
    
    if (deg_lambda <= 2*tt) {
        // Find roots of the error locator polynomial
        memcpy(reg, lambda, sizeof(reg));
        int count = 0;
        
        for (int i = 1; i <= nn; i++) {
            GF q = 1;
            for (int j = 1; j <= deg_lambda; j++) {
                if (reg[j] != GF_INFINITY) {
                    int iMod = reg[j] + j;
                    iMod %= nn;
                    if (iMod < 0) iMod += nn;
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
            // Form the error evaluator polynomial omega(x)
            for (int i = 0; i < 2*tt; i++) {
                GF tmp = 0;
                for (int j = 0; j <= deg_lambda && j <= i; j++) {
                    if (syndromes[i+1-j] != GF_INFINITY && lambda[j] != GF_INFINITY) {
                        int iMod = syndromes[i+1-j] + lambda[j];
                        iMod %= nn;
                        if (iMod < 0) iMod += nn;
                        tmp ^= Pow2Poly[iMod];
                    }
                }
                omega[i] = Poly2Pow[tmp];
            }
            
            // Compute error values and correct
            for (int i = 0; i < count; i++) {
                int pos = loc[i];
                
                // Compute numerator
                GF num = 0;
                for (int j = 0; j < deg_lambda; j++) {
                    if (omega[j] != GF_INFINITY) {
                        int iMod = omega[j] + j * root[i];
                        iMod %= nn;
                        if (iMod < 0) iMod += nn;
                        num ^= Pow2Poly[iMod];
                    }
                }
                
                // Compute denominator
                GF den = 0;
                for (int j = 1; j <= deg_lambda; j++) {
                    if (lambda[j] != GF_INFINITY) {
                        int iMod = lambda[j] + (j-1) * root[i];
                        iMod %= nn;
                        if (iMod < 0) iMod += nn;
                        den ^= Pow2Poly[iMod];
                    }
                }
                
                if (den != 0) {
                    int iMod = Poly2Pow[num] - Poly2Pow[den] + root[i] * (b0-1);
                    iMod %= nn;
                    if (iMod < 0) iMod += nn;
                    recd[pos] ^= Pow2Poly[iMod];
                }
            }
            return count;
        }
        return -2;  // Number of roots != degree of lambda
    }
    return -3;  // deg_lambda > tt
}
