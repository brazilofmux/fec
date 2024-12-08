#include "rs.h"
#include <string.h>

// Tables for converting between power form and polynomial form
static GF Pow2Poly[n], Poly2Pow[n];
#define PrimitivePolynomial  0x1D    // [x^8] + x^4 + x^3 + x^2 + 1

// Simple modulo function for GF operations
static inline void MOD_NN(int& x) {
    x = ((x) >= nn) ? ((x) - nn) : (x);
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
