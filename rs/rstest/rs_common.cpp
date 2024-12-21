#include "rs_common.h"

// Define the global tables used by both implementations
GF Pow2Poly[n];
GF Poly2Pow[n];

// Implementation of RSGenField from the original rs.cpp
static void RSGenField(void) {
    int i;
    int iPoly = 1;
    const int PrimitivePolynomial = 0x1D;    // [x^8] + x^4 + x^3 + x^2 + 1

    // Handle special infinite case
    Pow2Poly[n-1] = 0;
    Poly2Pow[0] = GF_INFINITY;

    // Handle Pow2Poly[0..n-2] and Poly2Pow[1..n-1]
    for (i = 0; i < n-1; i++) {
        Pow2Poly[i] = iPoly;
        Poly2Pow[iPoly] = i;
        iPoly <<= 1;
        if (iPoly > n-1) {
            iPoly ^= PrimitivePolynomial;
            iPoly &= 0xFF;
        }
    }
}

static bool initialized = false;

void RS_Init() {
    if (!initialized) {
        RSGenField();
        initialized = true;
    }
}