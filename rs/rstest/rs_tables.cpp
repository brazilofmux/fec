#include "rs_tables.h"

RS_TABLES& RS_TABLES::instance() {
    static RS_TABLES instance;
    return instance;
}

RS_TABLES::RS_TABLES() : is_initialized_(false) {
    pow2poly_ = new GF[n];
    poly2pow_ = new GF[n];
}

RS_TABLES::~RS_TABLES() {
    delete[] pow2poly_;
    delete[] poly2pow_;
}

void RS_TABLES::ensure_initialized() {
    if (!is_initialized_) {
        initialize_gf();
        is_initialized_ = true;
    }
}

void RS_TABLES::initialize_gf() {
    int iPoly = 1;
    const int PrimitivePolynomial = 0x1D;

    // Handle special infinite case
    pow2poly_[n - 1] = 0;
    poly2pow_[0] = GF_INFINITY;

    // Handle Pow2Poly[0..n-2] and Poly2Pow[1..n-1]
    for (int i = 0; i < n - 1; i++) {
        pow2poly_[i] = iPoly;
        poly2pow_[iPoly] = i;
        iPoly <<= 1;
        if (iPoly > n - 1) {
            iPoly ^= PrimitivePolynomial;
            iPoly &= 0xFF;
        }
    }
}
