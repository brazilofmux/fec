#include "rs_tables.h"

RS_TABLES& RS_TABLES::instance() {
    static RS_TABLES instance;
    return instance;
}

RS_TABLES::RS_TABLES()
    : is_initialized_(false)
    , pow2poly_(new GF[GF_N])
    , poly2pow_(new GF[GF_N])
    , mod_table_(nullptr)
{
}

RS_TABLES::~RS_TABLES() {
    delete[] pow2poly_;
    delete[] poly2pow_;
    delete[] mod_table_;
}

void RS_TABLES::ensure_initialized() {
    if (!is_initialized_) {
        initialize_gf();
        is_initialized_ = true;
    }
}

void RS_TABLES::initialize_gf() {
    initialize_mod_table();

    int iPoly = 1;
    const int PrimitivePolynomial = 0x1D;  // [x^8] + x^4 + x^3 + x^2 + 1

    // Handle special infinite case
    pow2poly_[GF_N - 1] = 0;
    poly2pow_[0] = GF_INFINITY;

    // Handle Pow2Poly[0..n-2] and Poly2Pow[1..n-1]
    for (int i = 0; i < GF_N - 1; i++) {
        pow2poly_[i] = iPoly;
        poly2pow_[iPoly] = i;
        iPoly <<= 1;
        if (iPoly > GF_N - 1) {
            iPoly ^= PrimitivePolynomial;
            iPoly &= 0xFF;
        }
    }
}
