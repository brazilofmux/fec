#include "rs_tables.h"

RS_TABLES& RS_TABLES::instance() {
    static RS_TABLES instance;
    return instance;
}

RS_TABLES::RS_TABLES()
    : is_initialized_(false)
    , pow2poly_(new GF[GF_N])
    , poly2pow_(new GF[GF_N])
    , split_lo_(new GF[GF_N * 16])
    , split_hi_(new GF[GF_N * 16])
    , mod_table_(nullptr)
{
}

RS_TABLES::~RS_TABLES() {
    delete[] pow2poly_;
    delete[] poly2pow_;
    delete[] split_lo_;
    delete[] split_hi_;
    delete[] mod_table_;
}

void RS_TABLES::ensure_initialized() {
    if (!is_initialized_) {
        initialize_gf();
        initialize_split_tables();
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

void RS_TABLES::initialize_split_tables() {
    // Requires pow2poly_/poly2pow_ already populated by initialize_gf().
    auto gf_mul = [this](int a, int b) -> GF {
        if (a == 0 || b == 0) return 0;
        int sum = poly2pow_[a] + poly2pow_[b];
        if (sum >= nn) sum -= nn;
        return pow2poly_[sum];
    };
    for (int s = 0; s < GF_N; s++) {
        for (int k = 0; k < 16; k++) {
            split_lo_[s * 16 + k] = gf_mul(s, k);
            split_hi_[s * 16 + k] = gf_mul(s, k << 4);
        }
    }
}
