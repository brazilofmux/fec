#include "rs_encoder_general_descending.h"
#include <memory>

RS_ENCODER_GENERAL_DESCENDING::RS_ENCODER_GENERAL_DESCENDING(int tt, int b0)
    : RS_ENCODER_BASE(tt, b0)
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[256 * (2 * tt_ + 1)])
{
    RSGenPoly();
    RSGenTable();
}

RS_ENCODER_GENERAL_DESCENDING::~RS_ENCODER_GENERAL_DESCENDING() {
    delete[] ptable;
}

void RS_ENCODER_GENERAL_DESCENDING::RSGenPoly() {
    // Initially, g(x) = (X+@^b0)
    gg[0] = pow2poly_[b0_];
    gg[1] = 1;

    for (int i = 2; i <= 2 * tt_; i++) {
        gg[i] = 1;

        // Multiply (gg[0]+gg[1]*x + ... +gg[i]x^i) by (@^(b0+i-1) + x)
        for (int j = i - 1; j > 0; j--) {
            if (gg[j] != 0) {
                int iMod = mod_nn(poly2pow_[gg[j]] + b0_ + i - 1);
                gg[j] = gg[j - 1] ^ pow2poly_[iMod];
            }
            else {
                gg[j] = gg[j - 1];
            }
        }

        // gg[0] can never be zero
        int iMod = mod_nn(poly2pow_[gg[0]] + b0_ + i - 1);
        gg[0] = pow2poly_[iMod];
    }

    // Convert gg[] to power form for quicker encoding
    for (int i = 0; i <= 2 * tt_; i++) {
        gg[i] = poly2pow_[gg[i]];
    }
}

void RS_ENCODER_GENERAL_DESCENDING::RSGenTable() {
    const int NON_SYM_LEN = 2 * tt_ + 1;
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow_[ch];
        for (int j = 0; j <= 2*tt_; j++) {
            if (gg[j] != GF_INFINITY && feedback != GF_INFINITY) {
                int iMod = mod_nn(gg[j] + feedback);
                ptable[ch * NON_SYM_LEN + j] = pow2poly_[iMod];
            }
            else {
                ptable[ch * NON_SYM_LEN + j] = 0;
            }
        }
    }
}

void RS_ENCODER_GENERAL_DESCENDING::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    const int NON_SYM_LEN = 2 * tt_ + 1;
    memset(bb, 0, 2 * tt_ * sizeof(GF));
    for (int i = kk_ - 1; i >= 0; i--) {
        const GF feedback = bb[2 * tt_ - 1] ^ data[i];
        GF* TableRow = ptable + NON_SYM_LEN * feedback;
        for (int j = 2 * tt_ - 1; j > 0; j--) {
            process_byte1(bb, TableRow, j);
        }
        bb[0] = TableRow[0];
    }
}
