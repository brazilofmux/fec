#include "rs_encoder_general.h"
#include <memory.h>

RS_ENCODER_GENERAL::RS_ENCODER_GENERAL(int tt, int b0)
    : RS_ENCODER_BASE(tt, b0)
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[256 * 2 * tt_])
{
    RSGenPoly();
    RSGenTable();
}

RS_ENCODER_GENERAL::~RS_ENCODER_GENERAL() {
    delete[] ptable;
}

void RS_ENCODER_GENERAL::RSGenPoly() {
    // Initially, g(x) = (X+@^b0)
    gg[0] = pow2poly_[b0_];
    gg[1] = 1;

    for (int i = 2; i <= 2 * tt_; i++) {
        gg[i] = 1;

        // Below multiply (gg[0]+gg[1]*x + ... +gg[i]x^i) by (@^(b0+i-1) + x)
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

void RS_ENCODER_GENERAL::RSGenTable() {
    // Generate table entries for each possible feedback value
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow_[ch];

        // Match the original implementation's indexing and access
        for (int j = 0; j < 2 * tt_; j++) {
            if (gg[j + 1] != GF_INFINITY && feedback != GF_INFINITY) {
                int iMod = mod_nn(gg[j + 1] + feedback);
                ptable[ch * 2 * tt_ + j] = pow2poly_[iMod];
            }
            else {
                ptable[ch * 2 * tt_ + j] = 0;
            }
        }
    }
}

void RS_ENCODER_GENERAL::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    for (int i = 0; i < kk_; i++) {
        GF feedback = bb[0] ^ data[i];
        GF* TableRow = ptable + 2 * tt_ * feedback;

        // Process pairs of bytes
        int j;
        for (j = 0; j < 2 * tt_ - 2; j += 2) {
            bb[j] = bb[j + 1] ^ TableRow[j];
            bb[j + 1] = bb[j + 2] ^ TableRow[j + 1];
        }

        // Handle final bytes
        bb[j] = bb[j + 1] ^ TableRow[j];
        bb[j + 1] = feedback;
    }
}
