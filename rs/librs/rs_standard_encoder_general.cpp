#include "rs_standard_encoder_general.h"
#include <cstring>

RS_STANDARD_ENCODER_GENERAL::RS_STANDARD_ENCODER_GENERAL(int tt, int b0)
    : RS_ENCODER_BASE(tt, b0)
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[256 * (2 * tt_ + 1)])
{
    RSGenPoly();
    RSGenTable();
}

RS_STANDARD_ENCODER_GENERAL::~RS_STANDARD_ENCODER_GENERAL() {
    delete[] ptable;
}

void RS_STANDARD_ENCODER_GENERAL::RSGenPoly() {
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

void RS_STANDARD_ENCODER_GENERAL::RSGenTable() {
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

void RS_STANDARD_ENCODER_GENERAL::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    const int NON_SYM_LEN = 2 * tt_ + 1;
    const int parity_bytes = 2 * tt_;

    memset(bb, 0, parity_bytes * sizeof(GF));

    for (int i = kk_ - 1; i >= 0; i--) {
        const GF feedback = bb[parity_bytes - 1] ^ data[i];

        // Fast path for zero feedback: just shift the registers.
        if (feedback == 0) {
            for (int j = parity_bytes - 1; j > 0; j--) {
                bb[j] = bb[j - 1];
            }
            bb[0] = 0;
            continue;
        }

        GF* TableRow = ptable + NON_SYM_LEN * feedback;

        // The blocks read one byte below where they write (bb[j] depends on
        // bb[j-1]), so a block may only be placed where its read window stays
        // inside the buffer: j-7 >= 1 for 8-byte blocks, j-3 >= 1 for 4-byte
        // blocks. bb[0] has no bb[-1] term and is written directly.
        int j = parity_bytes - 1;

        while (j >= 8) {
            process_block8(&bb[j - 7], &TableRow[j - 7]);
            j -= 8;
        }

        if (j >= 4) {
            process_block4(&bb[j - 3], &TableRow[j - 3]);
            j -= 4;
        }

        while (j > 0) {
            process_byte1(bb, TableRow, j);
            j--;
        }

        bb[0] = TableRow[0];
    }
}
