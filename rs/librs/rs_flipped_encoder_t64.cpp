#include "rs_flipped_encoder_t64.h"
#include <memory>

RS_FLIPPED_ENCODER_T64::RS_FLIPPED_ENCODER_T64(int b0)
    : RS_ENCODER_BASE(64, b0)  // tt is always 64
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[TABLE_SIZE])
{
    RSGenPoly();
    RSGenTable();
}

RS_FLIPPED_ENCODER_T64::~RS_FLIPPED_ENCODER_T64() {
    delete[] ptable;
}

void RS_FLIPPED_ENCODER_T64::RSGenPoly() {
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

void RS_FLIPPED_ENCODER_T64::RSGenTable() {
    // Generate table entries for each possible feedback value
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow_[ch];

        // For tt=64, we generate 128 coefficients per feedback value
        for (int j = 0; j < 128; j++) {
            if (gg[j + 1] != GF_INFINITY && feedback != GF_INFINITY) {
                int iMod = mod_nn(gg[j + 1] + feedback);
                ptable[ch * 128 + j] = pow2poly_[iMod];
            }
            else {
                ptable[ch * 128 + j] = 0;
            }
        }
    }
}

void RS_FLIPPED_ENCODER_T64::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    // Clear parity buffer
    memset(bb, 0, 2 * tt_);

    // Process all data bytes
    for (int i = 0; i < kk_; i++) {
        const GF feedback = bb[0] ^ data[i];
        GF* TableRow = ptable + 128 * feedback;

        // Process fifteen 8-byte blocks using 64-bit operations
        for (int j = 0; j <= 112; j += 8) {
            process_block8(bb + j, TableRow + j);
        }

        // Process one 4-byte block
        process_block4(bb + 120, TableRow + 120);

        // Process remaining bytes individually
        process_byte1(bb, TableRow, 124);
        process_byte1(bb, TableRow, 125);
        process_byte1(bb, TableRow, 126);

        // Final byte is feedback
        bb[127] = feedback;
    }
}
