#include "rs_flipped_encoder_general.h"
#include <memory>

RS_FLIPPED_ENCODER_GENERAL::RS_FLIPPED_ENCODER_GENERAL(int tt, int b0)
    : RS_ENCODER_BASE(tt, b0)
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[256 * 2 * tt_])
{
    RSGenPoly();
    RSGenTable();
}

RS_FLIPPED_ENCODER_GENERAL::~RS_FLIPPED_ENCODER_GENERAL() {
    delete[] ptable;
}

void RS_FLIPPED_ENCODER_GENERAL::RSGenPoly() {
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

void RS_FLIPPED_ENCODER_GENERAL::RSGenTable() {
    // Generate table entries for each possible feedback value
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow_[ch];

        // Generate coefficients for this feedback value
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

void RS_FLIPPED_ENCODER_GENERAL::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    // Clear parity buffer
    memset(bb, 0, 2 * tt_);

    // Pre-calculate block counts
    const int cnt = 2 * tt_ - 1;
    const int cDO16 = cnt >> 4;
    const int cDO8 = (cnt & 8) >> 3;
    const int cDO4 = (cnt & 4) >> 2;
    const int cDO2 = (cnt & 2) >> 1;
    const int cDO1 = cnt & 1;

    // Process all data bytes
    // The flipped encoder reads feedback from bb[0] (lowest index)
    // and processes the parity buffer in ascending order (0 to 2*tt-1)
    // This is the opposite of the standard encoder, which processes
    // the parity buffer in descending order (2*tt-1 down to 0)
    const GF* data_ptr = data;
    int kk_remaining = kk_;
    while (kk_remaining--) {
        const GF feedback = bb[0] ^ (*data_ptr++);
        GF* TableRow = ptable + 2 * tt_ * feedback;

        int j = 0;

        // Process 16-byte blocks
        int cnt16 = cDO16;
        while (cnt16--) {
            process_block8(bb + j, TableRow + j);
            process_block8(bb + j + 8, TableRow + j + 8);
            j += 16;
        }

        // Process 8-byte block if needed
        if (cDO8) {
            process_block8(bb + j, TableRow + j);
            j += 8;
        }

        // Process 4-byte block if needed
        if (cDO4) {
            process_block4(bb + j, TableRow + j);
            j += 4;
        }

        // Process 2 bytes if needed
        if (cDO2) {
            process_byte1(bb, TableRow, j);
            process_byte1(bb, TableRow, j + 1);
            j += 2;
        }

        // Process 1 byte if needed
        if (cDO1) {
            process_byte1(bb, TableRow, j);
            j++;
        }

        // Final feedback byte
        bb[j] = feedback;
    }
}
