#include "rs_standard_encoder_general.h"
#include <memory>

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

#if 0
void RS_STANDARD_ENCODER_GENERAL::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    const int NON_SYM_LEN = 2 * tt_ + 1;
    memset(bb, 0, 2 * tt_ * sizeof(GF));

    for (int i = kk_ - 1; i >= 0; i--) {
        const GF feedback = bb[2 * tt_ - 1] ^ data[i];
        GF* TableRow = ptable + NON_SYM_LEN * feedback;

        // Process bytes in blocks where possible
        int j = 2 * tt_ - 1;

        // Process 8-byte blocks
        while (j >= 8) {
            process_block8(&bb[j-7], &TableRow[j-7]);
            j -= 8;
        }

        // Process 4-byte blocks
        if (j >= 4) {
            process_block4(&bb[j-3], &TableRow[j-3]);
            j -= 4;
        }

        // Process remaining bytes individually
        while (j > 0) {
            process_byte1(bb, TableRow, j);
            j--;
        }

        bb[0] = TableRow[0];
    }
}
#elif 0
void RS_STANDARD_ENCODER_GENERAL::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    const int NON_SYM_LEN = 2 * tt_ + 1;
    memset(bb, 0, 2 * tt_ * sizeof(GF));

    for (int i = kk_ - 1; i >= 0; i--) {
        const GF feedback = bb[2 * tt_ - 1] ^ data[i];
        GF* TableRow = ptable + NON_SYM_LEN * feedback;

        // Specialized code paths based on the common Reed-Solomon code sizes
        switch (2 * tt_) {
            case 32: // t=16, common for many Reed-Solomon implementations
                process_block8(&bb[24], &TableRow[24]);
                process_block8(&bb[16], &TableRow[16]);
                process_block8(&bb[8], &TableRow[8]);
                process_block8(&bb[0], &TableRow[0]);
                break;

            case 16: // t=8, another common size
                process_block8(&bb[8], &TableRow[8]);
                process_block8(&bb[0], &TableRow[0]);
                break;

            case 10: // t=5, used in some applications
                process_block8(&bb[2], &TableRow[2]);
                process_byte1(bb, TableRow, 1);
                bb[0] = TableRow[0];
                break;

            default: // General case with dynamic unrolling
                int j = 2 * tt_ - 1;

                // Process 8-byte blocks
                while (j >= 8) {
                    process_block8(&bb[j-7], &TableRow[j-7]);
                    j -= 8;
                }

                // Process 4-byte blocks
                if (j >= 4) {
                    process_block4(&bb[j-3], &TableRow[j-3]);
                    j -= 4;
                }

                // Process remaining bytes individually
                while (j > 0) {
                    process_byte1(bb, TableRow, j);
                    j--;
                }

                bb[0] = TableRow[0];
                break;
        }
    }
}
#elif 1
void RS_STANDARD_ENCODER_GENERAL::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    const int NON_SYM_LEN = 2 * tt_ + 1;
    const int parity_bytes = 2 * tt_;

    // Initialize parity bytes to zero
    memset(bb, 0, parity_bytes * sizeof(GF));

    for (int i = kk_ - 1; i >= 0; i--) {
        // Calculate feedback
        const GF feedback = bb[parity_bytes - 1] ^ data[i];

        // Fast path for zero feedback (common case)
        if (feedback == 0) {
            // Just shift the registers
            for (int j = parity_bytes - 1; j > 0; j--) {
                bb[j] = bb[j - 1];
            }
            bb[0] = 0;
            continue;
        }

        GF* TableRow = ptable + NON_SYM_LEN * feedback;

        // Use optimized processing based on parity size and alignment
        // We need to handle alignment properly for the block operations

        if ((parity_bytes % 8) == 0 && parity_bytes >= 8) {
            // Perfectly aligned for 8-byte blocks
            for (int j = parity_bytes - 8; j >= 0; j -= 8) {
                process_block8(&bb[j], &TableRow[j]);
            }
        }
        else if ((parity_bytes % 4) == 0 && parity_bytes >= 4) {
            // Perfectly aligned for 4-byte blocks
            for (int j = parity_bytes - 4; j >= 0; j -= 4) {
                process_block4(&bb[j], &TableRow[j]);
            }
        }
        else {
            // Mixed strategy for non-aligned sizes
            int j = parity_bytes - 1;

            // Process 8-byte blocks where possible
            while (j >= 7) {
                process_block8(&bb[j-7], &TableRow[j-7]);
                j -= 8;
            }

            // Process 4-byte blocks where possible
            if (j >= 3) {
                process_block4(&bb[j-3], &TableRow[j-3]);
                j -= 4;
            }

            // Process remaining bytes individually
            while (j >= 0) {
                if (j > 0) {
                    process_byte1(bb, TableRow, j);
                } else {
                    bb[0] = TableRow[0];
                }
                j--;
            }
        }
    }
}
#endif
