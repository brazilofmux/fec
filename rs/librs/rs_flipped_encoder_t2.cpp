#include "rs_flipped_encoder_t2.h"
#include <cstring>
#include <memory>

RS_FLIPPED_ENCODER_T2::RS_FLIPPED_ENCODER_T2(int b0)
    : RS_ENCODER_BASE(2, b0)  // tt is always 2
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[TABLE_SIZE])
{
    RSGenPoly();
    RSGenTable();
}

RS_FLIPPED_ENCODER_T2::~RS_FLIPPED_ENCODER_T2() {
    delete[] ptable;
}

void RS_FLIPPED_ENCODER_T2::RSGenPoly() {
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

void RS_FLIPPED_ENCODER_T2::RSGenTable() {
    // Generate table entries for each possible feedback value
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow_[ch];

        // For tt=2, we generate four coefficients per feedback value
        for (int j = 0; j < 4; j++) {
            if (gg[j + 1] != GF_INFINITY && feedback != GF_INFINITY) {
                int iMod = mod_nn(gg[j + 1] + feedback);
                ptable[ch * 4 + j] = pow2poly_[iMod];
            }
            else {
                ptable[ch * 4 + j] = 0;
            }
        }
    }
}

void RS_FLIPPED_ENCODER_T2::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    // Initialize LFSR to 0
    uint32_t LFSR = 0;
    const GF* pdata = data;

    // Loads a 4-byte ptable row as one word. memcpy compiles to a single
    // load instruction at -O1+ while staying alignment- and aliasing-clean.
    auto table_row32 = [this](GF idx) {
        uint32_t row;
        memcpy(&row, ptable + 4 * idx, sizeof(row));
        return row;
    };

    // Process data in 4-byte chunks for most of the data
    int remaining = kk_;
    while (remaining >= 4) {
        uint32_t chunk;
        memcpy(&chunk, pdata, sizeof(chunk));  // XOR in next 4 bytes of data
        pdata += 4;
        LFSR ^= chunk;

        // Process each byte through the LFSR
        LFSR ^= table_row32((GF)LFSR) << 8;
        LFSR = rotr32(LFSR, 8);
        LFSR ^= table_row32((GF)LFSR) << 8;
        LFSR = rotr32(LFSR, 8);
        LFSR ^= table_row32((GF)LFSR) << 8;
        LFSR = rotr32(LFSR, 8);
        LFSR ^= table_row32((GF)LFSR) << 8;
        LFSR = rotr32(LFSR, 8);

        remaining -= 4;
    }

    // Store current LFSR state into output
    memcpy(bb, &LFSR, sizeof(LFSR));

    // Handle remaining bytes one at a time
    while (remaining > 0) {
        const GF feedback = bb[0] ^ *pdata++;
        bb[0] = bb[1] ^ ptable[4 * feedback];
        bb[1] = bb[2] ^ ptable[4 * feedback + 1];
        bb[2] = bb[3] ^ ptable[4 * feedback + 2];
        bb[3] = feedback;
        remaining--;
    }
}
