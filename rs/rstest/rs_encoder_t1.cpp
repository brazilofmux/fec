#include <memory>
#include "projdefs.h"
#include "rs_encoder_t1.h"

RS_ENCODER_T1::RS_ENCODER_T1(int b0)
    : RS_ENCODER_BASE(1, b0)  // tt is always 1
    , ptable(new GF[TABLE_SIZE])
    , gg(2 * tt_ + 1)  // Size based on tt_
{
    RSGenPoly();
    RSGenTable();
}

RS_ENCODER_T1::~RS_ENCODER_T1() {
    delete[] ptable;
}

void RS_ENCODER_T1::RSGenPoly() {
    // Initially, g(x) = (X+@^b0)
    gg[0] = Pow2Poly[b0_];
    gg[1] = 1;

    for (int i = 2; i <= 2 * tt_; i++) {
        gg[i] = 1;

        // Multiply (gg[0]+gg[1]*x + ... +gg[i]x^i) by (@^(b0+i-1) + x)
        for (int j = i - 1; j > 0; j--) {
            if (gg[j] != 0) {
                int iMod = mod_nn(Poly2Pow[gg[j]] + b0_ + i - 1);
                gg[j] = gg[j - 1] ^ Pow2Poly[iMod];
            }
            else {
                gg[j] = gg[j - 1];
            }
        }

        // gg[0] can never be zero
        int iMod = mod_nn(Poly2Pow[gg[0]] + b0_ + i - 1);
        gg[0] = Pow2Poly[iMod];
    }

    // Convert gg[] to power form for quicker encoding
    for (int i = 0; i <= 2 * tt_; i++) {
        gg[i] = Poly2Pow[gg[i]];
    }
}

void RS_ENCODER_T1::RSGenTable() {
    // Generate table entries for each possible feedback value
    for (int ch = 0; ch < 256; ch++) {
        int feedback = Poly2Pow[ch];

        // For tt=1, we generate two coefficients per feedback value
        for (int j = 0; j < 2 * tt_; j++) {
            if (gg[j + 1] != GF_INFINITY && feedback != GF_INFINITY) {
                int iMod = mod_nn(gg[j + 1] + feedback);
                ptable[ch * 2 * tt_ + j] = Pow2Poly[iMod];
            }
            else {
                ptable[ch * 2 * tt_ + j] = 0;
            }
        }
    }
}

void RS_ENCODER_T1::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    // Optimized encoder for tt=1 using 16-bit LFSR
    USHORT LFSR = 0;
    USHORT* pus = (USHORT*)data;
    const int pairs = (kk_ - 1) / 2;  // Process pairs of bytes

    for (int cnt = 0; cnt < pairs; cnt++) {
        LFSR ^= *pus++;
        LFSR ^= ptable[2 * (GF)(LFSR)] << 8;
        LFSR ^= ptable[2 * (GF)(LFSR >> 8)];
    }
    *(USHORT*)bb = LFSR;

    // Handle final byte if kk_ is odd
    if (kk_ % 2) {
        GF feedback = bb[0] ^ data[kk_ - 1];
        bb[0] = bb[1] ^ ptable[2 * feedback];
        bb[1] = feedback;
    }
}
