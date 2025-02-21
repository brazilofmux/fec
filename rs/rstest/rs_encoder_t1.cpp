#include <memory.h>
#include "projdefs.h"
#include "rs_encoder_t1.h"

RS_ENCODER_T1::RS_ENCODER_T1(int b0)
    : RS_ENCODER_BASE(1, b0)  // tt is always 1
    , ptable(new GF[TABLE_SIZE])
{
    RSGenTable();
}

RS_ENCODER_T1::~RS_ENCODER_T1() {
    delete[] ptable;
}

void RS_ENCODER_T1::RSGenTable() {
    const GF* pow2poly = get_pow2poly();
    const GF* poly2pow = get_poly2pow();

    // Generate the specialized table for tt=1
    // Each feedback value needs 2 entries
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow[ch];
        if (feedback != GF_INFINITY) {
            // First coefficient
            int iMod = mod_nn(feedback + b0_);
            ptable[ch * 2] = pow2poly[iMod];

            // Second coefficient
            iMod = mod_nn(feedback + b0_ + 1);
            ptable[ch * 2 + 1] = pow2poly[iMod];
        }
        else {
            ptable[ch * 2] = 0;
            ptable[ch * 2 + 1] = 0;
        }
    }
}

void RS_ENCODER_T1::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    // For tt=1, we can use a 16-bit LFSR implementation
    // This is highly optimized for the tt=1 case
    USHORT LFSR = 0;
    USHORT* pus = (USHORT*)data;

    // Process pairs of bytes using 16-bit operations
    // For tt=1, we know exactly how many iterations
    int cnt = 126;
    while (cnt--) {
        LFSR ^= *pus++;
        LFSR ^= ptable[2 * (GF)(LFSR)] << 8;
        LFSR ^= ptable[2 * (GF)(LFSR >> 8)];
    }
    *(USHORT*)bb = LFSR;

    // Handle the final bytes
    GF feedback = bb[0] ^ data[252];
    bb[0] = bb[1] ^ ptable[2 * feedback];
    bb[1] = feedback;
}
