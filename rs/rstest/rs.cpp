// Implementation of Reed-Solomon Error Correcting Codes including
// Errors only and Errors+Erasures decoding.
// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All right reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//
// Profiling on a AWS t2.large instance yield the following performance
// measured in microseconds per codeword:
//
//             ---Errors Only---  -Errors+Erasures-
//             ----Decoding-----  ----Decoding-----
// tt  Encode  w/errors  without  w/errors  without
//  1    0.60     1.20    -----    -----
//  2    0.50     2.80    -----    -----
//  3    2.10     5.10    -----    -----
//  4    2.40     4.10    -----    -----
//  5    2.80     6.20    -----    -----
//  6    2.40    -----    -----    -----
//  7    2.10    -----    -----    -----
//  8    2.20    10.70    -----    -----
// 16    1.60    -----    -----    -----
// 24    2.00    -----    -----    -----
// 32    1.50    -----    -----    -----    -----
// 40    2.00    -----    -----    -----
// 48    2.10    -----    -----    -----
// 56    2.10    80.90    -----    -----
// 64    1.30    93.00    -----    -----
//
// A codeword is 255 bytes. The number of parity bytes is 2*tt, and the
// number of data bytes is 255-2*tt.
//
// On Atom Chromebook:
//
//             ---Errors Only---  -Errors+Erasures-
//             ----Decoding-----  ----Decoding-----
// tt  Encode  w/errors  without  w/errors  without
//  1    0.75      0.89      0.92    -----    -----
//  2    0.71      1.66      1.62    -----    -----
//  3    1.62      2.26      2.26    -----    -----
//  4    1.56      3.02      3.03    -----    -----
//  5    1.68      3.93      3.93    -----    -----
//  6    1.77      4.77      4.79    -----    -----
//  7    1.66      5.91      5.88    -----    -----
//  8    1.60      6.95      6.96    -----    -----
// 16    1.41     14.95     14.96    -----    -----
// 24    1.90     22.91     22.86    -----    -----
// 32    1.62     30.98     31.13    -----    -----
// 40    1.81     38.86     39.66    -----    -----
// 48    1.73     46.54     46.60    -----    -----
// 56    1.77     54.91     54.77    -----    -----
// 64    1.76     62.49     63.00    -----    -----
//
// On AWS m5.xlarge:
//
//             ---Errors Only---  -Errors+Erasures-
//             ----Decoding-----  ----Decoding-----
// tt  Encode  w/errors  without  w/errors  without
//  1    0.80      0.97      0.98    -----    -----
//  2    0.70      2.77      2.81    -----    -----
//  3    1.48      3.54      3.54    -----    -----
//  4    1.93      4.05      4.05    -----    -----
//  5    1.59      4.52      4.52    -----    -----
//  6    2.04      5.22      5.21    -----    -----
//  7    1.53      5.95      6.15    -----    -----
//  8    1.93      6.70      6.67    -----    -----
// 16    1.75     15.17     15.30    -----    -----
// 24    1.69     24.42     24.43    -----    -----
// 32    1.39     33.65     33.68    -----    -----
// 40    1.47     43.04     42.97    -----    -----
// 48    1.45     52.12     52.14    -----    -----
// 56    1.32     61.59     61.58    -----    -----
// 64    1.29     70.59     70.63    -----    -----

#include <stdlib.h>
#include <memory.h>
#include <cstdio>
#include "projdefs.h"
#include "rs.h"
#include "rs_debug.h"

#define NO_PRINT

#define n   256         // n  = 2^8     size of the field
#define nn  255         // nn = 2^8-1   length of codeword
#define MAX_TT (127)    // maximum number of correctable errors (nn/2).
#define MIN_TT (1)      // minimum number of correctable errors (1).
#define MAX_KK (253)    // Maximum number of data symbols (nn-2*TT_MIN).
#define MIN_KK (1)      // Minimum number of data symbols (nn-2*TT_MAX).

GF RS_ENCODER::ModTable[4*nn+1];
GF *RS_ENCODER::pModTable;
#define MOD_NN(x) { x = pModTable[x]; }

RS_ENCODER::RS_ENCODER(int CorrectableErrors)
{
    tt = CorrectableErrors;
    kk = nn-2*tt;

    // To choose a symetrically generator polynomial, we
    // pick b0 = [(2^m-1)-(Dmin-2)]/2, for even Dmin-1.
    // and  b0 = -(Dmin-2)/2, for odd Dmin-1.
    //
    // This ensures that g(x) = x^2t*g(1/x), which makes the coefficients
    // of g(x) symmetric.
    //
    // For RS codes, Dmin-1 = 2t which is always even, so
    // b0 = [(2^m-1)-(Dmin-2)]/2 = [nn - (2t-1)]/2 = [nn-2t+1]/2
    //    = (kk+1)/2
    //
#ifdef SYMMETRICAL_GENERATOR
    b0 = (kk+1)/2;
#else
    b0 = 254;
#endif

    RSGenPoly();
    RSGenTable();
    RsVerify::verify_generator(gg, tt);
}

RS_ENCODER::~RS_ENCODER(void)
{
    delete ptable;
}

//
//   Obtain the generator polynomial of the tt-error correcting, length
// nn=(2^8-1) Reed Solomon code from the product of (X+@^(b0+i)),
// i = 0, ... ,(2*tt-1)
//
// Examples:
//      If b0 = 1, tt = 1. deg(g(x)) = 2*tt = 2, then
//          g(x) = (x+@) (x+@^2)
//
//      If b0 = 0, tt = 2. deg(g(x)) = 2*tt = 4, then
//          g(x) = (x+1) (x+@) (x+@^2) (x+@^3)
//
void RS_ENCODER::RSGenPoly(void)
{
    int i, j, iMod;

    // Initially, g(x) = (X+@^b0)
    //
    gg[0] = Pow2Poly[b0];
    gg[1] = 1;
    for (i = 2; i <= nn-kk; i++)
    {
        gg[i] = 1;

        // Below multiply (gg[0]+gg[1]*x + ... +gg[i]x^i) by (@^(b0+i-1) + x)
        //
        for (j = i-1; j > 0; j--)
        {
            if (gg[j] != 0)
            {
                iMod = Poly2Pow[gg[j]]+b0+i-1;
                MOD_NN(iMod);
                gg[j] = gg[j-1]^Pow2Poly[iMod];
            }
            else
            {
                gg[j] = gg[j-1];
            }
        }

        // gg[0] can never be zero
        //
        iMod = Poly2Pow[gg[0]]+b0+i-1;
        MOD_NN(iMod);
        gg[0] = Pow2Poly[iMod];
    }

    // Convert gg[] to power form for quicker encoding
    //
    for (i = 0; i <= 2*tt; i++)
    {
        gg[i] = Poly2Pow[gg[i]];
    }
}

#define SYM_GEN_ENC_LEN (tt-1)
#define SYM_GEN_LEN     (2*tt)
#define NON_SYM_LEN     (2*tt+1)

void RS_ENCODER::RSGenTable(void)
{
#ifdef SYMMETRICAL_GENERATOR
#ifdef SYMMETRICAL_ENCODE
    // Could be tt-1, but we may want a little bigger for alignment.
    //
    ptable = new GF[256*SYM_GEN_ENC_LEN];
#else
    // Could be 2*tt-1, but we may want a little bigger for alignment.
    //
    ptable = new GF[256*SYM_GEN_LEN];
#endif
#else
    // Could be 2*tt+1, but we may want a little bigger for alignment.
    //
    ptable = new GF[256*NON_SYM_LEN];
#endif

    // For every possible feedback
    //
    int ch, iMod;
    for (ch = 0; ch < 256; ch++)
    {
        int feedback = Poly2Pow[ch];
        // For about half the terms in symmetrical gg[]. That is, for terms
        // x^1, ..., x^tt. The coefficient for x^0 is always 1 for
        // symmetrical gg, and coefficients x^(tt+1), ..., x^(2*tt-1) are
        // equal to x^tt-1, ..., x^1, respectively. x^(2*tt) is always 1 even
        // non-symmetrical gg[], and doesn't correspond to any of the taps of
        // the feedback register, and so doesn't affect anything. x^tt is
        // unique.
        //
        int j;
#ifdef SYMMETRICAL_GENERATOR
#ifdef SYMMETRICAL_ENCODE
        for (j = 0; j < tt; j++)
#else
        for (j = 0; j < 2*tt; j++)
#endif
#else
        for (j = 0; j <= 2*tt; j++)
#endif
#ifdef SYMMETRICAL_GENERATOR
        {
            if (gg[j+1] != GF_INFINITY && feedback != GF_INFINITY)
            {
                iMod = gg[j+1]+feedback;
                MOD_NN(iMod);
                // table[ch][j] = Pow2Poly[iMod];
#ifdef SYMMETRICAL_ENCODE
                ptable[ch*SYM_GEN_ENC_LEN+j] = Pow2Poly[iMod];
#else
                ptable[ch*SYM_GEN_LEN+j] = Pow2Poly[iMod];
#endif
            }
            else
            {
                //table[ch][j] = 0;
#ifdef SYMMETRICAL_ENCODE
                ptable[ch*SYM_GEN_ENC_LEN+j] = 0;
#else
                ptable[ch*SYM_GEN_LEN+j] = 0;
#endif
            }
        }
#else
        {
            if (gg[j] != GF_INFINITY && feedback != GF_INFINITY)
            {
                iMod = gg[j]+feedback;
                MOD_NN(iMod);
                //table[ch][j] = Pow2Poly[iMod];
                ptable[ch*NON_SYM_LEN+j] = Pow2Poly[iMod];
            }
            else
            {
                //table[ch][j] = 0;
                ptable[ch*NON_SYM_LEN+j] = 0;
            }
        }
#endif
    }
}

//
// Take the string of symbols in data[i], i=0..(k-1) and encode systematically
// to produce 2*tt parity symbols in bb[0]..bb[2*tt-1]. data[] is input and
// bb[] is output in polynomial form. Encoding is done by using a table-driven
// feedback shift register with appropriate connections specified by the
// elements of gg[], which was generated above.
// Codeword is c(X) = data(X)*X^(nn-kk) + bb(X)
//

#if   defined(SYMMETRICAL_GENERATOR)  && defined(ASCENDING_ENCODE)

#define DO1(x) { bb[x] = bb[x+1]^TableRow[x]; }
#define DO4(x) { \
    *((UINT32 *)(bb+x)) = *((UINT32 *)(bb+x+1)) \
                        ^ (*((UINT32 *)(TableRow+x))); }
#define DO8(x) { \
    *((UINT64 *)(bb+x)) = *((UINT64 *)(bb+x+1)) \
                        ^ (*((UINT64 *)(TableRow+x))); }
#define DOFB(x) bb[x] = feedback;

#elif defined(SYMMETRICAL_GENERATOR)  && !defined(ASCENDING_ENCODE)

#define DO1(x) {  bb[x] = bb[x-1]^TableRow[x-1]; }
#define DO4(x) { \
    *((UINT32 *)(bb+x)) = *((UINT32 *)(bb+x-1)) \
                        ^ (*((UINT32 *)(TableRow+x-1))); }
#define DO8(x) { \
    *((UINT64 *)(bb+x)) = *((UINT64 *)(bb+x-1)) \
                        ^ (*((UINT64 *)(TableRow+x-1))); }
#define DOFB(x) bb[x] = feedback;

#elif !defined(SYMMETRICAL_GENERATOR) &&  defined(ASCENDING_ENCODE)
#error "Not possible."
#elif !defined(SYMMETRICAL_GENERATOR) && !defined(ASCENDING_ENCODE)

#define DO1(x) {  bb[x] = bb[x-1]^TableRow[x]; }
#define DO4(x) { \
    *((UINT32 *)(bb+x)) = *((UINT32 *)(bb+x-1)) \
                        ^ (*((UINT32 *)(TableRow+x))); }
#define DO8(x) { \
    *((UINT64 *)(bb+x)) = *((UINT64 *)(bb+x-1)) \
                        ^ (*((UINT64 *)(TableRow+x))); }
#define DOFB(x) bb[x] = TableRow[x];

#endif

void RS_ENCODER::RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT])
{
#if 0
    switch (tt)
    {
    case 1:  RSEncode1 (data, bb); break;
    case 2:  RSEncode2 (data, bb); break;
    case 4:  RSEncode4 (data, bb); break;
    case 8:  RSEncode8 (data, bb); break;
    case 14: RSEncode14(data, bb); break;
    case 16: RSEncode16(data, bb); break;
    case 32: RSEncode32(data, bb); break;
    case 64: RSEncode64(data, bb); break;
    default: RSEncodeGeneral(data, bb); break;
    }
#else
    RSEncodeGeneral(data, bb);
#endif
}

void RS_ENCODER::RSEncodeGeneral(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 2*tt);

#ifdef SYMMETRICAL_ENCODE // Symmetrical encoders
#ifdef ASCENDING_ENCODE // Symmetrical, ascending encoder

    for (int i = 0; i < kk; i++)
    {
        GF feedback = bb[0]^data[i];
        GF *TableRow = table[feedback];
        for (int j = 0; j < tt; j++)
        {
            bb[j] = bb[j+1]^TableRow[j];
        }
        for (j = tt; j < 2*tt-1; j++)
        {
            bb[j] = bb[j+1]^TableRow[2*tt-j-2];
        }
        bb[2*tt-1] = feedback;
    }

#else // ASCENDING_ENCODE: Symmetrical, descending encoder

    for (int i = kk-1; i >= 0; i--)
    {
        GF feedback = bb[2*tt-1]^data[i];
        GF *TableRow = table[feedback];
        for (int j = 2*tt-1; j > tt; j--)
        {
            bb[j] = bb[j-1]^TableRow[2*tt-j-1];
        }
        for (j = tt; j > 0; j--)
        {
            bb[j] = bb[j-1]^TableRow[j-1];
        }
        bb[0] = feedback;
    }

#endif // ASCENDING_ENCODE
#else // Non-Symmetrical encoders (faster, but larger lookup table
#ifdef ASCENDING_ENCODE // Non-Symmetrical, ascending encoders
#if defined(ASSEMBLY_ENCODE)

    if (14 <= this->tt)
    {
        if (32 <= this->tt)
        {
            if (64 == this->tt)
            {
                for (int i = 0; i < this->kk; i++)
                {
                    const GF feedback = bb[0] ^ data[i];
                    GF* TableRow = (this->ptable) + 2 * this->tt * feedback;

                    DO8(0); DO8(8); DO8(16); DO8(24);
                    DO8(32); DO8(40); DO8(48); DO8(56);
                    DO8(64); DO8(72); DO8(80); DO8(88);
                    DO8(96); DO8(104); DO8(112);

                    DO4(120);

                    DO1(124); DO1(125); DO1(126);

                    DOFB(127);
                }
                return;
            }
            else if (32 == this->tt)
            {
                for (int i = 0; i < this->kk; i++)
                {
                    const GF feedback = bb[0] ^ data[i];
                    GF* TableRow = (this->ptable) + 2 * this->tt * feedback;

                    DO8(0); DO8(8); DO8(16); DO8(24);
                    DO8(32); DO8(40); DO8(48);

                    DO4(56);

                    DO1(60); DO1(61); DO1(62);

                    DOFB(63);
                }
                return;
            }
        }
        else if (16 == this->tt)
        {
            for (int i = 0; i < this->kk; i++)
            {
                const GF feedback = bb[0] ^ data[i];
                GF* TableRow = (this->ptable) + 2 * this->tt * feedback;

                DO8(0); DO8(8); DO8(16);

                DO4(24);

                DO1(28); DO1(29); DO1(30);

                DOFB(31);
            }
            return;
        }
        else if (14 == this->tt)
        {
            for (int i = 0; i < this->kk; i++)
            {
                const GF feedback = bb[0] ^ data[i];
                GF* TableRow = (this->ptable) + 2 * this->tt * feedback;

                DO8(0); DO8(8); DO8(16);

                DO1(24); DO1(25); DO1(26);

                DOFB(27);
            }
            return;
        }
    }
    if (2 <= this->tt)
    {
        if (4 == this->tt)
        {
            for (int i = 0; i < this->kk; i++)
            {
                const GF feedback = bb[0] ^ data[i];
                GF* TableRow = (this->ptable) + 2 * this->tt * feedback;

                DO4(0);

                DO1(4); DO1(5); DO1(6);

                DOFB(7);
            }
            return;
        }
        if (2 == this->tt)
        {
            UINT32 LFSR = 0;
            int cnt = this->kk;
            int i = 0;
            while (cnt >= 4)
            {
                LFSR ^= *(UINT32*)(data + i);
                LFSR ^= (((UINT32*)(this->ptable))[(GF)(LFSR)]) << 8;
                LFSR = ROTR32(LFSR, 8);
                LFSR ^= (((UINT32*)(this->ptable))[(GF)(LFSR)]) << 8;
                LFSR = ROTR32(LFSR, 8);
                LFSR ^= (((UINT32*)(this->ptable))[(GF)(LFSR)]) << 8;
                LFSR = ROTR32(LFSR, 8);
                LFSR ^= (((UINT32*)(this->ptable))[(GF)(LFSR)]) << 8;
                LFSR = ROTR32(LFSR, 8);
                cnt -= sizeof(UINT32);
                i += sizeof(UINT32);
            }
            *(UINT32*)bb = LFSR;
            while (cnt)
            {
                const GF feedback = bb[0] ^ data[i];
                bb[0] = bb[1] ^ (this->ptable)[4 * feedback];
                bb[1] = bb[2] ^ (this->ptable)[4 * feedback + 1];
                bb[2] = bb[3] ^ (this->ptable)[4 * feedback + 2];
                bb[3] = feedback;
                cnt--;
                i++;
            }
            return;
        }
    }
    else if (1 == this->tt)
    {
        USHORT LFSR = 0;
        int cnt = this->kk;
        int i = 0;
        while (cnt >= 2)
        {
            LFSR ^= *(USHORT*)(data + i);
            LFSR ^= (this->ptable)[2 * (GF)(LFSR)] << 8;
            LFSR ^= (this->ptable)[2 * (GF)(LFSR >> 8)];
            cnt -= 2;
            i += 2;
        }
        *(USHORT*)bb = LFSR;
        if (cnt)
        {
            const GF feedback = bb[0] ^ data[i];
            bb[0] = bb[1] ^ (this->ptable)[2 * feedback];
            bb[1] = feedback;
        }
        return;
    }

    int cnt = 2*tt-1;
    int cDO16 = cnt >> 4;
    int cDO8 = (cnt & 8) >> 3;
    int cDO4 = (cnt & 4) >> 2;
    int cDO2 = (cnt & 2) >> 1;
    int cDO1 = cnt & 1;
    int kk_cnt = kk;
    while (kk_cnt--)
    {
        GF feedback = bb[0]^(*data++);
        GF *TableRow = ptable+2*tt*feedback;

        int j = 0;
        int cnt1 = cDO16;
        while (cnt1--)
        {
            DO8(j); DO8(j+8);
            j += 16;
        }
        if (cDO8)
        {
            DO8(j); j += 8;
        }
        if (cDO4)
        {
            DO4(j); j += 4;
        }
        if (cDO2)
        {
            DO1(j); DO1(j+1);
            j += 2;
        }
        if (cDO1)
        {
            DO1(j); j++;
        }
        DOFB(j);
    }

#else // DEFAULT non-symmetrical, ascending encoder

    for (int i = 0; i < kk; i++)
    {
        GF feedback = bb[0]^data[i];
        GF *TableRow = ptable+2*tt*feedback;
        int j;
        for (j = 0; j < 2*tt-2; j += 2)
        {
            DO1(j); DO1(j+1);
        }
        DO1(j);
        DOFB(j+1);
    }

#endif
#else // ASCENDING_ENCODE: Non-Symmetrical, descending encoders
#if defined(ASSEMBLY_ENCODE) && (tt == 64)

    for (int i = kk-1; i >= 0; i--)
    {
        GF feedback = data[i]^bb[2*tt-1];
        GF *TableRow = ptable+2*tt*feedback;

        DO1(127); DO1(126); DO1(125);

        DO4(121);

        DO8(113);DO8(105);DO8( 97);
        DO8( 89);DO8( 81);DO8( 73);DO8( 65);
        DO8( 57);DO8( 49);DO8( 41);DO8( 33);
        DO8( 25);DO8( 17);DO8(  9);DO8(  1);
        DOFB(0);
    }

#else  // DEFAULT Non-symmetrical, descending encoder

    for (int i = kk-1; i >= 0; i--)
    {
        GF feedback = bb[2*tt-1]^data[i];
        GF *TableRow = ptable+2*tt*feedback;
        for (int j = 2*tt-1; j > 0; j--)
        {
            DO1(j);
        }
        DOFB(0);
    }

#endif // ASSEMBLY_ENCODE
#endif // ASCENDING_ENCODE
#endif // SYMMETRICAL_ENCODE
}


// Special Case Encoders that are selected at run-time.
//
#if !defined(SYMMETRICAL_ENCODE) && defined(ASCENDING_ENCODE) && defined(ASSEMBLY_ENCODE)
void RS_ENCODER::RSEncode64(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 128);

    for (int i = 0; i < kk; i++)
    {
        GF *TableRow;
        GF feedback = bb[0]^data[i];
        TableRow = ptable+128*feedback;

        for (int j = 0; j <= 112; j += 8)
        {
            DO8(j);
        }

        DO4( 120);

        DO1(124); DO1(125); DO1(126);

        DOFB(127);
    }
}

void RS_ENCODER::RSEncode32(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 64);

    for (int i = 0; i < kk; i++)
    {
        GF *TableRow;
        GF feedback = bb[0]^data[i];
        TableRow = ptable+64*feedback;

        for (int j = 0; j <= 48; j += 8)
        {
            DO8(j);
        }

        DO4( 56);

        DO1(60); DO1(61); DO1(62);

        DOFB(63);
    }
}

void RS_ENCODER::RSEncode16(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 32);

    for (int i = 0; i < kk; i++)
    {
        GF *TableRow;
        GF feedback = bb[0]^data[i];
        TableRow = ptable+32*feedback;

        DO8(  0); DO8(  8); DO8( 16);

        DO4( 24);

        DO1(28); DO1(29); DO1(30);

        DOFB(31);
    }
}

void RS_ENCODER::RSEncode14(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 28);

    for (int i = 0; i < kk; i++)
    {
        GF *TableRow;
        GF feedback = bb[0]^data[i];
        TableRow = ptable+28*feedback;

        DO8(  0); DO8(  8); DO8( 16);

        DO1(24); DO1(25); DO1(26);

        DOFB(27);
    }
}

void RS_ENCODER::RSEncode8(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 16);

    for (int i = 0; i < kk; i++)
    {
        GF *TableRow;
        GF feedback = bb[0]^data[i];
        TableRow = ptable+16*feedback;

        DO8(  0);

        DO4(  8);

        DO1(12); DO1(13); DO1(14);

        DOFB(15);
    }
}

void RS_ENCODER::RSEncode4(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    memset(bb, 0, 8);

    for (int i = 0; i < kk; i++)
    {
        GF *TableRow;
        GF feedback = bb[0]^data[i];
        TableRow = ptable+8*feedback;

        DO4(0);

        DO1(4); DO1(5); DO1(6);

        DOFB(7);
     }
}

void RS_ENCODER::RSEncode2(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    UINT32 LFSR = 0;
    int cnt = 62;
    UINT32 *pul = (UINT32 *)data;
    UINT32 *plongtable = (UINT32 *)ptable;
    while (cnt--)
    {
        LFSR ^= *pul++;
        LFSR ^= plongtable[(GF)LFSR] << 8;
        LFSR  = ROTR32(LFSR,8);
        LFSR ^= plongtable[(GF)LFSR] << 8;
        LFSR  = ROTR32(LFSR,8);
        LFSR ^= plongtable[(GF)LFSR] << 8;
        LFSR  = ROTR32(LFSR,8);
        LFSR ^= plongtable[(GF)LFSR] << 8;
        LFSR  = ROTR32(LFSR,8);
    }
    *(UINT32 *)bb = LFSR;
    cnt = 3;
    data += 248;
    while (cnt--)
    {
        GF feedback = bb[0]^(*data++);
        bb[0] = bb[1]^ptable[4*feedback];
        bb[1] = bb[2]^ptable[4*feedback+1];
        bb[2] = bb[3]^ptable[4*feedback+2];
        bb[3] = feedback;
    }
}

void RS_ENCODER::RSEncode1(GF data[MAX_KK], GF bb[2*MAX_TT])
{
    USHORT LFSR = 0;
    USHORT *pus = (USHORT *)data;
    int cnt = 126;
    while (cnt--)
    {
        LFSR ^= *pus++;
#if 0
        LFSR ^= table[(GF)(LFSR)][0] << 8;
        LFSR ^= table[(GF)(LFSR >> 8)][0];
#else
        LFSR ^= ptable[2*(GF)(LFSR)] << 8;
        LFSR ^= ptable[2*(GF)(LFSR >> 8)];
#endif
    }
    *(USHORT *)bb = LFSR;

    GF feedback = bb[0]^data[252];
#if 0
    bb[0] = bb[1]^table[feedback][0];
#else
    bb[0] = bb[1]^ptable[2*feedback];
#endif
    bb[1] = feedback;
}
#endif


//    Performs ERRORS-ONLY decoding of RS codes. If decoding is successful,
// writes the codeword into recd[] itself. Otherwise recd[] is unaltered.  If
// channel caused no more than "tt" errors, the tranmitted codeword will be
// recovered.
//
int RS_ENCODER::RSDecode(GF recd[nn])
{
    RsVerify::verify_received_word(recd, nn);

    // RSDecode variables
    //
    GF root[n];     // must be big enough for all the errors (tt or 2*tt)
    GF err[n];      // must be at least nn
    GF reg[2*MAX_TT+2]; // must be at least tt+1 or 2*tt+1
    GF loc[2*MAX_TT];   // must be at least tt or 2*tt
    GF z[MAX_TT+2];     // must be at least tt+1
    GF elp[2*MAX_TT+2][2*MAX_TT];
    GF d[2*MAX_TT+2];
    GF s[2*MAX_TT+2];   // must be at least tt+1
    int l[2*MAX_TT+2];
    int u_lu[2*MAX_TT+2];

    int i, j, u, q;
    int count=0, syn_error=0;
    int iPowInit, iPow0, iMod;
    //int iPow1;

    // First form the syndromes; i.e., evaluate recd(x) at roots of g(x) namely
    // @^(b0+i), i = 0, ... ,(2*tt-1)
    //
    for (i = 0; i <= 2*tt; i++)
        s[i] = 0;

    iPowInit = 0;
    for (j = 0; j < nn; j++)
    {
        GF RECD = recd[j];
        if (RECD != 0)
        {
            iPow0 = iPowInit + Poly2Pow[RECD];
            MOD_NN(iPow0);
#if defined(ASSEMBLY_DECODE) && ((tt==64) || (tt==32) || (tt==16) || (tt==14) || (tt==8) || (tt==4) || (tt==2) || (tt==1))
            int cnt = 2*j;
            MOD_NN(cnt);
            iPow1 = iPow0 + j;
            MOD_NN(iPow1);

#define DSTEP2(x) { \
    s[x]   ^= Pow2Poly[iPow0]; \
    s[x+1] ^= Pow2Poly[iPow1]; \
    iPow0 += cnt; iPow1 += cnt; \
    MOD_NN(iPow0); MOD_NN(iPow1); \
}
#define DSTEP8(x) { DSTEP2(x) DSTEP2(x+2) DSTEP2(x+4) DSTEP2(x+6) }
#define DSTEP16(x) { DSTEP8(x) DSTEP8(x+8) }
#define DSTEP32(x) { DSTEP16(x) DSTEP16(x+16) }
#define DSTEP64(x) { DSTEP32(x) DSTEP32(x+32) }
#endif
#if defined(ASSEMBLY_DECODE) && (tt == 64)
            DSTEP64(1); DSTEP64(65);
#elif defined(ASSEMBLY_DECODE) && (tt == 32)
            DSTEP64(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 16)
            DSTEP32(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 14)
            DSTEP16(1); DSTEP8(17); DSTEP2(25); DSTEP2(27);
#elif defined(ASSEMBLY_DECODE) && (tt == 8)
            DSTEP16(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 4)
            DSTEP8(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 2)
            DSTEP2(1); DSTEP2(3);
#elif defined(ASSEMBLY_DECODE) && (tt == 1)
            DSTEP2(1);
#else // ASSEMBLY_DECODE
            for (i = 1; i <= 2*tt; i+=2)
            {
                MOD_NN(iPow0);
                s[i]   ^= Pow2Poly[iPow0];
                iPow0 += j;
                MOD_NN(iPow0);
                s[i+1] ^= Pow2Poly[iPow0];
                iPow0 += j;
            }
#endif
        }
        iPowInit += b0;
        MOD_NN(iPowInit);
    }

    for (i = 1; i <= 2*tt; i++)
    {
        // Set flag is non-zero syndrome => error
        //
        if (s[i] != 0)
            syn_error = 1;   /* set flag if non-zero syndrome => error */

        // Convert syndrome from polynomial form to power form
        //
        s[i] = Poly2Pow[s[i]];
    }
    RsVerify::verify_syndromes(s, tt);
    RsDebug::print_syndromes("Initial", s, tt);

#ifdef DECODER_DEBUG
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("The Syndromes of the received vector in polynomial-form are:\n");
    for (i = 1; i <= nn-kk; i++)
    {
        printf("s(%d) = %#x ", i, Pow2Poly[s[i]]);
    }
    printf("\n NOTE: The basis is {@^7,...,@,1}\n");
#endif


#ifdef DECODER_DEBUG
    printf("The Syndromes of the received vector in power-form are:\n");
    for (i = 1; i <= nn-kk; i++)
    {
        if (s[i] == GF_INFINITY)
        {
            printf("s(%d) = 0 ", i);
        }
        else
        {
            printf("s(%d) = @^%d ", i, s[i]);
        }
    }
    printf("\n---------------------------------------------------------------\n");
#endif

    if (syn_error)       /* if errors, try and correct */
    {
#ifdef DECODER_DEBUG
        printf("Received vector is not a codeword. Channel has caused at least one error\n");
#endif

//
// Compute the error location polynomial via the Berlekamp iterative algorithm,
// following the terminology of Lin and Costello :   d[u] is the 'mu'th
// discrepancy, where u = 'mu'+ 1 and 'mu' (the Greek letter!) is the step
// number ranging from -1 to 2*tt (see L&C),  l[u] is the degree of the elp at
// that step, and u_l[u] is the difference between the
// step number and the degree of the elp.
//
// The notation is the same as that in Lin and Costello's book; pages 155-156
// and 175.
//
/* initialise table entries */
        d[0] = 0;           /* power form */
        d[1] = s[1];        /* power form */
        elp[0][0] = 0;      /* power form */
        elp[1][0] = 1;      /* polynomial form */
        for (i = 1; i < nn-kk; i++)
        {
            elp[0][i] = GF_INFINITY;   /* power form */
            elp[1][i] = 0;   /* polynomial form */
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;

        do
        {
            u++;
            if (d[u]==GF_INFINITY)
            {
                l[u+1] = l[u];
                for (i=0; i<=l[u]; i++)
                {
                    elp[u+1][i] = elp[u][i];
                    elp[u][i] = Poly2Pow[elp[u][i]];
                }
            }
            else
            /* search for words with greatest u_lu[q] for which d[q]!=0 */
            {
                q = u-1;
                while ((d[q] == GF_INFINITY) && (q>0))
                    q--;

                /* have found first non-zero d[q]  */
                if (q>0)
                {
                    j=q;
                    do
                    {
                        j--;
                        if ((d[j] != GF_INFINITY) && (u_lu[q]<u_lu[j]))
                            q = j;
                    } while (j>0);
                }

                /* have now found q such that d[u]!=0 and u_lu[q] is maximum */
                /* store degree of new elp polynomial */
                if (l[u] > l[q]+u-q)
                    l[u+1] = l[u];
                else
                    l[u+1] = l[q]+u-q;

                /* form new elp(x) */
                for (i = 0; i < nn-kk; i++)
                    elp[u+1][i] = 0;
                for (i=0; i<=l[q]; i++)
                    if (elp[q][i] != GF_INFINITY)
                    {
                        iMod = d[u]-d[q]+elp[q][i];
                        MOD_NN(iMod);
                        elp[u+1][i+u-q] = Pow2Poly[iMod];
                    }
                for (i=0; i<=l[u]; i++)
                {
                    elp[u+1][i] ^= elp[u][i];
                    elp[u][i] = Poly2Pow[elp[u][i]];  /*convert old elp value to power*/
                }
            }
            u_lu[u+1] = u-l[u+1];

            /* form (u+1)th discrepancy */
            if (u<nn-kk)    /* no discrepancy computed on last iteration */
            {
                d[u+1] = Pow2Poly[s[u+1]];
                for (i = 1; i <= l[u+1]; i++)
                {
                    if ((s[u+1-i] != GF_INFINITY) && (elp[u+1][i] != 0))
                    {
                        iMod = s[u+1-i]+Poly2Pow[elp[u+1][i]];
                        MOD_NN(iMod);
                        d[u+1] ^= Pow2Poly[iMod];
                    }
                }
                d[u+1] = Poly2Pow[d[u+1]];    /* put d[u+1] into power form */
            }
            RsDebug::print_berlekamp_step(u, d[u], l[u], elp[u], l[u]);
        } while ((u < nn-kk) && (l[u+1] <= tt));

        u++;
        RsVerify::verify_lambda(elp[u], l[u]);

        if (l[u] <= tt)         /* can correct error */
        {
#ifdef DECODER_DEBUG
            printf("---------------------------------------------------------------\n");
            printf("The Final Error Locator Polynomial Lambda(x) in polynomial-form is: \n");
            for (i=0;i < l[u];i++)
            {
                printf("%#x x^%d + ",elp[u][i],i);
                if (i && (i % 6 == 0))
                    printf("\n"); /* 7 coefficients per line */
            }
            printf("%#x x^%d",elp[u][i],i);
            printf("\n");
#endif

            /* put elp into power form */
            for (i=0; i<=l[u]; i++)
                elp[u][i] = Poly2Pow[elp[u][i]];

#ifdef DECODER_DEBUG
            printf("The Final Error Locator Polynomial Lambda(x) in power-form is: \n");
            for (i=0;i < l[u];i++)
            {
                printf("%d x^%d + ",elp[u][i],i);
                if (i && (i % 6 == 0))
                    printf("\n"); /* 7 coefficients per line */
            }
            printf("%d x^%d",elp[u][i],i);
            printf("\n");
            printf("---------------------------------------------------------------\n");
#endif


            /* find roots of the error location polynomial */
            for (i=1; i <= l[u]; i++)
                reg[i] = elp[u][i];
            count = 0;
            for (i=1; i <= nn; i++)
            {
                q = 1;
                for (j = 1; j <= l[u]; j++)
                {
                    if (reg[j] != GF_INFINITY)
                    {
                        iMod = reg[j]+j;
                        MOD_NN(iMod);
                        q ^= Pow2Poly[iMod];
                        reg[j] = iMod;
                    }
                }
                if (!q)  /* store root and error location number indices */
                {
                    root[count] = i;
                    loc[count] = nn-i;
                    count++;
                }
            }
#ifndef NO_PRINT
            printf("Computed Error Locations\n");
            for (i=0;i < l[u];i++)
                printf("%d ", loc[i]);
            printf("\n");
            printf("(count, l[u]) = (%d, %d)\n", count, l[u]);
#endif
            if (count == l[u])    /* no. roots = degree of elp hence <= tt errors */
            {
#ifdef DECODER_DEBUG
                printf("\n No of roots of Lambda(x) found by Chien search equals with \n");
                printf("the degree of Lambda(x). Hence channel presumably caused <= %d errors\n",tt);
                printf("The roots found are (in power-form) : \n");
                for (i=0;i < count;i++)
                    printf("%d ",root[i]);
                printf("\n");
#endif
                /* form polynomial z(x) */
                for (i = 1; i <= l[u]; i++) /* Z[0] = 1 always - do not need */
                {
                    GF zi = Pow2Poly[s[i]] ^ Pow2Poly[elp[u][i]];
                    for (j=1; j<i; j++)
                    {
                        if ((s[j] != GF_INFINITY) && (elp[u][i-j] != GF_INFINITY))
                        {
                            iMod = elp[u][i-j] + s[j];
                            MOD_NN(iMod);
                            zi ^= Pow2Poly[iMod];
                        }
                    }
                    z[i] = Poly2Pow[zi];         /* put into power form */
                }

#ifdef DECODER_DEBUG
                printf("---------------------------------------------------------------\n");
                printf(" \n The Final Error Evaluator Polynomial Omega(x) in power-form is: \n");
                printf("0 x^0 + ");
                for (i=1;i <= l[u];i++)
                {
                    printf("%d x^%d + ",z[i],i);
                    if (i && (i % 6 == 0))
                        printf("\n"); /* 7 coefficients per line */
                }
                printf("%d x^%d",z[i],i);
                printf("\n");
#endif

#ifdef DECODER_DEBUG
                printf(" \n The Final Error Evaluator Polynomial Omega(x) in polynomial-form is: \n");
                printf("%#x x^%d + ",1,0);
                for (i=1;i <= l[u];i++)
                {
                    printf("%#x x^%d + ",z[i],i);
                    if (i && (i % 6 == 0))
                        printf("\n"); /* 7 coefficients per line */
                }
                printf("%#x x^%d",z[i],i);
                printf("\n");
                printf("---------------------------------------------------------------\n");
#endif

                // Evaluate errors at locations given by error location
                // numbers loc[i]
                //
                for (i = 0; i < nn; i++)
                    err[i] = 0;

                // Compute numerator of error term first
                //
                for (i=0; i < l[u]; i++)
                {
                    int jPow = root[i];
                    err[loc[i]] = 1;       /* accounts for z[0] */
                    for (j=1; j<=l[u]; j++)
                    {
                        if (z[j] != GF_INFINITY)
                        {
                            iMod = z[j]+jPow;
                            MOD_NN(iMod);
                            err[loc[i]] ^= Pow2Poly[iMod];
                        }
                        jPow += root[i]; // jPow = root[i]*j;
                        MOD_NN(jPow);
                    } /* z(x) evaluated at X(l)^(-1) */

                    // term X(l)^(1-b0)
                    //
                    if (err[loc[i]] != 0)
                    {
                        err[loc[i]] = Pow2Poly[(Poly2Pow[err[loc[i]]]+root[i]*(b0+nn-1))%nn];
                    }
                    if (err[loc[i]] != 0)
                    {
                        err[loc[i]] = Poly2Pow[err[loc[i]]];
                        q = 0;     /* form denominator of error term */
                        for (j = 0; j < l[u]; j++)
                            if (j != i)
                            {
                                iMod = loc[j]+root[i];
                                MOD_NN(iMod);
                                q += Poly2Pow[1^Pow2Poly[iMod]];
                                MOD_NN(q);
                            }
                        // q = q % nn;
                        iMod = err[loc[i]]-q;
                        MOD_NN(iMod);
                        err[loc[i]] = Pow2Poly[iMod];
                        recd[loc[i]] ^= err[loc[i]];  /*recd[i] must be in polynomial form */
                        RsDebug::print_error_location(loc[i], err[loc[i]], recd[loc[i]] ^ err[loc[i]], recd[loc[i]]);
                    }
                }
#ifdef DECODER_DEBUG
                printf("---------------------------------------------------------------\n");
                printf("The computed Error Value Bytes are (in polynomial-form) :\n");
                printf("Location (Error) \n");
                for (i=0;i < l[u];i++)
                {
                    printf("%d (%#x) ", loc[i], err[loc[i]]);
                }
                printf("\n");
                printf("---------------------------------------------------------------\n");
#endif
                return count; // Number of corrections.
            }
            else
            {    /* no. roots != degree of elp => >tt errors and cannot solve */
#ifdef DECODER_DEBUG
                printf("No of roots of Lambda(x) found by Chien search\n");
                printf("does not equal the degree of the Lambda(x).\n");
                printf("Hence channel has definitely caused >= %d errors\n",tt);
#endif
                return -2; // NOT ENOUGH ROOTS BY CHIEN SEARCH
            }
        }
        else
        {         /* elp has degree has degree >tt hence cannot solve */
#ifdef DECODER_DEBUG
            printf("Degree of the Lambda(x) > %d. \n",tt);
            printf("Hence channel has definitely caused >= %d errors\n",tt);
#endif
            return -3; // LAMDA DEGREE TOO HIGH
        }
#ifdef DECODER_DEBUG
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
    }
    else
    {
        /* no non-zero syndromes => no errors: output received codeword */
#ifdef DECODER_DEBUG
        printf("Received vector is a codeword. Channel has presumably caused no error\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
        return 0;
    }
}

// Performs ERRORS+ERASURES decoding of RS codes. If decoding is successful,
// writes the codeword into recd[] itself. Otherwise recd[] is unaltered
// except in the erased positions where it is set to zero. First "no_eras"
// erasures are declared by the calling program. Then, the maximum # of errors
// correctable is t_after_eras = floor((2*tt-no_eras)/2). If the number of
// channel errors is not greater than "t_after_eras" the transmitted codeword
// will be recovered. Details of algorithm can be found in R. Blahut's
// "Theory ... of Error-Correcting Codes".
//
int RS_ENCODER::RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras)
{
    // RSDecodeErasures variables
    //
    GF root[n];     // must be big enough for all the errors (tt or 2*tt)
    GF err[n];      // must be at least nn
    GF reg[2*MAX_TT+2]; // must be at least tt+1 or 2*tt+1
    GF loc[2*MAX_TT];   // must be at least tt or 2*tt
    GF s[2*MAX_TT+2];   // must be at least tt+1

    // RSDecodeErasures variables
    //
    GF phi[2*MAX_TT+1], tmp_pol[2*MAX_TT+1];
    GF lambda[2*MAX_TT+1];
    GF lambda_pr[2*MAX_TT+1];
    GF omega[2*MAX_TT+1];
    GF b[2*MAX_TT+1], T[2*MAX_TT+1];
    GF pres_root, pres_loc;
    GF num1, num2;
    GF q, den;
    GF discr_r;

    GF U;
    int deg_phi, deg_lambda = 0, L, deg_omega;
    int syn_error = 0;
    int i, j, r;
    int iMod;
    int count;

    // Maximum # ch errs correctable after "no_eras" erasures.
    //
    //int t_after_eras = (2*tt - no_eras)/2;

    // Compute erasure locator polynomial phi[x].
    //
    for (i=0;i < nn-kk+1;i++)
        tmp_pol[i] = 0;
    for (i=0;i < nn-kk+1;i++)
        phi[i] = 0;
    if (no_eras > 0)
    {
        phi[0] = 1; /* power form */
        phi[1] = Pow2Poly[eras_pos[0]];
        for (i=1;i < no_eras;i++)
        {
            U = eras_pos[i];
            for (j=1;j < i+2;j++)
            {
                if (phi[j-1] == 0)
                {
                    tmp_pol[j] = 0;
                }
                else
                {
                    iMod = U+Poly2Pow[phi[j-1]];
                    MOD_NN(iMod);
                    tmp_pol[j] = Pow2Poly[iMod];
                }
            }
            for (j=1;j < i+2;j++)
                phi[j] = phi[j]^tmp_pol[j];
        }
        /* put phi[x] in power form */
        for (i=0;i < nn-kk+1;i++)
            phi[i] = Poly2Pow[phi[i]];
#ifdef ERASURE_DEBUG
        /* find roots of the erasure location polynomial */
        for (i=1; i <= no_eras; i++)
            reg[i] = phi[i];
        count = 0;
        for (i=1; i <= nn; i++)
        {
            q = 1;
            for (j=1; j <= no_eras; j++)
            {
                if (reg[j] != GF_INFINITY)
                {
                    iMod = reg[j]+j;
                    MOD_NN(iMod);
                    q ^= Pow2Poly[iMod];
                    reg[j] = iMod;
                }
            }
            if (!q)        /* store root and error location number indices */
            {
                root[count] = i;
                loc[count] = nn-i;
                count++;
            }
        }
        if (count != no_eras)
        {
            printf("\n phi(x) is WRONG\n"); exit(1);
        }
#ifndef NO_PRINT
        printf("\n Erasure positions as determined by roots of Eras Loc Poly:\n");
        for (i = 0; i < count; i++)
            printf("%d ",loc[i]);
        printf("\n");
#endif
#endif
    }

    // First form the syndromes: i.e., evaluate recd(x) at roots of g(x) namely
    // @^(b0+i), i = 0, ... ,(2*tt-1)
    //
    for (i = 0; i <= 2*tt; i++)
        s[i] = 0;

    int iPowInit = 0, iPow0;
    //int iPow1;
    for (j = 0; j < nn; j++)
    {
        GF RECD = recd[j];
        if (RECD != 0)
        {
            iPow0 = iPowInit + Poly2Pow[RECD];
            MOD_NN(iPow0);
#if defined(ASSEMBLY_DECODE) && ((tt==64) || (tt==32) || (tt==16) || (tt==14) || (tt==8) || (tt==4) || (tt==2) || (tt==1))
            int cnt = 2*j;
            MOD_NN(cnt);
            iPow1 = iPow0 + j;
            MOD_NN(iPow1);

#define DSTEP2(x) { \
    s[x]   ^= Pow2Poly[iPow0]; \
    s[x+1] ^= Pow2Poly[iPow1]; \
    iPow0 += cnt; iPow1 += cnt; \
    MOD_NN(iPow0); MOD_NN(iPow1); \
}
#define DSTEP8(x) { DSTEP2(x) DSTEP2(x+2) DSTEP2(x+4) DSTEP2(x+6) }
#define DSTEP16(x) { DSTEP8(x) DSTEP8(x+8) }
#define DSTEP32(x) { DSTEP16(x) DSTEP16(x+16) }
#define DSTEP64(x) { DSTEP32(x) DSTEP32(x+32) }
#endif
#if defined(ASSEMBLY_DECODE) && (tt == 64)
            DSTEP64(1); DSTEP64(65);
#elif defined(ASSEMBLY_DECODE) && (tt == 32)
            DSTEP64(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 16)
            DSTEP32(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 14)
            DSTEP16(1); DSTEP8(17); DSTEP2(25); DSTEP2(27);
#elif defined(ASSEMBLY_DECODE) && (tt == 8)
            DSTEP16(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 4)
            DSTEP8(1);
#elif defined(ASSEMBLY_DECODE) && (tt == 2)
            DSTEP2(1); DSTEP2(3);
#elif defined(ASSEMBLY_DECODE) && (tt == 1)
            DSTEP2(1);
#else // ASSEMBLY_DECODE
            for (i = 1; i <= 2*tt; i+=2)
            {
                MOD_NN(iPow0);
                s[i]   ^= Pow2Poly[iPow0];
                iPow0 += j;
                MOD_NN(iPow0);
                s[i+1] ^= Pow2Poly[iPow0];
                iPow0 += j;
            }
#endif
        }
        iPowInit += b0;
        MOD_NN(iPowInit);
    }

    for (i = 1; i <= 2*tt; i++)
    {
        // Set flag is non-zero syndrome => error
        //
        if (s[i] != 0)
            syn_error = 1;   /* set flag if non-zero syndrome => error */

        // Convert syndrome from polynomial form to power form
        //
        s[i] = Poly2Pow[s[i]];
    }


    // If syndrome is zero, modified recd[] is a codeword, otherwise we perform
    // Berlekamp-Massey procedure to determine location of errors.
    //
    if (syn_error)
    {
        // Begin Berlekamp-Massey algorithm to determine error+erasure locator
        // polynomial
        r = no_eras;
        deg_phi = no_eras;
        L = no_eras;
        if (no_eras > 0)
        {
            /* Initialize lambda(x) and b(x) (in poly-form) to phi(x) */
            for (i=0;i < deg_phi+1; i++)
                lambda[i] = Pow2Poly[phi[i]];
            for (i=deg_phi+1;i < 2*tt+1;i++)
                lambda[i] = 0;
            deg_lambda = deg_phi;
            for (i=0; i < 2*tt+1; i++)
                b[i] = lambda[i];
        }
        else
        {
            lambda[0] = 1;
            for (i=1;i < 2*tt+1;i++)
                lambda[i] = 0;
            for (i=0;i < 2*tt+1;i++)
                b[i] = lambda[i];
        }
        while (++r <= 2*tt)
        {
            /* r is the step number */
            /* Compute discrepancy at the r-th step in poly-form */
            discr_r = 0;
            for (i=0;i < 2*tt+1;i++)
            {
                if ((lambda[i] != 0) && (s[r-i] != GF_INFINITY))
                {
                    iMod = Poly2Pow[lambda[i]]+s[r-i];
                    MOD_NN(iMod);
                    discr_r ^= Pow2Poly[iMod];
                }
            }
            if (discr_r == 0)
            {
                // 3 lines below: B(x) <-- x*B(x)
                //
                for (i = 2*tt; i > 0; i--)
                {
                    b[i] = b[i-1];
                }
                b[0] = 0;
            }
            else
            {
                /* 5 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
                T[0] = lambda[0];
                for (i=1;i < 2*tt+1;i++)
                {
                    if (b[i-1] == 0)
                    {
                        T[i] = lambda[i];
                    }
                    else
                    {
                        iMod = Poly2Pow[discr_r]+Poly2Pow[b[i-1]];
                        MOD_NN(iMod);
                        T[i] = lambda[i]^Pow2Poly[iMod];
                    }
                }

                if (2*L <= r+no_eras-1)
                {
                    L = r+no_eras-L;
                    /* 2 lines below: B(x) <-- inv(discr_r) * lambda(x) */
                    for (i = 0; i < 2*tt+1;i++)
                    {
                        if (lambda[i] == 0)
                        {
                            b[i] = 0;
                        }
                        else
                        {
                            iMod = Poly2Pow[lambda[i]]-Poly2Pow[discr_r];
                            MOD_NN(iMod);
                            b[i] = Pow2Poly[iMod];
                        }
                    }
                    for (i=0;i < 2*tt+1;i++)
                        lambda[i] = T[i];
                }
                else
                {
                    for (i=0;i < 2*tt+1;i++)
                        lambda[i] = T[i];
                    // 3 lines below: B(x) <-- x*B(x)
                    //
                    for (i = 2*tt; i > 0; i--)
                    {
                        b[i] = b[i-1];
                    }
                    b[0] = 0;
                }
            }
        }
        /* Put lambda(x) into power form */
        for (i=0; i < 2*tt+1; i++)
            lambda[i] = Poly2Pow[lambda[i]];
        RsVerify::verify_lambda(lambda, deg_lambda);

        /* Compute deg(lambda(x)) */
        deg_lambda = 2*tt;
        while ((lambda[deg_lambda] == GF_INFINITY) && (deg_lambda > 0))
            --deg_lambda;

        if (deg_lambda <= 2*tt)
        {
            // Find roots of the error+erasure locator polynomial. By
            // Chien Search
            //
            for (i=1; i < 2*tt+1; i++)
                reg[i] = lambda[i];
            count = 0; /* Number of roots of lambda(x) */
            for (i=1; i <= nn; i++)
            {
                q = 1;
                for (j=1; j <= deg_lambda; j++)
                    if (reg[j] != GF_INFINITY)
                    {
                        iMod = reg[j]+j;
                        MOD_NN(iMod);
                        q ^= Pow2Poly[iMod];
                        reg[j] = iMod;
                    }
                if (!q) /* store root (power-form) and error location number */
                {
                    root[count] = i;
                    loc[count] = nn-i;
                    count++;
                }
            }
#ifndef NO_PRINT
            printf("\n Final error positions:\t");
            for (i=0;i < count;i++)
                printf("%d ",loc[i]);
            printf("\n");
#endif
            // BUG!: If we didn't find as many roots as we expected, then
            // some of the erasures were false indicators. If deg_lamba > count
            // we may need to recompute lambda
            //
            if (deg_lambda == count)
            {
                // This is a correctable error
                // Compute err+eras evaluator poly in polynomial form.
                //    omega(x) = s(x)*lambda(x) (modulo x^(nn-kk))
                //
                for (i=0;i < 2*tt;i++)
                {
                    omega[i] = 0;
                    for (j=0;(j < deg_lambda+1) && (j < i+1);j++)
                    {
                        if ((s[i+1-j] != GF_INFINITY) && (lambda[j] != GF_INFINITY))
                        {
                            iMod = s[i+1-j]+lambda[j];
                            MOD_NN(iMod);
                            omega[i] ^= Pow2Poly[iMod];
                        }
                    }
                }
                omega[2*tt] = 0;

                // Compute lambda_pr(x) = formal derivative of lambda(x) in
                // poly-form
                //
                for (i = 0; i < 2*tt; i += 2)
                {
                    lambda_pr[i+1] = 0;
                    lambda_pr[i] = Pow2Poly[lambda[i+1]];
                }
                lambda_pr[2*tt] = 0;

                /* Compute deg(omega(x)) */
                deg_omega = 2*tt;
                while ((omega[deg_omega] == 0) && (deg_omega > 0))
                    --deg_omega;

                // Compute error values in polynomial form.
                // num1 = omega(inv(X(l))),
                // num2 = inv(X(l))^(b0-1) and den = lambda_pr(inv(X(l))) all
                // in polynomial form.
                //
                for (j=0;j < count;j++)
                {
                    pres_root = root[j];
                    pres_loc = loc[j];
                    num1 = 0;
                    iPow0 = 0;
                    for (i=0;i < deg_omega+1;i++)
                    {
                        if (omega[i] != 0)
                        {
                            iMod = Poly2Pow[omega[i]]+iPow0;
                            MOD_NN(iMod);
                            num1 ^= Pow2Poly[iMod];
                        }
                        iPow0 += pres_root;
                        MOD_NN(iPow0);
                    }
                    iMod = b0-1;
                    MOD_NN(iMod);
                    iMod = (iMod*pres_root)%nn;
                    num2 = Pow2Poly[iMod];
                    den = 0;
                    iPow0 = 0;
                    for (i = 0; i < deg_lambda+1; i++)
                    {
                        if (lambda_pr[i] != 0)
                        {
                            iMod = Poly2Pow[lambda_pr[i]]+iPow0;
                            MOD_NN(iMod);
                            den ^= Pow2Poly[iMod];
                        }
                        iPow0 += pres_root;
                        MOD_NN(iPow0);
                    }
                    if (den == 0)
                    {
                        //printf("\n ERROR: denominator = 0\n");
                    }
                    err[pres_loc] = 0;
                    if (num1 != 0)
                    {
                        iMod = Poly2Pow[num1]+Poly2Pow[num2]-Poly2Pow[den];
                        MOD_NN(iMod);
                        err[pres_loc] = Pow2Poly[iMod];
                    }
                }
                // Correct word by subtracting out error bytes.
                //
                for (j=0;j < count;j++)
                    recd[loc[j]] ^= err[loc[j]];

                return count;
            }
            else
            {
                /* deg(lambda) unequal to number of roots => uncorrectable error detected */
                //printf("deg_lambda(%d) != #root(%d)\n", deg_lambda, count);
                return -2;
            }
        }
        else
            /* deg(lambda) > 2*tt => uncorrectable error detected */
            return -3;
    }
    else
    {
        /* no non-zero syndromes => no errors: output received codeword */
        return 0;
    }
}

bool RS_ENCODER::bInitialized = false;

void RS_ENCODER::Init(void)
{
    if (!bInitialized)
    {
        bInitialized = true;

        // Create Modulus Table
        //
	pModTable = ModTable+nn;
        for (int i = -nn; i <= 3*nn; i++)
        {
            pModTable[i] = (nn+i)%nn;
        }
    }
}
