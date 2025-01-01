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
#include <iostream>
#include "projdefs.h"
#include "rs.h"
#include "RsVerification.h"

#define NO_PRINT

#define n   256         // n  = 2^8     size of the field
#define nn  255         // nn = 2^8-1   length of codeword
#define MAX_TT (127)    // maximum number of correctable errors (nn/2).
#define MIN_TT (1)      // minimum number of correctable errors (1).
#define MAX_KK (253)    // Maximum number of data symbols (nn-2*TT_MIN).
#define MIN_KK (1)      // Minimum number of data symbols (nn-2*TT_MAX).

GF RS_ENCODER::ModTable[4*nn+1];
GF *RS_ENCODER::pModTable;
#define MOD_NN(x) (pModTable[(x)])

RS_ENCODER::RS_ENCODER(int CorrectableErrors) : tt(CorrectableErrors), kk(nn - 2 * CorrectableErrors) {
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
    b0 = (kk + 1) / 2;  // Symmetric generator polynomial
#else
    b0 = 254;
#endif
    RSGenPoly();
    RSGenTable();
    RsVerification::verify_generator(gg, tt);
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
                iMod = MOD_NN(iMod);
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
        iMod = MOD_NN(iMod);
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
    int ch, iMod;
    ptable = new GF[256*SYM_GEN_ENC_LEN];
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
        for (j = 0; j < tt; j++)
        {
            if (gg[j+1] != GF_INFINITY && feedback != GF_INFINITY)
            {
                iMod = gg[j+1]+feedback;
                iMod = MOD_NN(iMod);
                ptable[ch*SYM_GEN_ENC_LEN+j] = Pow2Poly[iMod];
            }
            else
            {
                ptable[ch*SYM_GEN_ENC_LEN+j] = 0;
            }
        }
    }
#else
    // Could be 2*tt-1, but we may want a little bigger for alignment.
    //
    int ch, iMod;
    ptable = new GF[256*SYM_GEN_LEN];
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
        for (j = 0; j < 2*tt; j++)
        {
            if (gg[j+1] != GF_INFINITY && feedback != GF_INFINITY)
            {
                iMod = gg[j+1]+feedback;
                iMod = MOD_NN(iMod);
                ptable[ch*SYM_GEN_LEN+j] = Pow2Poly[iMod];
            }
            else
            {
                ptable[ch*SYM_GEN_LEN+j] = 0;
            }
        }
    }
#endif
#else
    // Could be 2*tt+1, but we may want a little bigger for alignment.
    //
    int ch, iMod;
    ptable = new GF[256*NON_SYM_LEN];
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
        for (j = 0; j <= 2*tt; j++)
        {
            if (gg[j] != GF_INFINITY && feedback != GF_INFINITY)
            {
                iMod = gg[j]+feedback;
                iMod = MOD_NN(iMod);
                ptable[ch*NON_SYM_LEN+j] = Pow2Poly[iMod];
            }
            else
            {
                ptable[ch*NON_SYM_LEN+j] = 0;
            }
        }
    }
#endif
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
        GF *TableRow = ptable+SYM_GEN_ENC_LEN*feedback;
        for (int j = 0; j < tt; j++)
        {
            bb[j] = bb[j+1]^TableRow[j];
        }
        for (int j = tt; j < 2*tt-1; j++)
        {
            bb[j] = bb[j+1]^TableRow[2*tt-j-2];
        }
        bb[2*tt-1] = feedback;
    }

#else // ASCENDING_ENCODE: Symmetrical, descending encoder

    for (int i = kk-1; i >= 0; i--)
    {
        GF feedback = bb[2*tt-1]^data[i];
        GF *TableRow = ptable+SYM_GEN_ENC_LEN*feedback;
        for (int j = 2*tt-1; j > tt; j--)
        {
            bb[j] = bb[j-1]^TableRow[2*tt-j-1];
        }
        for (int j = tt; j > 0; j--)
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

int RS_ENCODER::RSDecode(GF recd[nn]) {
    // Verify received word
    RsVerification::verify_received_word(recd, nn);

    // Step 1: Calculate syndromes
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    // Check for errors
    bool has_error = false;
    for (int i = 1; i <= 2 * tt; i++) {
        if (syndromes[i] != 0) {
            has_error = true;
        }
    }

    RsVerification::verify_syndromes(syndromes, tt);
    RsVerification::print_syndromes("Initial", syndromes, tt);

    if (!has_error) {
        return 0; // No errors detected
    }

    // Step 2: Compute error locator polynomial using Berlekamp-Massey
    std::vector<GF> lambda(2 * tt + 1, 0);
    lambda[0] = 1; // Start with identity polynomial
    berlekamp_massey(syndromes, lambda, 0);
    int deg_lambda = convert_to_index_and_get_degree(lambda);
    if (deg_lambda > 2 * tt)
    {
        return RS_ERROR_LAMBDA_ERROR;
    }

    // Verify error locator polynomial
    RsVerification::verify_lambda(lambda, deg_lambda);

    // Step 3: Find roots of the error locator polynomial using Chien Search
    std::vector<GF> root(2 * tt);
    std::vector<GF> loc(2 * tt);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);
    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    // Step 4: Compute evaluator polynomial omega
    std::vector<GF> omega(2 * tt + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    // Step 5: Apply Forney’s Algorithm for error correction
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0)
    {
        return -1;
    }
    return count;
}

void RS_ENCODER::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt + 1, 0);

    int iPowInit = 0;
    for (int j = 0; j < nn; j++) {
        GF RECD = recd[j];
        if (RECD != 0) {
            int iPow0 = iPowInit + Poly2Pow[RECD];
            iPow0 = MOD_NN(iPow0);

            for (int i = 1; i <= 2 * tt; i++) {
                iPow0 = MOD_NN(iPow0);
                syndromes[i] ^= Pow2Poly[iPow0];
                iPow0 += j;
            }
        }
        iPowInit += b0;
        iPowInit = MOD_NN(iPowInit);
    }

    for (int i = 1; i <= 2 * tt; i++) {
        syndromes[i] = Poly2Pow[syndromes[i]];
    }
}

static constexpr GF A0 = nn;

int RS_ENCODER::construct_erasure_locator(std::vector<GF>& lambda, const int* eras_pos, int no_eras) {
    if (no_eras == 0) {
        lambda.assign(lambda.size(), 0);
        lambda[0] = 1;
        return 0; // No erasures
    }

    lambda.assign(lambda.size(), 0);
    lambda[0] = 1;

    for (int i = 0; i < no_eras; i++) {
        GF u = eras_pos[i];
        for (int j = no_eras; j > 0; j--) {
            if (lambda[j - 1] != A0) {
                lambda[j] ^= Pow2Poly[MOD_NN(u + Poly2Pow[lambda[j - 1]])];
            }
        }
    }

    return no_eras; // Degree of the erasure locator polynomial
}

void RS_ENCODER::berlekamp_massey(const std::vector<GF>& syndromes, std::vector<GF>& lambda, int no_eras) {
    std::vector<GF> b(2 * tt + 1, 0);
    for (int i = 0; i < 2 * tt + 1; i++)
        b[i] = Poly2Pow[lambda[i]];

    std::vector<GF> t(2 * tt + 1, 0);

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    int r = no_eras;
    int el = no_eras;
    while (++r <= 2*tt) { /* r is the step number */
        /* Compute discrepancy at the r-th step in poly-form */
        int discr_r = 0;
        for (int i = 0; i < r; i++){
            if ((lambda[i] != 0) && (syndromes[r - i] != A0)) {
                discr_r ^= Pow2Poly[MOD_NN(Poly2Pow[lambda[i]] + syndromes[r - i])];
            }
        }
        RsVerification::print_berlekamp_step(r, discr_r, r, lambda, r);

        discr_r = Poly2Pow[discr_r]; /* Index form */
        if (discr_r == A0) {
            /* 2 lines below: B(x) <-- x*B(x) */
            std::copy_backward(b.begin(), b.begin() + (2 * tt), b.begin() + (2 * tt + 1));
            b[0] = A0;
        } else {
            /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
            t[0] = lambda[0];
            for (int i = 0 ; i < 2*tt; i++) {
                if(b[i] != A0)
                    t[i+1] = lambda[i+1] ^ Pow2Poly[MOD_NN(discr_r + b[i])];
                else
                    t[i+1] = lambda[i+1];
            }
            if (2 * el <= r + no_eras - 1) {
                el = r + no_eras - el;
                /*
                 * 2 lines below: B(x) <-- inv(discr_r) *
                 * lambda(x)
                 */
                for (int i = 0; i <= 2*tt; i++)
                    b[i] = (lambda[i] == 0) ? A0 : MOD_NN(Poly2Pow[lambda[i]] - discr_r + nn);
            } else {
                /* 2 lines below: B(x) <-- x*B(x) */
                std::copy_backward(b.begin(), b.begin() + (2 * tt), b.begin() + (2 * tt + 1));
                b[0] = A0;
            }
            std::copy(t.begin(), t.begin() + (2 * tt + 1), lambda.begin());
        }
    }
}

int RS_ENCODER::convert_to_index_and_get_degree(std::vector<GF>& poly) {
    int degree = 0;
    for (size_t i = 0; i < poly.size(); i++) {
        poly[i] = Poly2Pow[poly[i]];
        if (poly[i] != A0)
            degree = static_cast<int>(i);
    }
    return degree;
}

int RS_ENCODER::chien_search(const std::vector<GF>& lambda, int deg_lambda, std::vector<GF>& root, std::vector<GF>& loc, int& count) {
    std::vector<GF> reg(2 * tt + 1);
    memcpy(&reg[1], &lambda[1], 2 * tt * sizeof(GF));

    count = 0;
    for (int i = 1; i <= nn; i++) {
        GF q = 1;
        for (int j = deg_lambda; j > 0; j--) {
            if (reg[j] != A0) {
                reg[j] = MOD_NN(reg[j] + j);
                q ^= Pow2Poly[reg[j]];
            }
        }
        if (q == 0) {
            root[count] = i;
            loc[count] = nn - i;
            count++;
        }
    }

    return count;
}

int RS_ENCODER::compute_omega(const std::vector<GF>& syndromes, const std::vector<GF>& lambda, int deg_lambda, std::vector<GF>& omega) {
    int deg_omega = 0;
    for (int i = 0; i < 2 * tt; i++) {
        GF tmp = 0;
        int j = (deg_lambda < i) ? deg_lambda : i;
        for ( ; j >= 0; j--) {
            if ((syndromes[i + 1 - j] != A0) && (lambda[j] != A0))
                tmp ^= Pow2Poly[MOD_NN(syndromes[i + 1 - j] + lambda[j])];
        }
        if (tmp != 0)
            deg_omega = i;
        omega[i] = Poly2Pow[tmp];
    }
    omega[2*tt] = A0;
    return deg_omega;
}

// Forney’s Algorithm: Compute error magnitudes and correct errors
int RS_ENCODER::forney_correction(const std::vector<GF>& omega, int deg_omega, const std::vector<GF>& lambda, int deg_lambda,
                                       const std::vector<GF>& root, int count, const std::vector<GF>& loc, GF data[nn]) {
    for (int j = count - 1; j >= 0; j--) {
        GF num1 = 0;
        for (int i = deg_omega; i >= 0; i--) {
            if (omega[i] != A0)
                num1 ^= Pow2Poly[MOD_NN(omega[i] + i * root[j])];
        }
        GF num2 = Pow2Poly[MOD_NN(root[j] * (b0 - 1) + nn)];
        GF den = 0;

        /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
        for (int i = std::min(deg_lambda, 2 * tt - 1) & ~1; i >= 0; i -=2) {
            if (lambda[i + 1] != A0)
                den ^= Pow2Poly[MOD_NN(lambda[i+1] + i * root[j])];
        }

        if (den == 0) {
            return -1;
        }

        /* Apply error to data */
        if (num1 != 0) {
            data[loc[j]] ^= Pow2Poly[MOD_NN(Poly2Pow[num1] + Poly2Pow[num2] + nn - Poly2Pow[den])];
        }
    }
    return 0;
}

// RSDecodeErasures: Port of Phil Karn's Berlekamp-Massey implementation
// Brain-dead port for erasure decoding with no erasures initially
int RS_ENCODER::RSDecodeErasures(GF data[nn], int eras_pos[], int no_eras) {
    RsVerification::verify_received_word(data, nn);

    std::vector<GF> lambda(2 * tt + 1, 0);
    std::vector<GF> b(2 * tt + 1, 0);
    std::vector<GF> omega(2 * tt + 1);
    std::vector<GF> root(2 * tt);
    std::vector<GF> reg(2 * tt + 1);
    std::vector<GF> loc(2 * tt);

    /* first form the syndromes; i.e., evaluate recd(x) at roots of g(x)
     * namely @**(B0+i), i = 0, ... ,(NN-KK-1)
     */
    std::vector<GF> syndromes;
    calculate_syndromes(data, syndromes);

    // Check for errors
    bool syn_error = false;
    for (int i = 1; i <= 2 * tt; i++) {
        if (syndromes[i] != 0) {
            syn_error = true;
        }
    }

    RsVerification::verify_syndromes(syndromes, tt);
    RsVerification::print_syndromes("Initial", syndromes, tt);

    if (!syn_error) {
        /*
         * if syndrome is zero, data[] is a codeword and there are no
         * errors to correct. So return data[] unmodified
         */
        return 0;
    }

    construct_erasure_locator(lambda, eras_pos, no_eras);
    berlekamp_massey(syndromes, lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * tt)
    {
        return RS_ERROR_LAMBDA_ERROR;
    }

    RsVerification::verify_lambda(lambda, deg_lambda);

    /* Find roots of the error+erasure locator polynomial. By Chien Search */
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**(NN-KK)). in index form. Also find deg(omega).
     */
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(B0-1) and den = lambda_pr(inv(X(l))) all in poly-form
     */
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, data) < 0)
    {
        return -1;
    }
    return count;
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
