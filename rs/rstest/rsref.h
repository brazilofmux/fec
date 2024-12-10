#ifndef RSREF_H
#define RSREF_H

// Implementation of Reed-Solomon Error Correcting Codes including
// Errors only and Errors+Erasures decoding.
// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All right reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

#define n   256         // n  = 2^8     size of the field
#define nn  255         // nn = 2^8-1   length of codeword
#define MAX_TT (127)    // maximum number of correctable errors (nn/2).
#define MIN_TT (1)      // minimum number of correctable errors (1).
#define MAX_KK (253)    // Maximum number of data symbols (nn-2*TT_MIN).
#define MIN_KK (1)      // Minimum number of data symbols (nn-2*TT_MAX).
typedef unsigned char GF;
extern GF Pow2Poly[n], Poly2Pow[n];
//static const GF GF_INFINITY = nn;

#define DLLExport   __declspec( dllexport )


// Encoder Options: SYMMETRICAL_ENCODE will use half the memory, but
// will execute slightly more slowly. It only works when
// SYMMETRICAL_GENERATOR is defined (i.e., requires a well-chosen b0).
//
// ASSEMBLY_ENCODE will use whatever means to speed execution (including
// assembly) in order to speed execution if a special case of the code is
// available.
//
// ASCENDING_ENCODE changes the 2tt-1 to 0 encoder into a 0 to 2tt-1 encoder.
// The parity is exactly the same, but some processors may benefit from
// accessing ascending memory addresses rather than descending addresses.
// This option requires SYMMETRICAL_GENERATOR (well chosen b0), but doesn't
// not require the SYMMETRICAL_ENCODE flag.
//

#define SYMMETRICAL_GENERATOR
#define ASSEMBLY_ENCODE
#define ASSEMBLY_DECODE

#ifdef SYMMETRICAL_GENERATOR
//#define SYMMETRICAL_ENCODE
#define ASCENDING_ENCODE
#endif

#define RS_ERROR_CHIEN_SEARCH -2 // Roots of lamda (by Chien Search) not equal to degree of lamda.
#define RS_ERROR_LAMDA_ERROR  -3 // Degree of lamda too high to yield necessary number of roots.

class RS_ENCODER_REF
{
public:
    RS_ENCODER_REF(int CorrectableErrors);
    void RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]);
    int  RSDecode(GF recd[nn]);
    int  RSDecodeErasures
        (
            GF recd[nn],
            int eras_pos[2*MAX_TT],
            int no_eras
        );
    ~RS_ENCODER_REF(void);

private:

    void RSGenPoly(void);
    void RSGenTable(void);

    // Special case encoders.
    //
    void RSEncodeGeneral(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode64(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode32(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode16(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode14(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode8(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode4(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode2(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode1(GF data[MAX_KK], GF bb[2*MAX_TT]);

    GF *ptable;

    // Generator polynomial along with the power of the first root. RS codes only
    // require that the root are sequential not that they start at @^0. This
    // flexibility is useful in that the generator polynomial can be made
    // symmetrical by picking a certain starting power, b0. A symetrical generator
    // polynomial, gg, simplifies hardware implementations, and can reduce the
    // tables-sizes of table-driven software implementations by half. It will
    // probably slow software down a little to make use of this symmetry unless
    // a smaller table is the difference between fitting or not fitting into a
    // limited cache.
    //
    int b0;             // g(x) has roots @^b0, @^(b0+1), ... ,@^(b0+2*tt-1)
    GF gg[2*MAX_TT+1];  // number of non-zero coefficients is 2*tt+1.


    int tt;   // Number of errors that can be corrected.
    int kk;   // Number of data symbols per codeword.
};

void RSRef_Init(void);

#endif // RSREF_H
