#ifndef RS_H
#define RS_H

#include "rs_common.h"

// Implementation of Reed-Solomon Error Correcting Codes including
// Errors only and Errors+Erasures decoding.
// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All right reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.

#define DLLExport   __declspec( dllexport )

// Encoder Options
#define SYMMETRICAL_GENERATOR
#define ASSEMBLY_ENCODE
#define ASSEMBLY_DECODE

#ifdef SYMMETRICAL_GENERATOR
//#define SYMMETRICAL_ENCODE
#define ASCENDING_ENCODE
#endif

#define RS_ERROR_NONE          0
#define RS_ERROR_CHIEN_SEARCH -2 // Roots of lamda (by Chien Search) not equal to degree of lamda.
#define RS_ERROR_LAMDA_ERROR  -3 // Degree of lamda too high to yield necessary number of roots.

class RS_ENCODER
{
public:
    static void Init(void);

    RS_ENCODER(int CorrectableErrors);
    void RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]);
    int  RSDecode(GF recd[nn]);
    int  RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras);
    ~RS_ENCODER(void);

    static bool bInitialized;

private:
    static GF ModTable[4*nn+1];
    static GF *pModTable;

    GF *ptable;
    void RSGenTable(void);
    void RSGenPoly(void);

    // Special case encoders
    void RSEncodeGeneral(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode64(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode32(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode16(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode14(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode8(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode4(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode2(GF data[MAX_KK], GF bb[2*MAX_TT]);
    void RSEncode1(GF data[MAX_KK], GF bb[2*MAX_TT]);

    int b0;             // g(x) has roots @^b0, @^(b0+1), ... ,@^(b0+2*tt-1)
    GF gg[2*MAX_TT+1]; // number of non-zero coefficients is 2*tt+1
    int tt;            // Number of errors that can be corrected
    int kk;            // Number of data symbols per codeword
};

#endif // RS_H
