#ifndef RS_H
#define RS_H

#include <vector>
#include "rs_common.h"

// Implementation of Reed-Solomon Error Correcting Codes including
// Errors only and Errors+Erasures decoding.
// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All right reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.

#define DLLExport   __declspec( dllexport )

// Encoder Options:
//  - SYMMETRICAL_GENERATOR can only be used when b0 = (kk+1)/2
#define SYMMETRICAL_GENERATOR
#define ASSEMBLY_ENCODE

#ifdef SYMMETRICAL_GENERATOR
// The ASCENDING_ENCODE encoders depend on a symmetrical generator polynomial.
#define ASCENDING_ENCODE
#endif

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

    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes);
    int construct_erasure_locator(std::vector<GF>& lambda, const int* eras_pos, int no_eras);
    void berlekamp_massey(const std::vector<GF>& syndromes, std::vector<GF>& lambda, int no_eras);
    int convert_to_index_and_get_degree(std::vector<GF>& poly);
    int chien_search(const std::vector<GF>& lambda, int deg_lambda, std::vector<GF>& root, std::vector<GF>& loc, int& count);
    int compute_omega(const std::vector<GF>& syndromes, const std::vector<GF>& lambda, int deg_lambda, std::vector<GF>& omega);
    int forney_correction(const std::vector<GF>& omega, int deg_omega, const std::vector<GF>& lambda, int deg_lambda,
                                       const std::vector<GF>& root, int count, const std::vector<GF>& loc, GF data[nn]);

    GF gg[2*MAX_TT+1]; // number of non-zero coefficients is 2*tt+1
    int tt;  // Number of correctable errors
    int kk;  // Number of data symbols per codeword
    int b0;  // Start position for generator polynomial roots
};

#endif // RS_H
