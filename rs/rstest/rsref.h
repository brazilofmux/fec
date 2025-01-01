#ifndef RSREF_H
#define RSREF_H

#include <vector>
#include "rs_common.h"

// Implementation of Reed-Solomon Error Correcting Codes including
// Errors only and Errors+Erasures decoding.
// This is a reference implementation focused on clarity rather than performance.

class RS_ENCODER_REF {
public:
    RS_ENCODER_REF(int CorrectableErrors);
    ~RS_ENCODER_REF() = default;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]);
    int RSDecode(GF recd[nn]);
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras);

private:
    void RSGenPoly();
    void berlekamp_massey(const std::vector<GF>& syndromes, int tt, std::vector<GF>& lambda, int& lambda_deg);
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes);

    std::vector<GF> generator_poly; // Generator polynomial coefficients in power form
    int tt;  // Number of correctable errors
    int kk;  // Number of data symbols per codeword
    int b0;  // Start position for generator polynomial roots
};

#endif // RSREF_H
