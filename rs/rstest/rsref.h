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
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes);
    int construct_erasure_locator(std::vector<GF>& lambda, const int* eras_pos, int no_eras);
    int convert_to_index_and_get_degree(std::vector<GF>& poly);
    int chien_search(const std::vector<GF> lambda, int deg_lambda, std::vector<GF>& root, std::vector<GF>& loc, int& count);
    void forney_correction(const std::vector<GF>& omega, int deg_omega, const std::vector<GF>& lambda, int deg_lambda,
                                       const std::vector<GF>& root, int count, const std::vector<GF>& loc, GF data[nn]);

    std::vector<GF> generator_poly; // Generator polynomial coefficients in power form
    int tt;  // Number of correctable errors
    int kk;  // Number of data symbols per codeword
    int b0;  // Start position for generator polynomial roots
};

#endif // RSREF_H
