#ifndef RS_DECODER_GENERAL_H
#define RS_DECODER_GENERAL_H

#include "rs_decoder_base.h"
#include <vector>

class RS_DECODER_GENERAL : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_GENERAL(int tt);

    int RSDecode(GF recd[nn]) override;
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) override;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes);
    int construct_erasure_locator(std::vector<GF>& lambda, const int* eras_pos, int no_eras);
    void berlekamp_massey(const std::vector<GF>& syndromes, std::vector<GF>& lambda, int no_eras);
    int convert_to_index_and_get_degree(std::vector<GF>& poly);
    int chien_search(const std::vector<GF>& lambda, int deg_lambda,
        std::vector<GF>& root, std::vector<GF>& loc, int& count);
    int compute_omega(const std::vector<GF>& syndromes, const std::vector<GF>& lambda,
        int deg_lambda, std::vector<GF>& omega);
    int forney_correction(const std::vector<GF>& omega, int deg_omega,
        const std::vector<GF>& lambda, int deg_lambda,
        const std::vector<GF>& root, int count,
        const std::vector<GF>& loc, GF data[nn]);

    const int b0_;  // Starting position for generator roots

private:
    const GF* pow2poly_;
    const GF* poly2pow_;
};

#endif
