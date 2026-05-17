#ifndef RS_DECODER_BASE_H
#define RS_DECODER_BASE_H

#include <chrono>
#include <vector>
#include "rs_common.h"
#include "rs_tables.h"

class RS_DECODER_BASE {
public:
    virtual ~RS_DECODER_BASE() = default;

    virtual int RSDecode(GF recd[nn]);
    virtual int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras);

    int get_tt() const { return tt_; }
    int get_kk() const { return kk_; }

    struct DecodeProfile {
        double total_us = 0.0;
        double syndrome_us = 0.0;
        double berlekamp_massey_us = 0.0;
        double chien_search_us = 0.0;
        double forney_us = 0.0;
        int  errors_found = 0;
        bool success = false;
    };

    virtual DecodeProfile profile_decode(GF recd[nn]);

protected:
    RS_DECODER_BASE(int tt, int b0) : tt_(tt), kk_(nn - 2 * tt), b0_(b0),
        pow2poly_(RS_TABLES::instance().get_pow2poly()),
        poly2pow_(RS_TABLES::instance().get_poly2pow()) {
    }

    virtual void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) = 0;

    int run_pipeline(std::vector<GF>& syndromes,
                     const int* eras_pos, int no_eras,
                     GF recd[nn], DecodeProfile* profile);

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

    const int tt_;
    const int kk_;
    const int b0_;
    const GF* pow2poly_;
    const GF* poly2pow_;

    static inline GF mod_nn(int x) { return RS_TABLES::instance().mod_nn(x); }
    static inline int mod_nn_full(int x) {
        while (x >= nn) {
            x -= nn;
            x = (x >> 8) + (x & nn);
        }
        return x;
    }
};

#endif // RS_DECODER_BASE_H
