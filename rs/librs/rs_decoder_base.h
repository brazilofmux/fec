#ifndef RS_DECODER_BASE_H
#define RS_DECODER_BASE_H

#include <chrono>
#include <vector>
#include "rs_common.h"
#include "rs_tables.h"

class RS_DECODER_BASE {
public:
    virtual ~RS_DECODER_BASE() = default;

    // Concrete decode entry points. Subclasses customize only calculate_syndromes;
    // the shared pipeline (BM -> Chien -> Forney) lives in the base class. The
    // functions are virtual for the rare case a subclass needs to replace the
    // whole pipeline, but the default implementation is what every Direct/T1/
    // GENERAL decoder uses today.
    virtual int RSDecode(GF recd[nn]);
    virtual int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras);

    int get_tt() const { return tt_; }
    int get_kk() const { return kk_; }

    // Lightweight profiling helper — returns wall time (microseconds) spent in each
    // major phase for a normal (non-erasure) decode of a codeword that has errors.
    // Useful for identifying the next bottleneck once syndromes are fast.
    //
    // Note: recd is modified in place (matches RSDecode's signature).
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

    // Virtual method for syndrome calculation - the one thing subclasses customize.
    virtual void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) = 0;

    // Shared post-syndrome pipeline: BM -> Chien -> Forney. Used by both
    // RSDecode/RSDecodeErasures and profile_decode. Pass eras_pos=nullptr,
    // no_eras=0 for the errors-only case. If profile != nullptr, the per-phase
    // wall times are recorded into it (caller is responsible for syndrome_us
    // and total_us).
    int run_pipeline(std::vector<GF>& syndromes,
                     const int* eras_pos, int no_eras,
                     GF recd[nn], DecodeProfile* profile);

    // Shared decoding methods. Polynomials are raw arrays sized 2*tt+1
    // (bounded by 2*MAX_TT+1) so the pipeline can keep them on the stack.
    int construct_erasure_locator(GF* lambda, int len, const int* eras_pos, int no_eras);
    void berlekamp_massey(const GF* syndromes, GF* lambda, int no_eras);
    int convert_to_index_and_get_degree(GF* poly, int len);
    int chien_search(const GF* lambda, int deg_lambda,
        GF* root, GF* loc, int& count);
    int compute_omega(const GF* syndromes, const GF* lambda,
        int deg_lambda, GF* omega);
    int forney_correction(const GF* omega, int deg_omega,
        const GF* lambda, int deg_lambda,
        const GF* root, int count,
        const GF* loc, GF data[nn]);

    // Protected data
    const int tt_;
    const int kk_;
    const int b0_;
    const GF* pow2poly_;
    const GF* poly2pow_;

    // Static utilities
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
