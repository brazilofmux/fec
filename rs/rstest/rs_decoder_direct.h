#ifndef RS_DECODER_DIRECT_H
#define RS_DECODER_DIRECT_H

#include "rs_decoder_base.h"
#include <vector>

// Experimental decoder that computes syndromes by direct evaluation rather
// than via per-syndrome LFSRs. The precomputed power table is laid out so
// that, for each input byte position j, the 2*tt syndrome multipliers are
// contiguous in memory — enabling a single fetch (and, in a follow-up, SIMD
// multiply-accumulate) per input byte.
class RS_DECODER_DIRECT : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_DIRECT(int tt, int b0);
    ~RS_DECODER_DIRECT() = default;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;

    // log_table_[j * row_stride_ + (i - 1)] = (j * (b0 + i - 1)) mod nn
    // i.e. the discrete log (base alpha) of the multiplier applied to recd[j]
    // for syndrome i. Row-major: contiguous along i for fixed j.
    std::vector<GF> log_table_;
    int row_stride_;
};

#endif // RS_DECODER_DIRECT_H
