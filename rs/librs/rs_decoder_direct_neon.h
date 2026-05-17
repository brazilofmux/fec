#ifndef RS_DECODER_DIRECT_NEON_H
#define RS_DECODER_DIRECT_NEON_H

#include "rs_decoder_base.h"
#include <vector>

// Direct-evaluation syndrome decoder with NEON SIMD inner loop on aarch64.
// On non-aarch64 builds falls back to a scalar implementation that uses the
// same chunked precomputed table (so the timing is comparable).
class RS_DECODER_DIRECT_NEON : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_DIRECT_NEON(int tt, int b0);
    ~RS_DECODER_DIRECT_NEON() = default;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;

    // power_table_chunked_[c*(nn*16) + j*16 + lane] = alpha^(j*(b0 + c*16 + lane))
    // when c*16+lane < 2*tt; otherwise 0 (multiplying by zero contributes nothing).
    // Layout: for fixed chunk c, j varies over contiguous 16-byte slabs.
    std::vector<GF> power_table_chunked_;
    int n_chunks_;  // = ceil(2*tt / 16)

    // Split-nibble multiply tables come from RS_TABLES (shared across instances).
    const GF* split_lo_;
    const GF* split_hi_;
};

#endif // RS_DECODER_DIRECT_NEON_H
