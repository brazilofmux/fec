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

    int RSDecode(GF recd[nn]) override;
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) override;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;

    // power_table_chunked_[c*(nn*16) + j*16 + lane] = alpha^(j*(b0 + c*16 + lane))
    // when c*16+lane < 2*tt; otherwise 0 (multiplying by zero contributes nothing).
    // Layout: for fixed chunk c, j varies over contiguous 16-byte slabs.
    std::vector<GF> power_table_chunked_;
    int n_chunks_;  // = ceil(2*tt / 16)

    // Split-nibble GF(256) multiply tables (poly 0x11D), 256 scalars * 32 B = 8 KB.
    //   split_lo_[s*16 + k] = s * k         (poly form)
    //   split_hi_[s*16 + k] = s * (k << 4)  (poly form)
    // TODO: these depend only on the GF poly, not on (tt, b0). Hoist into
    // RS_TABLES so all NEON decoder instances share one 8 KB copy.
    std::vector<GF> split_lo_;
    std::vector<GF> split_hi_;
};

#endif // RS_DECODER_DIRECT_NEON_H
