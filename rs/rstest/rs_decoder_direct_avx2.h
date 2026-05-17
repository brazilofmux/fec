#ifndef RS_DECODER_DIRECT_AVX2_H
#define RS_DECODER_DIRECT_AVX2_H

#include "rs_decoder_base.h"
#include <vector>

// Direct-evaluation syndrome decoder with AVX2 SIMD inner loop on x86-64.
// Processes 32 syndrome lanes per chunk using 256-bit YMM registers.
// On non-AVX2 builds or non-x86 hosts, falls back to a scalar implementation
// using the same chunked table layout so timing comparisons remain fair.
class RS_DECODER_DIRECT_AVX2 : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_DIRECT_AVX2(int tt, int b0);
    ~RS_DECODER_DIRECT_AVX2() = default;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;

    // Whether this instance may use AVX2 at runtime (requires compile-time support + CPU detection).
    bool avx2_usable_;

    // Number of 32-lane chunks needed to cover 2*tt syndromes.
    int n_chunks_;

    // power_table_chunked_[c * (nn * 32) + j * 32 + lane] = alpha^(j * (b0 + c*32 + lane))
    // when c*32 + lane < 2*tt; otherwise 0.
    // Stored in polynomial form for the split-nibble multiply step.
    std::vector<GF> power_table_chunked_;

    // Split-nibble multiply tables come from RS_TABLES (shared across instances).
    // 16-byte rows; the SIMD path loads one row with _mm_loadu_si128 and broadcasts
    // it to 256 bits.
    const GF* split_lo_;
    const GF* split_hi_;
};

#endif // RS_DECODER_DIRECT_AVX2_H
