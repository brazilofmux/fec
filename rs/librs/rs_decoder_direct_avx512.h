#ifndef RS_DECODER_DIRECT_AVX512_H
#define RS_DECODER_DIRECT_AVX512_H

#include "rs_decoder_base.h"
#include <vector>

// Direct-evaluation syndrome decoder with AVX-512 SIMD inner loop on x86-64.
// Processes 64 syndrome lanes per chunk using 512-bit ZMM registers.
// Requires AVX-512F + AVX-512BW (for _mm512_shuffle_epi8).
// On non-AVX512 builds or non-x86 hosts, falls back to a scalar implementation
// using the same chunked table layout so timing comparisons remain fair.
class RS_DECODER_DIRECT_AVX512 : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_DIRECT_AVX512(int tt, int b0);
    ~RS_DECODER_DIRECT_AVX512() = default;

    int RSDecode(GF recd[nn]) override;
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) override;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;

    // Whether this instance may use AVX-512 at runtime.
    bool avx512_usable_;

    // Number of 64-lane chunks needed to cover 2*tt syndromes.
    int n_chunks_;

    // power_table_chunked_[c * (nn * 64) + j * 64 + lane] = alpha^(j * (b0 + c*64 + lane))
    // when c*64 + lane < 2*tt; otherwise 0.
    std::vector<GF> power_table_chunked_;

    // Split-nibble GF(256) multiply tables, 256 scalars * 64 B = 16 KB each.
    // Layout matches the 64-lane chunk width for convenient broadcast + shuffle.
    std::vector<GF> split_lo_;
    std::vector<GF> split_hi_;
};

#endif // RS_DECODER_DIRECT_AVX512_H
