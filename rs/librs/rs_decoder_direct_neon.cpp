#include <algorithm>
#include <cstddef>
#include <vector>
#include "rs_decoder_direct_neon.h"

#if defined(__aarch64__)
#include <arm_neon.h>
#endif

RS_DECODER_DIRECT_NEON::RS_DECODER_DIRECT_NEON(int tt, int b0)
    : RS_DECODER_BASE(tt, b0),
      split_lo_(RS_TABLES::instance().get_split_lo()),
      split_hi_(RS_TABLES::instance().get_split_hi()) {
    const int width = 2 * tt;
    n_chunks_ = (width + 15) / 16;
    power_table_chunked_.assign(static_cast<size_t>(n_chunks_) * nn * 16, 0);

    for (int c = 0; c < n_chunks_; c++) {
        for (int j = 0; j < nn; j++) {
            for (int lane = 0; lane < 16; lane++) {
                const int i = c * 16 + lane;  // syndrome index (0-based)
                if (i < width) {
                    const int e = mod_nn_full(j * (b0 + i));
                    power_table_chunked_[static_cast<size_t>(c) * nn * 16 + j * 16 + lane]
                        = pow2poly_[e];
                }
            }
        }
    }
}

void RS_DECODER_DIRECT_NEON::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);
    const int width = 2 * tt_;

#if defined(__aarch64__)
    alignas(16) GF S[2 * MAX_TT + 16] = {0};
    const uint8x16_t mask_lo = vdupq_n_u8(0x0F);

    for (int c = 0; c < n_chunks_; c++) {
        uint8x16_t Sreg = vdupq_n_u8(0);
        const GF* base = &power_table_chunked_[static_cast<size_t>(c) * nn * 16];
        for (int j = 0; j < nn; j++) {
            const GF s = recd[j];
            const uint8x16_t lo = vld1q_u8(&split_lo_[s * 16]);
            const uint8x16_t hi = vld1q_u8(&split_hi_[s * 16]);
            const uint8x16_t x  = vld1q_u8(&base[j * 16]);
            const uint8x16_t xl = vandq_u8(x, mask_lo);
            const uint8x16_t xh = vshrq_n_u8(x, 4);
            Sreg = veorq_u8(Sreg, veorq_u8(vqtbl1q_u8(lo, xl), vqtbl1q_u8(hi, xh)));
        }
        vst1q_u8(&S[c * 16], Sreg);
    }

    for (int i = 0; i < width; i++) {
        syndromes[i + 1] = poly2pow_[S[i]];
    }
#else
    // Portable scalar fallback so the file compiles on non-aarch64 hosts.
    std::vector<GF> S(width, 0);
    for (int j = 0; j < nn; j++) {
        const GF s = recd[j];
        if (s == 0) continue;
        const int log_s = poly2pow_[s];
        for (int c = 0; c < n_chunks_; c++) {
            const GF* base = &power_table_chunked_[static_cast<size_t>(c) * nn * 16 + j * 16];
            const int lanes = std::min(16, width - c * 16);
            for (int lane = 0; lane < lanes; lane++) {
                const GF p = base[lane];
                if (p != 0) {
                    int e = poly2pow_[p] + log_s;
                    if (e >= nn) e -= nn;
                    S[c * 16 + lane] ^= pow2poly_[e];
                }
            }
        }
    }
    for (int i = 0; i < width; i++) {
        syndromes[i + 1] = poly2pow_[S[i]];
    }
#endif
}
