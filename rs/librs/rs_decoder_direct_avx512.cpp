#include <algorithm>
#include <cstddef>
#include <vector>
#include "rs_decoder_direct_avx512.h"

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#if defined(__AVX512F__) && defined(__AVX512BW__)
#define RS_HAVE_AVX512_COMPILE 1
#endif
#endif

namespace {
#if defined(RS_HAVE_AVX512_COMPILE)
// Runtime CPU detection for AVX-512 (F + BW for the byte shuffle we rely on).
static bool cpu_has_avx512() {
#if defined(__GNUC__) || defined(__clang__)
    // We require both AVX512F (foundation) and AVX512BW (byte/word ops, including shuffle_epi8).
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw");
#elif defined(_MSC_VER)
    // Conservative: assume the binary was built with the right /arch:AVX512 and the CPU is new enough.
    return true;
#else
    return false;
#endif
}
#endif
} // namespace

RS_DECODER_DIRECT_AVX512::RS_DECODER_DIRECT_AVX512(int tt, int b0)
    : RS_DECODER_BASE(tt, b0), avx512_usable_(false),
      split_lo_(RS_TABLES::instance().get_split_lo()),
      split_hi_(RS_TABLES::instance().get_split_hi()) {

#if defined(RS_HAVE_AVX512_COMPILE)
    avx512_usable_ = cpu_has_avx512();
#else
    avx512_usable_ = false;
#endif

    const int width = 2 * tt;
    // 64 lanes per AVX-512 chunk
    n_chunks_ = (width + 63) / 64;

    // Allocate chunked power table: n_chunks * nn * 64 bytes
    power_table_chunked_.assign(static_cast<size_t>(n_chunks_) * nn * 64, 0);

    for (int c = 0; c < n_chunks_; c++) {
        for (int j = 0; j < nn; j++) {
            for (int lane = 0; lane < 64; lane++) {
                const int i = c * 64 + lane;  // syndrome index (0-based)
                if (i < width) {
                    const int e = mod_nn_full(j * (b0 + i));
                    power_table_chunked_[static_cast<size_t>(c) * nn * 64 + j * 64 + lane]
                        = pow2poly_[e];
                }
            }
        }
    }
}

void RS_DECODER_DIRECT_AVX512::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);
    const int width = 2 * tt_;

#if defined(RS_HAVE_AVX512_COMPILE)
    if (avx512_usable_) {
        // AVX-512 fast path — 64 lanes per chunk.
        alignas(64) GF S[2 * MAX_TT + 64] = {0};

        for (int c = 0; c < n_chunks_; c++) {
            __m512i Sreg = _mm512_setzero_si512();
            const GF* base = &power_table_chunked_[static_cast<size_t>(c) * nn * 64];

            for (int j = 0; j < nn; j++) {
                const GF s = recd[j];

                // Load the 16-entry low/high nibble tables and broadcast them to 512 bits.
                // We only need the first 16 bytes; the rest of the 64-byte row is padding.
                const __m128i lo16 = _mm_loadu_si128((const __m128i*)&split_lo_[s * 16]);
                const __m128i hi16 = _mm_loadu_si128((const __m128i*)&split_hi_[s * 16]);

                // Broadcast 128-bit table to full 512-bit register (AVX512F).
                const __m512i lo = _mm512_broadcast_i128(lo16);
                const __m512i hi = _mm512_broadcast_i128(hi16);

                // Load 64 bytes of the power table row for this input position.
                const __m512i x = _mm512_loadu_si512((const void*)&base[j * 64]);

                // Split into low/high nibbles (0x0F mask + logical shift).
                const __m512i mask_lo = _mm512_set1_epi8(0x0F);
                const __m512i xl = _mm512_and_si512(x, mask_lo);
                const __m512i xh = _mm512_srli_epi16(_mm512_andnot_si512(mask_lo, x), 4);

                // 64 parallel 4-bit table lookups using the AVX-512BW shuffle.
                const __m512i contrib_lo = _mm512_shuffle_epi8(lo, xl);
                const __m512i contrib_hi = _mm512_shuffle_epi8(hi, xh);

                Sreg = _mm512_xor_si512(Sreg, _mm512_xor_si512(contrib_lo, contrib_hi));
            }

            _mm512_storeu_si512((__m512i*)&S[c * 64], Sreg);
        }

        for (int i = 0; i < width; i++) {
            syndromes[i + 1] = poly2pow_[S[i]];
        }
        return;
    }
#endif

    // Portable scalar fallback (identical algorithm).
    std::vector<GF> S(width, 0);
    for (int j = 0; j < nn; j++) {
        const GF s = recd[j];
        if (s == 0) continue;
        const int log_s = poly2pow_[s];
        for (int c = 0; c < n_chunks_; c++) {
            const GF* base = &power_table_chunked_[static_cast<size_t>(c) * nn * 64 + j * 64];
            const int lanes = std::min(64, width - c * 64);
            for (int lane = 0; lane < lanes; lane++) {
                const GF p = base[lane];
                if (p != 0) {
                    int e = poly2pow_[p] + log_s;
                    if (e >= nn) e -= nn;
                    S[c * 64 + lane] ^= pow2poly_[e];
                }
            }
        }
    }
    for (int i = 0; i < width; i++) {
        syndromes[i + 1] = poly2pow_[S[i]];
    }
}

