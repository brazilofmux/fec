#include <algorithm>
#include <cstddef>
#include <cstring>
#include <vector>
#include "rs_decoder_direct_avx2.h"

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#if defined(__AVX2__)
#define RS_HAVE_AVX2_COMPILE 1
#endif
#endif

namespace {
#if defined(RS_HAVE_AVX2_COMPILE)
// Runtime CPU detection for AVX2 (and F16C, etc. — we only need AVX2 here).
static bool cpu_has_avx2() {
    // __builtin_cpu_supports is a GCC/Clang builtin; MSVC has its own path.
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_cpu_supports("avx2");
#elif defined(_MSC_VER)
    // On MSVC we assume the binary was compiled with /arch:AVX2 or we use
    // IsProcessorFeaturePresent — but for simplicity we rely on the compile
    // flag implying the CPU is new enough. A proper cpuid probe can be added.
    return true;
#else
    return false;
#endif
}
#endif

// Scalar GF multiply used only during table construction (identical to NEON version).
inline GF gf_mul(GF a, GF b, const GF* pow2poly, const GF* poly2pow) {
    if (a == 0 || b == 0) return 0;
    int sum = poly2pow[a] + poly2pow[b];
    if (sum >= nn) sum -= nn;
    return pow2poly[sum];
}
} // namespace

RS_DECODER_DIRECT_AVX2::RS_DECODER_DIRECT_AVX2(int tt, int b0)
    : RS_DECODER_BASE(tt, b0), avx2_usable_(false) {

#if defined(RS_HAVE_AVX2_COMPILE)
    avx2_usable_ = cpu_has_avx2();
#else
    avx2_usable_ = false;
#endif

    const int width = 2 * tt;
    // 32 lanes per AVX2 chunk
    n_chunks_ = (width + 31) / 32;

    // Allocate chunked power table: n_chunks * nn * 32 bytes
    power_table_chunked_.assign(static_cast<size_t>(n_chunks_) * nn * 32, 0);

    for (int c = 0; c < n_chunks_; c++) {
        for (int j = 0; j < nn; j++) {
            for (int lane = 0; lane < 32; lane++) {
                const int i = c * 32 + lane;  // syndrome index (0-based)
                if (i < width) {
                    const int e = mod_nn_full(j * (b0 + i));
                    power_table_chunked_[static_cast<size_t>(c) * nn * 32 + j * 32 + lane]
                        = pow2poly_[e];
                }
            }
        }
    }

    // Split-nibble tables sized for 32 lanes (we only need 16 distinct multipliers,
    // but 32 makes the AVX2 broadcast convenient and matches the chunk width).
    split_lo_.assign(256 * 32, 0);
    split_hi_.assign(256 * 32, 0);
    for (int s = 0; s < 256; s++) {
        for (int k = 0; k < 16; k++) {   // low nibble multipliers 0..15
            split_lo_[s * 32 + k] = gf_mul(static_cast<GF>(s), static_cast<GF>(k),
                                           pow2poly_, poly2pow_);
            split_hi_[s * 32 + k] = gf_mul(static_cast<GF>(s), static_cast<GF>(k << 4),
                                           pow2poly_, poly2pow_);
        }
        // Duplicate the 16-entry table into the upper half so we can treat the whole
        // 32-byte row uniformly if desired (not strictly required with the broadcast below).
        for (int k = 16; k < 32; k++) {
            split_lo_[s * 32 + k] = split_lo_[s * 32 + (k - 16)];
            split_hi_[s * 32 + k] = split_hi_[s * 32 + (k - 16)];
        }
    }
}

void RS_DECODER_DIRECT_AVX2::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);
    const int width = 2 * tt_;

#if defined(RS_HAVE_AVX2_COMPILE)
    if (avx2_usable_) {
        // AVX2 fast path — 32 lanes per chunk, split-nibble GF multiply via shuffle.
        alignas(32) GF S[2 * MAX_TT + 32] = {0};

        for (int c = 0; c < n_chunks_; c++) {
            __m256i Sreg = _mm256_setzero_si256();
            const GF* base = &power_table_chunked_[static_cast<size_t>(c) * nn * 32];

            for (int j = 0; j < nn; j++) {
                const GF s = recd[j];
                // Broadcast the 16-entry (actually 32-byte padded) split tables for this s
                // into 256-bit registers so we can use a single 32-byte shuffle per table.
                const __m128i lo16 = _mm_loadu_si128((const __m128i*)&split_lo_[s * 32]);
                const __m128i hi16 = _mm_loadu_si128((const __m128i*)&split_hi_[s * 32]);
                const __m256i lo   = _mm256_broadcastsi128_si256(lo16);
                const __m256i hi   = _mm256_broadcastsi128_si256(hi16);

                // Load 32 bytes of the power table row for this input position
                const __m256i x = _mm256_loadu_si256((const __m256i*)&base[j * 32]);

                // Split into low/high nibbles
                const __m256i mask_lo = _mm256_set1_epi8(0x0F);
                const __m256i xl = _mm256_and_si256(x, mask_lo);
                const __m256i xh = _mm256_srli_epi16(_mm256_andnot_si256(mask_lo, x), 4);

                // Table lookups: each byte in xl/xh is an index 0..15 into the 16-entry table
                const __m256i contrib_lo = _mm256_shuffle_epi8(lo, xl);
                const __m256i contrib_hi = _mm256_shuffle_epi8(hi, xh);

                Sreg = _mm256_xor_si256(Sreg, _mm256_xor_si256(contrib_lo, contrib_hi));
            }

            _mm256_storeu_si256((__m256i*)&S[c * 32], Sreg);
        }

        for (int i = 0; i < width; i++) {
            syndromes[i + 1] = poly2pow_[S[i]];
        }
        return;
    }
#endif

    // Portable scalar fallback (identical algorithm, just not vectorized).
    // This keeps behavior identical on non-AVX2 platforms and during bring-up.
    std::vector<GF> S(width, 0);
    for (int j = 0; j < nn; j++) {
        const GF s = recd[j];
        if (s == 0) continue;
        const int log_s = poly2pow_[s];
        for (int c = 0; c < n_chunks_; c++) {
            const GF* base = &power_table_chunked_[static_cast<size_t>(c) * nn * 32 + j * 32];
            const int lanes = std::min(32, width - c * 32);
            for (int lane = 0; lane < lanes; lane++) {
                const GF p = base[lane];
                if (p != 0) {
                    int e = poly2pow_[p] + log_s;
                    if (e >= nn) e -= nn;
                    S[c * 32 + lane] ^= pow2poly_[e];
                }
            }
        }
    }
    for (int i = 0; i < width; i++) {
        syndromes[i + 1] = poly2pow_[S[i]];
    }
}

int RS_DECODER_DIRECT_AVX2::RSDecode(GF recd[nn]) {
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    bool has_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0) { has_error = true; break; }
    }
    if (!has_error) return 0;

    std::vector<GF> lambda(2 * tt_ + 1, 0);
    lambda[0] = 1;
    berlekamp_massey(syndromes, lambda, 0);
    int deg_lambda = convert_to_index_and_get_degree(lambda);
    if (deg_lambda > 2 * tt_) return RS_ERROR_LAMBDA_ERROR;

    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);
    if (deg_lambda != root_count) return RS_ERROR_CHIEN_SEARCH;

    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }
    return count;
}

int RS_DECODER_DIRECT_AVX2::RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    bool syn_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0) { syn_error = true; break; }
    }
    if (!syn_error) return 0;

    std::vector<GF> lambda(2 * tt_ + 1, 0);
    construct_erasure_locator(lambda, eras_pos, no_eras);
    berlekamp_massey(syndromes, lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda);
    if (deg_lambda > 2 * tt_) return RS_ERROR_LAMBDA_ERROR;

    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);
    if (deg_lambda != root_count) return RS_ERROR_CHIEN_SEARCH;

    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }
    return count;
}
