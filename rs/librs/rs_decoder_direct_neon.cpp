#include <algorithm>
#include <cstddef>
#include <cstring>
#include <vector>
#include "rs_decoder_direct_neon.h"

#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace {
// One-time scalar GF multiply (poly 0x11D) for table construction.
inline GF gf_mul(GF a, GF b, const GF* pow2poly, const GF* poly2pow) {
    if (a == 0 || b == 0) return 0;
    int sum = poly2pow[a] + poly2pow[b];
    if (sum >= nn) sum -= nn;
    return pow2poly[sum];
}
}

RS_DECODER_DIRECT_NEON::RS_DECODER_DIRECT_NEON(int tt, int b0)
    : RS_DECODER_BASE(tt, b0) {
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

    split_lo_.assign(256 * 16, 0);
    split_hi_.assign(256 * 16, 0);
    for (int s = 0; s < 256; s++) {
        for (int k = 0; k < 16; k++) {
            split_lo_[s * 16 + k] = gf_mul(static_cast<GF>(s), static_cast<GF>(k),
                                           pow2poly_, poly2pow_);
            split_hi_[s * 16 + k] = gf_mul(static_cast<GF>(s), static_cast<GF>(k << 4),
                                           pow2poly_, poly2pow_);
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

int RS_DECODER_DIRECT_NEON::RSDecode(GF recd[nn]) {
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

int RS_DECODER_DIRECT_NEON::RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
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
