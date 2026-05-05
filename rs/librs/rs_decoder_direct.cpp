#include <cstddef>
#include <vector>
#include "rs_decoder_direct.h"

RS_DECODER_DIRECT::RS_DECODER_DIRECT(int tt, int b0)
    : RS_DECODER_BASE(tt, b0), row_stride_(2 * tt) {
    // Closed-form syndrome:
    //   S_i = sum_{j=0}^{nn-1} recd[j] * alpha^(j * (b0 + i - 1))
    // Store the *log* (exponent) of each multiplier rather than its polynomial
    // form, so the inner loop is one fewer table lookup per element.
    log_table_.resize(static_cast<size_t>(nn) * row_stride_);
    for (int j = 0; j < nn; j++) {
        GF* row = &log_table_[j * row_stride_];
        for (int i = 1; i <= 2 * tt; i++) {
            row[i - 1] = static_cast<GF>(mod_nn_full(j * (b0 + i - 1)));
        }
    }
}

void RS_DECODER_DIRECT::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);

    const int width = 2 * tt_;
    std::vector<GF> S(width, 0);

    for (int j = 0; j < nn; j++) {
        const GF s = recd[j];
        if (s == 0) continue;
        const int log_s = poly2pow_[s];
        const GF* row = &log_table_[j * row_stride_];
        for (int i = 0; i < width; i++) {
            S[i] ^= pow2poly_[mod_nn(row[i] + log_s)];
        }
    }

    for (int i = 0; i < width; i++) {
        syndromes[i + 1] = poly2pow_[S[i]];
    }
}

int RS_DECODER_DIRECT::RSDecode(GF recd[nn]) {
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    bool has_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0) {
            has_error = true;
            break;
        }
    }

    if (!has_error) {
        return 0;
    }

    std::vector<GF> lambda(2 * tt_ + 1, 0);
    lambda[0] = 1;
    berlekamp_massey(syndromes, lambda, 0);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * tt_) {
        return RS_ERROR_LAMBDA_ERROR;
    }

    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }

    return count;
}

int RS_DECODER_DIRECT::RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    bool syn_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0 && !syn_error) {
            syn_error = true;
        }
    }

    if (!syn_error) {
        return 0;
    }

    std::vector<GF> lambda(2 * tt_ + 1, 0);
    construct_erasure_locator(lambda, eras_pos, no_eras);

    berlekamp_massey(syndromes, lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * tt_) {
        return RS_ERROR_LAMBDA_ERROR;
    }

    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }

    return count;
}
