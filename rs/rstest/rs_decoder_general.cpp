#include <vector>
#include "rs_decoder_general.h"
#include <cstdlib>
#include <stdio.h>

RS_DECODER_GENERAL::RS_DECODER_GENERAL(int tt, int b0) : RS_DECODER_BASE(tt, b0) { }

#if 0
// A more mathematical version that doesn't perform well:
//
void RS_DECODER_GENERAL::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);

    for (int i = 1; i <= 2 * tt_; i++) {
        for (int j = 0; j < nn; j++) {
            if (recd[j] != 0) {
                syndromes[i] ^= pow2poly_[mod_nn_full(j * (b0_ + i - 1) + poly2pow_[recd[j]])];
            }
        }
    }

    for (int i = 1; i <= 2 * tt_; i++) {
        syndromes[i] = poly2pow_[syndromes[i]];
    }
}
#endif

#if 1
// Claude 3.7 Sonnet version:
//
void RS_DECODER_GENERAL::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);

    // Pre-calculate all powers for each received byte position
    std::vector<int> base_powers(nn);
    for (int j = 0; j < nn; j++) {
        base_powers[j] = mod_nn_full(j * b0_);
    }

    // Process non-zero bytes only
    for (int j = 0; j < nn; j++) {
        if (recd[j] != 0) {
            const int log_recd = poly2pow_[recd[j]];
            int power = mod_nn_full(base_powers[j] + log_recd);

            for (int i = 1; i <= 2 * tt_; i++) {
                syndromes[i] ^= pow2poly_[power];
                power = mod_nn_full(power + j);
            }
        }
    }

    // Convert to index form at the end
    for (int i = 1; i <= 2 * tt_; i++) {
        syndromes[i] = poly2pow_[syndromes[i]];
    }
}
#elif 1
// My 25+ year old version:
//
void RS_DECODER_GENERAL::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);

    int iPowInit = 0;
    for (int j = 0; j < nn; j++) {
        GF RECD = recd[j];
        if (RECD != 0) {
            int iPow0 = iPowInit + poly2pow_[RECD];
            iPow0 = mod_nn(iPow0);

            for (int i = 1; i <= 2 * tt_; i++) {
                iPow0 = mod_nn(iPow0);
                syndromes[i] ^= pow2poly_[iPow0];
                iPow0 += j;
            }
        }
        iPowInit += b0_;
        iPowInit = mod_nn(iPowInit);
    }

    for (int i = 1; i <= 2 * tt_; i++) {
        syndromes[i] = poly2pow_[syndromes[i]];
    }
}
#endif

int RS_DECODER_GENERAL::RSDecode(GF recd[nn]) {
    // Step 1: Calculate syndromes
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    // Check for errors
    bool has_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0) {
            has_error = true;
            break;
        }
    }

    if (!has_error) {
        return 0; // No errors detected
    }

    // Step 2: Compute error locator polynomial using Berlekamp-Massey
    std::vector<GF> lambda(2 * tt_ + 1, 0);
    lambda[0] = 1; // Start with identity polynomial
    berlekamp_massey(syndromes, lambda, 0);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * tt_) {
        return RS_ERROR_LAMBDA_ERROR;
    }

    // Step 3: Find roots of the error locator polynomial using Chien Search
    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    // Step 4: Compute evaluator polynomial omega
    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    // Step 5: Apply Forney's Algorithm for error correction
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }

    return count;
}

int RS_DECODER_GENERAL::RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
    // Step 1: Calculate syndromes
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    // Check for errors
    bool syn_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0 && !syn_error) {
            syn_error = true;
        }
    }

    if (!syn_error) {
        return 0;  // No errors to correct
    }

    // Step 2: Initialize lambda as erasure locator polynomial
    std::vector<GF> lambda(2 * tt_ + 1, 0);
    construct_erasure_locator(lambda, eras_pos, no_eras);

    // Step 3: Use Berlekamp-Massey to find error+erasure locator polynomial
    berlekamp_massey(syndromes, lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * tt_) {
        return RS_ERROR_LAMBDA_ERROR;
    }

    // Step 4: Find roots using Chien Search
    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    // Step 5: Compute evaluator polynomial omega
    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    // Step 6: Apply Forney Algorithm
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }

    return count;
}
