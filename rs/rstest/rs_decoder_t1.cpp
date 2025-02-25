#include <algorithm>
#include <memory>
#include <vector>
#include "rs_decoder_t1.h"

RS_DECODER_T1::RS_DECODER_T1(int b0) : RS_DECODER_BASE(1, b0) {
    // Precompute powers for α^b0 and α^(b0+1) in constructor
    for (int i = 0; i < 256; i++) {
        precomp_power1[i] = (i == 0) ? 0 : pow2poly_[mod_nn(poly2pow_[i] + b0_)];
        precomp_power2[i] = (i == 0) ? 0 : pow2poly_[mod_nn(poly2pow_[i] + b0_ + 1)];
    }
}

void RS_DECODER_T1::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);
    GF lfsr1 = 0;
    GF lfsr2 = 0;

    // Roll the loop back up - cleaner, potentially easier for compiler to optimize
    for (int j = nn - 1; j >= 0; j--) {
        lfsr1 = precomp_power1[lfsr1] ^ recd[j];
        lfsr2 = precomp_power2[lfsr2] ^ recd[j];
    }

    syndromes[1] = poly2pow_[lfsr1];
    syndromes[2] = poly2pow_[lfsr2];
}

int RS_DECODER_T1::RSDecode(GF recd[nn]) {
    // Step 1: Calculate syndromes
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    // Check for errors
    bool has_error = false;
    for (int i = 1; i <= 2 * 1; i++) {
        if (syndromes[i] != 0) {
            has_error = true;
            break;
        }
    }

    if (!has_error) {
        return 0; // No errors detected
    }

    // Step 2: Compute error locator polynomial using Berlekamp-Massey
    std::vector<GF> lambda(2 * 1 + 1, 0);
    lambda[0] = 1; // Start with identity polynomial
    berlekamp_massey(syndromes, lambda, 0);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * 1) {
        return RS_ERROR_LAMBDA_ERROR;
    }

    // Step 3: Find roots of the error locator polynomial using Chien Search
    std::vector<GF> root(2 * 1);
    std::vector<GF> loc(2 * 1);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    // Step 4: Compute evaluator polynomial omega
    std::vector<GF> omega(2 * 1 + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    // Step 5: Apply Forney's Algorithm for error correction
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }

    return count;
}

int RS_DECODER_T1::RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
    // Step 1: Calculate syndromes
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    // Check for errors
    bool syn_error = false;
    for (int i = 1; i <= 2 * 1; i++) {
        if (syndromes[i] != 0 && !syn_error) {
            syn_error = true;
        }
    }

    if (!syn_error) {
        return 0;  // No errors to correct
    }

    // Step 2: Initialize lambda as erasure locator polynomial
    std::vector<GF> lambda(2 * 1 + 1, 0);
    construct_erasure_locator(lambda, eras_pos, no_eras);

    // Step 3: Use Berlekamp-Massey to find error+erasure locator polynomial
    berlekamp_massey(syndromes, lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda);

    if (deg_lambda > 2 * 1) {
        return RS_ERROR_LAMBDA_ERROR;
    }

    // Step 4: Find roots using Chien Search
    std::vector<GF> root(2 * 1);
    std::vector<GF> loc(2 * 1);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);

    if (deg_lambda != root_count) {
        return RS_ERROR_CHIEN_SEARCH;
    }

    // Step 5: Compute evaluator polynomial omega
    std::vector<GF> omega(2 * 1 + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);

    // Step 6: Apply Forney Algorithm
    if (forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd) < 0) {
        return -1;
    }

    return count;
}
