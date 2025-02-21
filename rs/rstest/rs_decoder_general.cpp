// rs_decoder_general.cpp
#include "rs_decoder_general.h"
#include <algorithm>
#include <memory.h>
#include <vector>

RS_DECODER_GENERAL::RS_DECODER_GENERAL(int tt)
    : RS_DECODER_BASE(tt)
    , b0_((nn - 2 * tt + 1) / 2)
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
{
}

int RS_DECODER_GENERAL::RSDecode(GF recd[nn]) {
    // Step 1: Calculate syndromes
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);

    // Check for errors
    bool has_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0 && !has_error) {
            has_error = true;
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

void RS_DECODER_GENERAL::berlekamp_massey(const std::vector<GF>& syndromes,
    std::vector<GF>& lambda, int no_eras) {
    std::vector<GF> b(2 * tt_ + 1);
    std::vector<GF> t(2 * tt_ + 1);

    // Initialize b vector
    std::transform(lambda.begin(), lambda.begin() + 2 * tt_ + 1, b.begin(),
        [this](GF val) { return poly2pow_[val]; });

    int r = no_eras;
    int el = no_eras;
    const int max_iter = 2 * tt_;

    while (++r <= max_iter) {
        // Calculate discrepancy
        GF discr_r = 0;
        const int r_bound = r;
        for (int i = 0; i < r_bound; i++) {
            const GF lambda_i = lambda[i];
            const GF syndrome_r_i = syndromes[r - i];

            if (lambda_i && syndrome_r_i != GF_INFINITY) {
                const int sum = mod_nn_full(poly2pow_[lambda_i] + syndrome_r_i);
                discr_r ^= pow2poly_[sum];
            }
        }

        discr_r = poly2pow_[discr_r];

        if (discr_r == GF_INFINITY) {
            // Shift b vector right by 1
            std::copy_backward(b.begin(), b.end() - 1, b.end());
            b[0] = GF_INFINITY;
        }
        else {
            // Calculate temporary polynomial t
            t[0] = lambda[0];
            for (int i = 0; i < max_iter; i++) {
                t[i + 1] = (b[i] != GF_INFINITY) ?
                    lambda[i + 1] ^ pow2poly_[mod_nn_full(discr_r + b[i])] :
                    lambda[i + 1];
            }

            if (2 * el <= r + no_eras - 1) {
                el = r + no_eras - el;
                // Update b polynomial
                for (int i = 0; i <= max_iter; i++) {
                    b[i] = (lambda[i] == 0) ?
                        GF_INFINITY :
                        mod_nn_full(poly2pow_[lambda[i]] - discr_r + nn);
                }
            }
            else {
                // Shift b vector right by 1
                std::copy_backward(b.begin(), b.end() - 1, b.end());
                b[0] = GF_INFINITY;
            }

            // Update lambda with t
            std::copy(t.begin(), t.end(), lambda.begin());
        }
    }
}

int RS_DECODER_GENERAL::convert_to_index_and_get_degree(std::vector<GF>& poly) {
    int degree = 0;
    for (size_t i = 0; i < poly.size(); i++) {
        poly[i] = poly2pow_[poly[i]];
        if (poly[i] != GF_INFINITY)
            degree = static_cast<int>(i);
    }
    return degree;
}

int RS_DECODER_GENERAL::chien_search(const std::vector<GF>& lambda, int deg_lambda,
    std::vector<GF>& root, std::vector<GF>& loc, int& count) {
    std::vector<GF> reg(2 * tt_ + 1);
    memcpy(&reg[1], &lambda[1], 2 * tt_ * sizeof(GF));

    count = 0;
    for (int i = 1; i <= nn; i++) {
        GF q = 1;
        for (int j = deg_lambda; j > 0; j--) {
            if (reg[j] != GF_INFINITY) {
                reg[j] = mod_nn(reg[j] + j);
                q ^= pow2poly_[reg[j]];
            }
        }
        if (q == 0) {
            root[count] = i;
            loc[count] = nn - i;
            count++;
        }
    }

    return count;
}

int RS_DECODER_GENERAL::compute_omega(const std::vector<GF>& syndromes,
    const std::vector<GF>& lambda,
    int deg_lambda, std::vector<GF>& omega) {
    int deg_omega = 0;
    for (int i = 0; i < 2 * tt_; i++) {
        GF tmp = 0;
        int j = (deg_lambda < i) ? deg_lambda : i;
        for (; j >= 0; j--) {
            if ((syndromes[i + 1 - j] != GF_INFINITY) && (lambda[j] != GF_INFINITY))
                tmp ^= pow2poly_[mod_nn(syndromes[i + 1 - j] + lambda[j])];
        }
        if (tmp != 0)
            deg_omega = i;
        omega[i] = poly2pow_[tmp];
    }
    omega[2 * tt_] = GF_INFINITY;
    return deg_omega;
}

int RS_DECODER_GENERAL::forney_correction(const std::vector<GF>& omega, int deg_omega,
    const std::vector<GF>& lambda, int deg_lambda,
    const std::vector<GF>& root, int count,
    const std::vector<GF>& loc, GF data[nn]) {
    for (int j = count - 1; j >= 0; j--) {
        GF num1 = 0;
        for (int i = deg_omega; i >= 0; i--) {
            if (omega[i] != GF_INFINITY)
                num1 ^= pow2poly_[mod_nn_full(omega[i] + i * root[j])];
        }
        GF num2 = pow2poly_[mod_nn_full(root[j] * (b0_ - 1) + nn)];
        GF den = 0;

        /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
        for (int i = std::min(deg_lambda, 2 * tt_ - 1) & ~1; i >= 0; i -= 2) {
            if (lambda[i + 1] != GF_INFINITY)
                den ^= pow2poly_[mod_nn_full(lambda[i + 1] + i * root[j])];
        }

        if (den == 0) {
            return -1;
        }

        /* Apply error to data */
        if (num1 != 0) {
            data[loc[j]] ^= pow2poly_[mod_nn(poly2pow_[num1] + poly2pow_[num2] + nn - poly2pow_[den])];
        }
    }
    return 0;
}

int RS_DECODER_GENERAL::construct_erasure_locator(std::vector<GF>& lambda, const int* eras_pos, int no_eras) {
    if (no_eras == 0) {
        lambda.assign(lambda.size(), 0);
        lambda[0] = 1;
        return 0; // No erasures
    }

    lambda.assign(lambda.size(), 0);
    lambda[0] = 1;

    for (int i = 0; i < no_eras; i++) {
        GF u = eras_pos[i];
        for (int j = no_eras; j > 0; j--) {
            if (lambda[j - 1] != GF_INFINITY) {
                lambda[j] ^= pow2poly_[mod_nn(u + poly2pow_[lambda[j - 1]])];
            }
        }
    }

    return no_eras; // Degree of the erasure locator polynomial
}

// Main erasure decode implementation
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
