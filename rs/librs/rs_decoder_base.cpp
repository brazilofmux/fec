#include <algorithm>
#include <cstring>
#include <vector>
#include "rs_decoder_base.h"

void RS_DECODER_BASE::berlekamp_massey(const std::vector<GF>& syndromes,
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

int RS_DECODER_BASE::convert_to_index_and_get_degree(std::vector<GF>& poly) {
    int degree = 0;
    for (size_t i = 0; i < poly.size(); i++) {
        poly[i] = poly2pow_[poly[i]];
        if (poly[i] != GF_INFINITY)
            degree = static_cast<int>(i);
    }
    return degree;
}

int RS_DECODER_BASE::chien_search(const std::vector<GF>& lambda, int deg_lambda,
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

int RS_DECODER_BASE::compute_omega(const std::vector<GF>& syndromes,
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

int RS_DECODER_BASE::forney_correction(const std::vector<GF>& omega, int deg_omega,
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

int RS_DECODER_BASE::construct_erasure_locator(std::vector<GF>& lambda, const int* eras_pos, int no_eras) {
    if (no_eras == 0) {
        lambda.assign(lambda.size(), 0);
        lambda[0] = 1;
        return 0; // No erasures
    }

    // Initialize lambda
    lambda.assign(lambda.size(), 0);
    lambda[0] = 1;

    // First erasure position directly sets lambda[1]
    lambda[1] = pow2poly_[eras_pos[0]];

    // Process remaining erasure positions
    for (int i = 1; i < no_eras; i++) {
        GF u = eras_pos[i];
        for (int j = i + 1; j > 0; j--) {
            GF tmp = poly2pow_[lambda[j - 1]];
            if (tmp != GF_INFINITY) {
                lambda[j] ^= pow2poly_[mod_nn(u + tmp)];
            }
        }
    }

    return no_eras; // Degree of the erasure locator polynomial
}

int RS_DECODER_BASE::run_pipeline(std::vector<GF>& syndromes,
                                  const int* eras_pos, int no_eras,
                                  GF recd[nn], DecodeProfile* profile) {
    using clock = std::chrono::high_resolution_clock;
    auto elapsed_us = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double, std::micro>(b - a).count();
    };

    // Early-out on a clean codeword.
    bool has_error = false;
    for (int i = 1; i <= 2 * tt_; i++) {
        if (syndromes[i] != 0) { has_error = true; break; }
    }
    if (!has_error) return 0;

    // Berlekamp-Massey (optionally seeded from the erasure locator).
    auto t_bm0 = clock::now();
    std::vector<GF> lambda(2 * tt_ + 1, 0);
    if (no_eras > 0 && eras_pos != nullptr) {
        construct_erasure_locator(lambda, eras_pos, no_eras);
    } else {
        lambda[0] = 1;
    }
    berlekamp_massey(syndromes, lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda);
    if (profile) profile->berlekamp_massey_us = elapsed_us(t_bm0, clock::now());
    if (deg_lambda > 2 * tt_) return RS_ERROR_LAMBDA_ERROR;

    // Chien search.
    auto t_chien0 = clock::now();
    std::vector<GF> root(2 * tt_);
    std::vector<GF> loc(2 * tt_);
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);
    if (profile) profile->chien_search_us = elapsed_us(t_chien0, clock::now());
    if (deg_lambda != root_count) return RS_ERROR_CHIEN_SEARCH;

    // Omega + Forney correction.
    auto t_forney0 = clock::now();
    std::vector<GF> omega(2 * tt_ + 1);
    int deg_omega = compute_omega(syndromes, lambda, deg_lambda, omega);
    int forney_rc = forney_correction(omega, deg_omega, lambda, deg_lambda, root, count, loc, recd);
    if (profile) profile->forney_us = elapsed_us(t_forney0, clock::now());
    if (forney_rc < 0) return -1;

    return count;
}

int RS_DECODER_BASE::RSDecode(GF recd[nn]) {
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);
    return run_pipeline(syndromes, nullptr, 0, recd, nullptr);
}

int RS_DECODER_BASE::RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);
    return run_pipeline(syndromes, eras_pos, no_eras, recd, nullptr);
}

RS_DECODER_BASE::DecodeProfile RS_DECODER_BASE::profile_decode(GF recd[nn]) {
    using clock = std::chrono::high_resolution_clock;
    DecodeProfile p;
    auto t_start = clock::now();

    auto t_syn0 = clock::now();
    std::vector<GF> syndromes;
    calculate_syndromes(recd, syndromes);
    p.syndrome_us = std::chrono::duration<double, std::micro>(clock::now() - t_syn0).count();

    int rc = run_pipeline(syndromes, nullptr, 0, recd, &p);

    p.total_us     = std::chrono::duration<double, std::micro>(clock::now() - t_start).count();
    p.errors_found = (rc >= 0) ? rc : -1;
    p.success      = (rc >= 0);
    return p;
}

