#include <algorithm>
#include <cstring>
#include <vector>
#include "rs_decoder_base.h"

// Scratch polynomials are sized for the largest supported code so every
// pipeline buffer fits on the stack (255 bytes each) — no heap traffic in
// the decode path.
static const int POLY_MAX = 2 * MAX_TT + 1;

void RS_DECODER_BASE::berlekamp_massey(const GF* syndromes, GF* lambda, int no_eras) {
    // lambda is in polynomial form with lambda[i] == 0 for i > deg(lambda);
    // b is kept in index (log) form with GF_INFINITY marking zero terms.
    // deg_l / deg_b are upper bounds on the degrees, kept tight so the O(t^2)
    // discrepancy/update loops only touch live coefficients (O(t*errors)).
    GF b[POLY_MAX];
    GF t[POLY_MAX];
    const int max_iter = 2 * tt_;

    int deg_l = no_eras;
    for (int i = 0; i <= deg_l; i++) {
        b[i] = poly2pow_[lambda[i]];
    }
    std::fill(b + deg_l + 1, b + max_iter + 1, GF_INFINITY);
    int deg_b = deg_l;

    int r = no_eras;
    int el = no_eras;

    while (++r <= max_iter) {
        // Discrepancy: sum of lambda[i] * S[r-i]. Terms beyond deg(lambda)
        // are zero; terms with i >= r would index S[0], which doesn't exist.
        const int di_max = (deg_l < r - 1) ? deg_l : (r - 1);
        GF discr_r = 0;
        for (int i = 0; i <= di_max; i++) {
            const GF lambda_i = lambda[i];
            const GF syndrome_r_i = syndromes[r - i];

            if (lambda_i && syndrome_r_i != GF_INFINITY) {
                int sum = poly2pow_[lambda_i] + syndrome_r_i;
                if (sum >= nn) sum -= nn;
                discr_r ^= pow2poly_[sum];
            }
        }

        discr_r = poly2pow_[discr_r];

        if (discr_r == GF_INFINITY) {
            // b *= x
            if (deg_b < max_iter) deg_b++;
            memmove(b + 1, b, deg_b * sizeof(GF));
            b[0] = GF_INFINITY;
        }
        else {
            // t = lambda + discr_r * x * b, degree bounded by
            // max(deg(lambda), deg(b) + 1).
            const int ub = std::min(std::max(deg_l, deg_b + 1), max_iter);
            t[0] = lambda[0];
            for (int i = 0; i < ub; i++) {
                const GF bi = b[i];
                if (bi != GF_INFINITY) {
                    int e = discr_r + bi;
                    if (e >= nn) e -= nn;
                    t[i + 1] = lambda[i + 1] ^ pow2poly_[e];
                }
                else {
                    t[i + 1] = lambda[i + 1];
                }
            }

            if (2 * el <= r + no_eras - 1) {
                el = r + no_eras - el;
                // b = lambda / discr_r (index form)
                for (int i = 0; i <= deg_l; i++) {
                    b[i] = (lambda[i] == 0) ?
                        GF_INFINITY :
                        mod_nn(poly2pow_[lambda[i]] - discr_r + nn);
                }
                if (deg_b > deg_l) {
                    std::fill(b + deg_l + 1, b + deg_b + 1, GF_INFINITY);
                }
                deg_b = deg_l;
            }
            else {
                // b *= x
                if (deg_b < max_iter) deg_b++;
                memmove(b + 1, b, deg_b * sizeof(GF));
                b[0] = GF_INFINITY;
            }

            // lambda = t. Entries above ub are unchanged zeros, so only the
            // live prefix needs copying; then re-tighten the degree bound.
            memcpy(lambda, t, (ub + 1) * sizeof(GF));
            deg_l = ub;
            while (deg_l > 0 && lambda[deg_l] == 0) deg_l--;
        }
    }
}

int RS_DECODER_BASE::convert_to_index_and_get_degree(GF* poly, int len) {
    int degree = 0;
    for (int i = 0; i < len; i++) {
        poly[i] = poly2pow_[poly[i]];
        if (poly[i] != GF_INFINITY)
            degree = i;
    }
    return degree;
}

int RS_DECODER_BASE::chien_search(const GF* lambda, int deg_lambda,
    GF* root, GF* loc, int& count) {
    // Compact the live (non-zero) coefficients once: zero terms stay zero
    // across all evaluation points, so the inner loop can be branch-free.
    int idx[POLY_MAX];
    int regv[POLY_MAX];
    int nterms = 0;
    for (int j = 1; j <= deg_lambda; j++) {
        if (lambda[j] != GF_INFINITY) {
            idx[nterms] = j;
            regv[nterms] = lambda[j];
            nterms++;
        }
    }

    count = 0;
    for (int i = 1; i <= nn; i++) {
        GF q = 1;
        for (int k = 0; k < nterms; k++) {
            int v = regv[k] + idx[k];
            if (v >= nn) v -= nn;
            regv[k] = v;
            q ^= pow2poly_[v];
        }
        if (q == 0) {
            root[count] = static_cast<GF>(i);
            loc[count] = static_cast<GF>(nn - i);
            // A degree-d polynomial has at most d roots — stop once all found.
            if (++count == deg_lambda) break;
        }
    }

    return count;
}

int RS_DECODER_BASE::compute_omega(const GF* syndromes, const GF* lambda,
    int deg_lambda, GF* omega) {
    int deg_omega = 0;
    for (int i = 0; i < 2 * tt_; i++) {
        GF tmp = 0;
        const int j_max = (deg_lambda < i) ? deg_lambda : i;
        for (int j = 0; j <= j_max; j++) {
            const GF s = syndromes[i + 1 - j];
            const GF lj = lambda[j];
            if (s != GF_INFINITY && lj != GF_INFINITY) {
                int e = s + lj;
                if (e >= nn) e -= nn;
                tmp ^= pow2poly_[e];
            }
        }
        if (tmp != 0)
            deg_omega = i;
        omega[i] = poly2pow_[tmp];
    }
    omega[2 * tt_] = GF_INFINITY;
    return deg_omega;
}

int RS_DECODER_BASE::forney_correction(const GF* omega, int deg_omega,
    const GF* lambda, int deg_lambda,
    const GF* root, int count,
    const GF* loc, GF data[nn]) {
    for (int j = count - 1; j >= 0; j--) {
        const int rj = root[j];

        // omega evaluated at the root: exponent i*rj stepped incrementally
        // instead of multiplied and folded each term.
        GF num1 = 0;
        int e = 0;
        for (int i = 0; i <= deg_omega; i++) {
            if (omega[i] != GF_INFINITY)
                num1 ^= pow2poly_[mod_nn(omega[i] + e)];
            e += rj;
            if (e >= nn) e -= nn;
        }
        GF num2 = pow2poly_[mod_nn_full(rj * (b0_ - 1) + nn)];

        /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
        GF den = 0;
        const int i_top = std::min(deg_lambda, 2 * tt_ - 1) & ~1;
        int step2 = rj + rj;
        if (step2 >= nn) step2 -= nn;
        int e2 = 0;
        for (int i = 0; i <= i_top; i += 2) {
            if (lambda[i + 1] != GF_INFINITY)
                den ^= pow2poly_[mod_nn(lambda[i + 1] + e2)];
            e2 += step2;
            if (e2 >= nn) e2 -= nn;
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

int RS_DECODER_BASE::construct_erasure_locator(GF* lambda, int len, const int* eras_pos, int no_eras) {
    memset(lambda, 0, len * sizeof(GF));
    lambda[0] = 1;

    if (no_eras == 0) {
        return 0; // No erasures
    }

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

    const int width = 2 * tt_;

    // Early-out on a clean codeword. Syndromes are in index (pow) form, so a
    // zero syndrome is GF_INFINITY — a value of 0 means the syndrome is 1.
    bool has_error = false;
    for (int i = 1; i <= width; i++) {
        if (syndromes[i] != GF_INFINITY) { has_error = true; break; }
    }
    if (!has_error) return 0;

    // Berlekamp-Massey (optionally seeded from the erasure locator).
    auto t_bm0 = clock::now();
    GF lambda[POLY_MAX];
    if (no_eras > 0 && eras_pos != nullptr) {
        construct_erasure_locator(lambda, width + 1, eras_pos, no_eras);
    } else {
        memset(lambda, 0, (width + 1) * sizeof(GF));
        lambda[0] = 1;
    }
    berlekamp_massey(syndromes.data(), lambda, no_eras);
    int deg_lambda = convert_to_index_and_get_degree(lambda, width + 1);
        if (profile) profile->berlekamp_massey_us = elapsed_us(t_bm0, clock::now());
    if (deg_lambda > width) return RS_ERROR_LAMBDA_ERROR;

    // Chien search.
    auto t_chien0 = clock::now();
    GF root[POLY_MAX];
    GF loc[POLY_MAX];
    int count = 0;
    int root_count = chien_search(lambda, deg_lambda, root, loc, count);
        if (profile) profile->chien_search_us = elapsed_us(t_chien0, clock::now());
    if (deg_lambda != root_count) return RS_ERROR_CHIEN_SEARCH;

    // Omega + Forney correction.
    auto t_forney0 = clock::now();
    GF omega[POLY_MAX];
    int deg_omega = compute_omega(syndromes.data(), lambda, deg_lambda, omega);
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
    // Validate erasure inputs before any work: positions index the GF tables
    // directly in construct_erasure_locator, and a count above 2*tt would
    // walk the locator polynomial past its bounds.
    if (no_eras < 0 || no_eras > 2 * tt_)
        return RS_ERROR_INVALID_ERASURES;
    if (no_eras > 0) {
        if (eras_pos == nullptr)
            return RS_ERROR_INVALID_ERASURES;
        for (int i = 0; i < no_eras; i++) {
            if (static_cast<unsigned int>(eras_pos[i]) >= nn)
                return RS_ERROR_INVALID_ERASURES;
        }
    }

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
