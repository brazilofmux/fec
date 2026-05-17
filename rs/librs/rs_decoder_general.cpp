#include <vector>
#include "rs_decoder_general.h"
#include <cstdlib>
#include <stdio.h>

RS_DECODER_GENERAL::RS_DECODER_GENERAL(int tt, int b0) : RS_DECODER_BASE(tt, b0) {
    // Precompute powers for each syndrome
    precomp_powers.resize(2 * tt_ + 1);
    for (int i = 1; i <= 2 * tt_; i++) {
        precomp_powers[i].resize(256);
        for (int j = 0; j < 256; j++) {
            precomp_powers[i][j] = (j == 0) ? 0 : pow2poly_[mod_nn(poly2pow_[j] + b0_ + i - 1)];
        }
    }
}

void RS_DECODER_GENERAL::calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) {
    syndromes.assign(2 * tt_ + 1, 0);

    // Initialize LFSR states
    std::vector<GF> lfsr(2 * tt_ + 1, 0);

    // Process each byte through all LFSRs
    for (int j = nn - 1; j >= 0; j--) {
        for (int i = 1; i <= 2 * tt_; i++) {
            lfsr[i] = precomp_powers[i][lfsr[i]] ^ recd[j];
        }
    }

    // Convert to index form
    for (int i = 1; i <= 2 * tt_; i++) {
        syndromes[i] = poly2pow_[lfsr[i]];
    }
}
