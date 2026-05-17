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
