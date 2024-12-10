// Verification helper for RS implementation
#include <stdint.h>
#include <cstdio>
#include "rs.h"
#include "rs_debug.h"

void RsVerify::init(VerifyResults* results) {
    if (!table_initialized) {
        const uint32_t polynomial = 0xEDB88320;
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t c = i;
            for (size_t j = 0; j < 8; j++) {
                if (c & 1) {
                    c = polynomial ^ (c >> 1);
                } else {
                    c >>= 1;
                }
            }
            crc32_table[i] = c;
        }
        table_initialized = true;
    }
    current_results = results;
}

uint32_t RsVerify::calculate_crc32(const void* data, size_t length) {
    if (!table_initialized) init(nullptr);

    const uint8_t* buf = static_cast<const uint8_t*>(data);
    uint32_t crc = 0xFFFFFFFF;

    for (size_t i = 0; i < length; i++) {
        crc = crc32_table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    }

    return crc ^ 0xFFFFFFFF;
}

void RsVerify::verify_tables() {
    if (current_results) {
        uint32_t pow2poly_hash = calculate_crc32(Pow2Poly, sizeof(Pow2Poly));
        uint32_t poly2pow_hash = calculate_crc32(Poly2Pow, sizeof(Poly2Pow));
        current_results->pow2poly_hash = pow2poly_hash;
        current_results->poly2pow_hash = poly2pow_hash;
    }
}

void RsVerify::verify_generator(const GF* gg, int tt) {
    if (current_results) {
        uint32_t gg_hash = calculate_crc32(gg, (2*tt + 1) * sizeof(GF));
        current_results->generator_hash = gg_hash;
    }
}

void RsVerify::verify_syndromes(const GF* syndromes, int tt) {
    if (current_results) {
        uint32_t syndrome_hash = calculate_crc32(syndromes, (2*tt + 1) * sizeof(GF));
        current_results->syndrome_hash = syndrome_hash;
    }
}

void RsVerify::verify_lambda(const GF* lambda, int deg_lambda) {
    if (current_results) {
        uint32_t lambda_hash = calculate_crc32(lambda, (deg_lambda + 1) * sizeof(GF));
        current_results->lambda_hash = lambda_hash;
    }
}

void RsVerify::verify_received_word(const GF* recd, int nn_arg) {
    if (current_results) {
        uint32_t recd_hash = calculate_crc32(recd, nn_arg * sizeof(GF));
        current_results->received_word_hash = recd_hash;
    }
}

uint32_t RsVerify::crc32_table[256];
bool RsVerify::table_initialized = false;
RsVerify::VerifyResults * RsVerify::current_results = nullptr;
