#include <stdio.h>
#include "RsVerification.h"

bool RsVerification::debug_enabled = false;
RsVerification::VerificationResults* RsVerification::current_results = nullptr;
uint32_t RsVerification::crc32_table[256];
bool RsVerification::table_initialized = false;

void RsVerification::init(bool debug_enabled) {
    RsVerification::debug_enabled = debug_enabled;
    init_crc_table();
}

void RsVerification::setResults(VerificationResults* results) {
    current_results = results;
    if (current_results) {
        current_results->clear();
    }
}

void RsVerification::init_crc_table() {
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
}

uint32_t RsVerification::calculate_crc32(const void* data, size_t length) {
    if (!table_initialized) init_crc_table();

    const uint8_t* buf = static_cast<const uint8_t*>(data);
    uint32_t crc = 0xFFFFFFFF;

    for (size_t i = 0; i < length; i++) {
        crc = crc32_table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    }

    return crc ^ 0xFFFFFFFF;
}

void RsVerification::print_syndromes(const char* prefix, GF* syndromes, int tt) {
    if (!debug_enabled || !current_results) return;

    char buf[256];
    snprintf(buf, sizeof(buf), "\n%s Syndromes:\n", prefix);
    current_results->syndromes.push_back(buf);

    for (int i = 1; i <= 2*tt; i++) {
        snprintf(buf, sizeof(buf), "s[%d] = %02X (power form: %d)\n",
                i, Pow2Poly[syndromes[i]], syndromes[i]);
        current_results->syndromes.push_back(buf);
    }
    current_results->syndromes.push_back("\n");
}

void RsVerification::print_berlekamp_step(int step, GF d, int L, GF* lambda, int deg_lambda) {
    if (!debug_enabled || !current_results) return;

    std::string step_info;
    char buf[512];

    snprintf(buf, sizeof(buf),
             "Berlekamp Step %d:\n  Discrepancy d = %02X (power: %d)\n  L = %d\n  deg_lambda = %d",
             step, Pow2Poly[d], d, L, deg_lambda);
    step_info = buf;

    step_info += "\n  Lambda coefficients (power form):";
    for (int i = 0; i <= deg_lambda; i++) {
        snprintf(buf, sizeof(buf), " %d", lambda[i]);
        step_info += buf;
    }

    current_results->berlekamp_steps.push_back(step_info);
}

void RsVerification::print_error_location(int loc, GF error_value, GF original, GF corrected) {
    if (!debug_enabled || !current_results) return;

    char buf[256];
    snprintf(buf, sizeof(buf),
             "Error at position %d:\n  Error value: %02X\n  Original value: %02X\n  Corrected to: %02X",
             loc, error_value, original, corrected);

    current_results->error_locations.push_back(buf);
}

void RsVerification::verify_generator(const GF* gg, int tt) {
    if (!current_results) return;
    current_results->generator_hash = calculate_crc32(gg, (2*tt + 1) * sizeof(GF));
}

void RsVerification::verify_syndromes(const GF* syndromes, int tt) {
    if (!current_results) return;
    current_results->syndrome_hash = calculate_crc32(syndromes, (2*tt + 1) * sizeof(GF));
}

void RsVerification::verify_lambda(const GF* lambda, int deg_lambda) {
    if (!current_results) return;
    current_results->lambda_hash = calculate_crc32(lambda, (deg_lambda + 1) * sizeof(GF));
}

void RsVerification::verify_received_word(const GF* recd, int nn_arg) {
    if (!current_results) return;
    current_results->received_word_hash = calculate_crc32(recd, nn_arg * sizeof(GF));
}

bool RsVerification::compare_results(const VerificationResults& left,
                              const VerificationResults& right,
                              bool verbose) {
    bool match = true;

    // Helper to show full vector comparison
    auto compare_vectors = [verbose](const char* name,
                            const std::vector<std::string>& left,
                            const std::vector<std::string>& right) -> bool {
        bool vectors_match = true;
        if (left.size() != right.size()) {
            printf("\n%s count mismatch:\n", name);
            printf("  Left size:  %zu\n", left.size());
            printf("  Right size: %zu\n", right.size());
            vectors_match = false;
        }

        // Always show contents if verbose or mismatch
        if (!vectors_match || verbose) {
            printf("\n%s contents:\n", name);
            size_t max_size = std::max(left.size(), right.size());
            for (size_t i = 0; i < max_size; i++) {
                if (i < left.size() && i < right.size()) {
                    if (left[i] != right[i]) {
                        printf("\nEntry %zu differs:\n", i);
                        printf("  Left:  %s\n", left[i].c_str());
                        printf("  Right: %s\n", right[i].c_str());
                        vectors_match = false;
                    } else if (verbose) {
                        printf("Entry %zu matches: %s\n", i, left[i].c_str());
                    }
                } else {
                    printf("\nEntry %zu:\n", i);
                    if (i < left.size()) printf("  Left only:  %s\n", left[i].c_str());
                    if (i < right.size()) printf("  Right only: %s\n", right[i].c_str());
                    vectors_match = false;
                }
            }
        }
        return vectors_match;
    };

    printf("\n=== Detailed State Comparison ===\n");

    // Compare all vectors
    match &= compare_vectors("Syndromes", left.syndromes, right.syndromes);
    match &= compare_vectors("Berlekamp Steps", left.berlekamp_steps, right.berlekamp_steps);
    match &= compare_vectors("Error Locations", left.error_locations, right.error_locations);

    // Compare hashes
    printf("\nHash Comparisons:\n");
    struct {
        const char* name;
        uint32_t left;
        uint32_t right;
    } hashes[] = {
        {"Generator", left.generator_hash, right.generator_hash},
        {"Received Word", left.received_word_hash, right.received_word_hash},
        {"Syndrome", left.syndrome_hash, right.syndrome_hash},
        {"Lambda", left.lambda_hash, right.lambda_hash}
    };

    for (const auto& hash : hashes) {
        if (hash.left != hash.right) {
            printf("  %s hash mismatch:\n", hash.name);
            printf("    Left:  0x%08X\n", hash.left);
            printf("    Right: 0x%08X\n", hash.right);
            match = false;
        } else if (verbose) {
            printf("  %s hash matches: 0x%08X\n", hash.name, hash.left);
        }
    }

    printf("\n=== End State Comparison ===\n");
    return match;
}
