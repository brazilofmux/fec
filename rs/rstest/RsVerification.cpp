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

    // Compare syndromes
    if (left.syndromes.size() != right.syndromes.size()) {
        printf("ERROR: Syndrome count mismatch - left: %zu, right: %zu\n",
               left.syndromes.size(), right.syndromes.size());
        match = false;
    } else {
        for (size_t i = 0; i < left.syndromes.size(); i++) {
            if (left.syndromes[i] != right.syndromes[i]) {
                printf("ERROR: Syndrome %zu mismatch:\n", i);
                printf("  left: %s\n", left.syndromes[i].c_str());
                printf("  right: %s\n", right.syndromes[i].c_str());
                match = false;
            } else if (verbose) {
                printf("Syndrome %zu matches: %s\n", i, left.syndromes[i].c_str());
            }
        }
    }

    // Compare Berlekamp steps
    if (left.berlekamp_steps.size() != right.berlekamp_steps.size()) {
        printf("ERROR: Berlekamp step count mismatch - left: %zu, right: %zu\n",
               left.berlekamp_steps.size(), right.berlekamp_steps.size());
        match = false;
    } else {
        for (size_t i = 0; i < left.berlekamp_steps.size(); i++) {
            if (left.berlekamp_steps[i] != right.berlekamp_steps[i]) {
                printf("ERROR: Berlekamp step %zu mismatch:\n", i);
                printf("  left: %s\n", left.berlekamp_steps[i].c_str());
                printf("  right: %s\n", right.berlekamp_steps[i].c_str());
                match = false;
            } else if (verbose) {
                printf("Berlekamp step %zu matches: %s\n", i, 
                       left.berlekamp_steps[i].c_str());
            }
        }
    }

    // Compare error locations
    if (left.error_locations.size() != right.error_locations.size()) {
        printf("ERROR: Error location count mismatch - left: %zu, right: %zu\n",
               left.error_locations.size(), right.error_locations.size());
        match = false;
    } else {
        for (size_t i = 0; i < left.error_locations.size(); i++) {
            if (left.error_locations[i] != right.error_locations[i]) {
                printf("ERROR: Error location %zu mismatch:\n", i);
                printf("  left: %s\n", left.error_locations[i].c_str());
                printf("  right: %s\n", right.error_locations[i].c_str());
                match = false;
            } else if (verbose) {
                printf("Error location %zu matches: %s\n", i, 
                       left.error_locations[i].c_str());
            }
        }
    }

    // Compare verification hashes
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
            printf("ERROR: %s hash mismatch:\n", hash.name);
            printf("  left: 0x%08X\n", hash.left);
            printf("  right: 0x%08X\n", hash.right);
            match = false;
        } else if (verbose) {
            printf("%s hash matches: 0x%08X\n", hash.name, hash.left);
        }
    }

    return match;
}
