// RsVerification.h
#ifndef RS_VERIFICATION_H
#define RS_VERIFICATION_H

#include <string>
#include <vector>
#include <cstdint>
#include "rs.h"

class RsVerification {
public:
    // Combined debug and verification results
    struct VerificationResults {
        // Debug information
        std::vector<std::string> syndromes;
        std::vector<std::string> berlekamp_steps;
        std::vector<std::string> error_locations;

        // Verification hashes
        uint32_t generator_hash;
        uint32_t received_word_hash;
        uint32_t syndrome_hash;
        uint32_t lambda_hash;

        // Clear all results
        void clear() {
            syndromes.clear();
            berlekamp_steps.clear();
            error_locations.clear();
            generator_hash = 0;
            received_word_hash = 0;
            syndrome_hash = 0;
            lambda_hash = 0;
        }
    };

    static void init(bool debug_enabled = true);
    static void setResults(VerificationResults* results);

    // Debug functions
    static void print_syndromes(const char* prefix, const GF* syndromes, int tt);
    static void print_berlekamp_step(int step, GF d, int L, const GF* lambda, int deg_lambda);
    static void print_error_location(int loc, GF error_value, GF original, GF corrected);

    // Verification functions
    static void verify_generator(const GF* gg, int tt);
    static void verify_syndromes(const GF* syndromes, int tt);
    static void verify_lambda(const GF* lambda, int deg_lambda);
    static void verify_received_word(const GF* recd, int nn_arg);

    // Overloads for std::vector
    static void verify_generator(const std::vector<GF>& gg, int tt) {
        verify_generator(gg.data(), tt);
    }
    static void verify_syndromes(const std::vector<GF>& syndromes, int tt) {
        verify_syndromes(syndromes.data(), tt);
    }
    static void verify_lambda(const std::vector<GF>& lambda, int deg_lambda) {
        verify_lambda(lambda.data(), deg_lambda);
    }
    static void print_syndromes(const char* prefix, const std::vector<GF>& syndromes, int tt) {
        print_syndromes(prefix, syndromes.data(), tt);
    }
    static void print_berlekamp_step(int step, GF d, int L, const std::vector<GF>& lambda, int deg_lambda) {
        print_berlekamp_step(step, d, L, lambda.data(), deg_lambda);
    }

    // CRC32 calculation helper
    static uint32_t calculate_crc32(const void* data, size_t length);

    // Comparison helpers
    static bool compare_results(const VerificationResults& left,
                              const VerificationResults& right,
                              bool verbose = false);

private:
    static bool debug_enabled;
    static VerificationResults* current_results;
    static uint32_t crc32_table[256];
    static bool table_initialized;

    static void init_crc_table();
};

#endif // RS_VERIFICATION_H
