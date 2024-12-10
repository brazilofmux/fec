// rs_debug.h
#ifndef RS_DEBUG_H
#define RS_DEBUG_H

#include <stdint.h>
#include <cstdio>
#include <string>
#include <vector>
#include "rs.h" // For GF type and other RS definitions
#include "rsref.h"

class RsDebug {
public:
    struct DebugResults {
        std::vector<std::string> syndromes;
        std::vector<std::string> berlekamp_steps;
        std::vector<std::string> error_locations;
    };
    
    static void init(bool enabled, DebugResults* results);
    static void print_syndromes(const char* prefix, GF* syndromes, int tt);
    static void print_berlekamp_step(int step, GF d, int L, GF* lambda, int deg_lambda);
    static void print_error_location(int loc, GF error_value, GF original, GF corrected);

private:
    static bool debug_enabled;
    static DebugResults* current_results;
};

class RsVerify {
public:
    struct VerifyResults {
        uint32_t pow2poly_hash;
        uint32_t poly2pow_hash;
        uint32_t generator_hash;
        uint32_t received_word_hash;
        uint32_t syndrome_hash;
        uint32_t lambda_hash;
    };

    static void init(VerifyResults* results);
    static uint32_t calculate_crc32(const void* data, size_t length);
    static void verify_tables();
    static void verify_generator(const GF* gg, int tt);
    static void verify_syndromes(const GF* syndromes, int tt);
    static void verify_lambda(const GF* lambda, int deg_lambda);
    static void verify_received_word(const GF* recd, int nn_arg);
private:
    static uint32_t crc32_table[256];
    static bool table_initialized;
    static VerifyResults* current_results;
};

#endif // RS_DEBUG_H
