// rs_debug.h
#ifndef RS_DEBUG_H
#define RS_DEBUG_H

#include <stdint.h>
#include "rs.h" // For GF type and other RS definitions

class RsDebug {
public:
    static void init(bool enabled);
    static void print_syndromes(const char* prefix, GF* syndromes, int tt);
    static void print_berlekamp_step(int step, GF d, int L, GF* lambda, int deg_lambda);
    static void print_error_location(int loc, GF error_value, GF original, GF corrected);

private:
    static bool debug_enabled;
};

class RsVerify {
public:
    static void init();
    static uint32_t calculate_crc32(const void* data, size_t length);
    static void verify_tables();
    static void verify_generator(const GF* gg, int tt);
    static void verify_syndromes(const GF* syndromes, int tt);
    static void verify_lambda(const GF* lambda, int deg_lambda);
    static void verify_received_word(const GF* recd, int nn_arg);
private:
    static uint32_t crc32_table[256];
    static bool table_initialized;
};

#endif // RS_DEBUG_H
