#include <stdio.h>
#include "rs.h"
#include "rs_debug.h"

void RsDebug::init(bool enabled, DebugResults* results) {
    debug_enabled = enabled;
    current_results = results;
}

void RsDebug::print_syndromes(const char* prefix, GF* syndromes, int tt) {
    if (!debug_enabled || current_results == nullptr) return;

    char buf[256];
    snprintf(buf, sizeof(buf), "\n%s Syndromes:\n", prefix);
    current_results->syndromes.push_back(buf);
    for (int i = 1; i <= 2*tt; i++) {
        snprintf(buf, sizeof(buf), "s[%d] = %02X (power form: %d)\n", i, Pow2Poly[syndromes[i]], syndromes[i]);
        current_results->syndromes.push_back(buf);
    }
    current_results->syndromes.push_back("\n");
}

void RsDebug::print_berlekamp_step(int step, GF d, int L, GF* lambda, int deg_lambda) {
    if (!debug_enabled || current_results == nullptr) return;

    std::string step_info;
    char buf[512];

    // Format main step info
    snprintf(buf, sizeof(buf), "Berlekamp Step %d:\n  Discrepancy d = %02X (power: %d)\n  L = %d\n  deg_lambda = %d",
             step, Pow2Poly[d], d, L, deg_lambda);
    step_info = buf;

    // Add lambda coefficients
    step_info += "\n  Lambda coefficients (power form):";
    for (int i = 0; i <= deg_lambda; i++) {
        snprintf(buf, sizeof(buf), " %d", lambda[i]);
        step_info += buf;
    }

    current_results->berlekamp_steps.push_back(step_info);
}

void RsDebug::print_error_location(int loc, GF error_value, GF original, GF corrected) {
    if (!debug_enabled || current_results == nullptr) return;

    char buf[256];
    snprintf(buf, sizeof(buf), "Error at position %d:\n  Error value: %02X\n  Original value: %02X\n  Corrected to: %02X",
             loc, error_value, original, corrected);

    current_results->error_locations.push_back(buf);
}

bool RsDebug::debug_enabled = false;
RsDebug::DebugResults * RsDebug::current_results = nullptr;
