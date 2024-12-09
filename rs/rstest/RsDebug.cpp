#include <stdio.h>
#include "rs.h"
#include "rs_debug.h"

void RsDebug::init(bool enabled) {
    debug_enabled = enabled;
}
    
void RsDebug::print_syndromes(const char* prefix, GF* syndromes, int tt) {
    if (!debug_enabled) return;
    
    printf("\n%s Syndromes:\n", prefix);
    for (int i = 1; i <= 2*tt; i++) {
        printf("s[%d] = %02X (power form: %d)\n", 
               i, 
               Pow2Poly[syndromes[i]], 
               syndromes[i]);
    }
    printf("\n");
}

void RsDebug::print_berlekamp_step(int step, GF d, int L, GF* lambda, int deg_lambda) {
    if (!debug_enabled) return;
    
    printf("\nBerlekamp Step %d:\n", step);
    printf("  Discrepancy d = %02X (power: %d)\n", Pow2Poly[d], d);
    printf("  L = %d\n", L);
    printf("  deg_lambda = %d\n", deg_lambda);
    printf("  Lambda coefficients (power form):");
    for (int i = 0; i <= deg_lambda; i++) {
        printf(" %d", lambda[i]);
    }
    printf("\n");
}

void RsDebug::print_error_location(int loc, GF error_value, GF original, GF corrected) {
    if (!debug_enabled) return;
    
    printf("Error at position %d:\n", loc);
    printf("  Error value: %02X\n", error_value);
    printf("  Original value: %02X\n", original);
    printf("  Corrected to: %02X\n", corrected);
}

bool RsDebug::debug_enabled = false;
