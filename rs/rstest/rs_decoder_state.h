#ifndef RS_DECODER_STATE_H
#define RS_DECODER_STATE_H

#include <memory.h>
#include "rs_common.h"  // For GF type and constants

// Structure to hold intermediate state during Reed-Solomon decoding
// Forms are explicitly defined for each field to ensure consistency
// across all implementations
struct rs_decoder_state {
    // Syndrome stage
    // syndromes[]: Power form after calculation
    // syndromes[i] represents S_i, the i-th syndrome
    GF syndromes[2*MAX_TT+1];
    bool syn_error;  // True if any non-zero syndrome found

    // Berlekamp-Massey stage
    // lambda[]: Power form during calculation
    // lambda[] is the error locator polynomial coefficients
    GF lambda[2*MAX_TT+1];
    int deg_lambda;  // Degree of lambda(x)

    // Chien search stage
    // root[]: Power form - locations where lambda(α^(-i)) = 0
    // Each root[i] represents j where α^j is a root
    GF root[MAX_TT];
    // loc[]: Position form - actual error positions in codeword
    // loc[i] = nn - root[i]
    GF loc[MAX_TT];
    int count;  // Number of roots/locations found

    // Erasure handling
    // phi[]: Power form during calculation
    // phi[] is the erasure locator polynomial coefficients
    GF phi[2*MAX_TT+1];
    int deg_phi;  // Degree of phi(x)

    // Initialize all state to 0
    void clear() {
        memset(this, 0, sizeof(rs_decoder_state));
    }
};

#endif // RS_DECODER_STATE_H
