#ifndef RS_COMMON_H
#define RS_COMMON_H

// Common definitions for Reed-Solomon implementations
typedef unsigned char GF;  // Galois Field element type

#define n   256         // n  = 2^8     size of the field
#define nn  255         // nn = 2^8-1   length of codeword
#define MAX_TT (127)    // maximum number of correctable errors (nn/2)
#define MIN_TT (1)      // minimum number of correctable errors (1)
#define MAX_KK (253)    // Maximum number of data symbols (nn-2*TT_MIN)
#define MIN_KK (1)      // Minimum number of data symbols (nn-2*TT_MAX)

// Common tables used by both implementations
extern GF Pow2Poly[n];  // Power to polynomial form conversion
extern GF Poly2Pow[n];  // Polynomial to power form conversion

static const GF GF_INFINITY = nn;

// Common initialization function
void RS_Init();

#define RS_ERROR_FORNEY_DENOMINATOR (-1)  // Forney algorithm failure: Denominator Lambda'(X_l) is zero, preventing error magnitude calculation (division by zero). Indicates an uncorrectable error pattern.
#define RS_ERROR_CHIEN_SEARCH       (-2)  // Chien search failure: Number of roots found for Lambda(x) does not match its degree.
#define RS_ERROR_LAMBDA_ERROR       (-3)  // Berlekamp-Massey failure: Degree of computed Lambda(x) exceeds the maximum possible (2*tt), cannot be correct.

#endif // RS_COMMON_H