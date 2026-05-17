#ifndef RS_COMMON_H
#define RS_COMMON_H

// Common definitions for Reed-Solomon implementations
typedef unsigned char GF;  // Galois Field element type

#define GF_N  256       // GF_N = 2^8   size of the Galois field
#define nn  255         // nn = 2^8-1   length of codeword
#define MAX_TT (127)    // maximum number of correctable errors (nn/2)
#define MIN_TT (1)      // minimum number of correctable errors (1)
#define MAX_KK (253)    // Maximum number of data symbols (nn-2*TT_MIN)
#define MIN_KK (1)      // Minimum number of data symbols (nn-2*TT_MAX)

static const GF GF_INFINITY = nn;

#define RS_ERROR_CHIEN_SEARCH -2 // Roots of lamda (by Chien Search) not equal to degree of lamda.
#define RS_ERROR_LAMBDA_ERROR -3 // Degree of lamda too high to yield necessary number of roots.

#endif // RS_COMMON_H
