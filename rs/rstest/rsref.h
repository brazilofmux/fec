#ifndef RSREF_H
#define RSREF_H

#include "rs_common.h"

// Implementation of Reed-Solomon Error Correcting Codes including
// Errors only and Errors+Erasures decoding.
// This is a reference implementation focused on clarity rather than performance.

class RS_ENCODER_REF {
public:
    // Core interface
    RS_ENCODER_REF(int CorrectableErrors);
    ~RS_ENCODER_REF(void);

    void RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]);
    int RSDecode(GF recd[nn]);
    int RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras);

private:
    // Core algorithms
    static void RSGenField(void);
    void RSGenPoly(void);

    // Generator polynomial parameters
    int b0;             // g(x) has roots @^b0, @^(b0+1), ... ,@^(b0+2*tt-1)
    GF gg[2*MAX_TT+1]; // Generator polynomial coefficients
    int tt;            // Number of errors that can be corrected
    int kk;            // Number of data symbols per codeword
};

#endif // RSREF_H
