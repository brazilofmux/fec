#ifndef RS_ENCODER_GENERAL_H
#define RS_ENCODER_GENERAL_H

#include "rs_encoder_base.h"
#include <vector>

class RS_ENCODER_GENERAL : public RS_ENCODER_BASE {
public:
    RS_ENCODER_GENERAL(int tt, int b0);
    ~RS_ENCODER_GENERAL();

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;     // Generator polynomial coefficients
    GF* ptable;             // Lookup table for encoding
};

#endif
