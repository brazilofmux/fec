#ifndef RS_ENCODER_T2_H
#define RS_ENCODER_T2_H

#include "rs_encoder_base.h"

class RS_ENCODER_T2 : public RS_ENCODER_BASE {
public:
    RS_ENCODER_T2(int b0) : RS_ENCODER_BASE(2, b0) {}
    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override {}
};

#endif
