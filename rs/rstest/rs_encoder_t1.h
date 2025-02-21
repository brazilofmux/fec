#ifndef RS_ENCODER_T1_H
#define RS_ENCODER_T1_H

#include "rs_encoder_base.h"

class RS_ENCODER_T1 : public RS_ENCODER_BASE {
public:
    RS_ENCODER_T1(int b0) : RS_ENCODER_BASE(1, b0) {}
    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override {}
};

#endif
