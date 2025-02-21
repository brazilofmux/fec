#ifndef RS_ENCODER_GENERAL_H
#define RS_ENCODER_GENERAL_H

#include "rs_encoder_base.h"

class RS_ENCODER_GENERAL : public RS_ENCODER_BASE {
public:
    RS_ENCODER_GENERAL(int tt, int b0) : RS_ENCODER_BASE(tt, b0) {}
    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override {}
};

#endif
