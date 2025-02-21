#ifndef RS_DECODER_T1_H
#define RS_DECODER_T1_H

#include "rs_decoder_base.h"

class RS_DECODER_T1 : public RS_DECODER_BASE {
public:
    RS_DECODER_T1(int tt) : RS_DECODER_BASE(tt) {}
    int RSDecode(GF recd[nn]) override { return -1; }
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) override { return -1; }
};

#endif
