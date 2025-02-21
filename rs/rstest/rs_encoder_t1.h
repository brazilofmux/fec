#ifndef RS_ENCODER_T1_H
#define RS_ENCODER_T1_H

#include "rs_encoder_base.h"

class RS_ENCODER_T1 final : public RS_ENCODER_BASE {
public:
    explicit RS_ENCODER_T1(int b0);
    ~RS_ENCODER_T1() override;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenTable();

    static constexpr int TABLE_SIZE = 512;  // 256 * 2 for tt=1
    GF* ptable;
};

#endif
