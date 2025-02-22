#ifndef RS_ENCODER_T1_H
#define RS_ENCODER_T1_H

#include <vector>
#include "rs_encoder_base.h"

class RS_ENCODER_T1 final : public RS_ENCODER_BASE {
public:
    explicit RS_ENCODER_T1(int b0);
    ~RS_ENCODER_T1() override;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;
    static constexpr int TABLE_SIZE = 256 * 2;
    GF* ptable;
};

#endif
