#ifndef RS_DECODER_BASE_H
#define RS_DECODER_BASE_H

#include "rs_common.h"

class RS_DECODER_BASE {
public:
    virtual ~RS_DECODER_BASE() = default;

    virtual int RSDecode(GF recd[nn]) = 0;
    virtual int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) = 0;

    int get_tt() const { return tt_; }
    int get_kk() const { return kk_; }

protected:
    RS_DECODER_BASE(int tt) : tt_(tt), kk_(nn - 2 * tt) {}

    const int tt_;
    const int kk_;

    static const GF* get_pow2poly() { return RS_TABLES::instance().get_pow2poly(); }
    static const GF* get_poly2pow() { return RS_TABLES::instance().get_poly2pow(); }
};

#endif // RS_DECODER_BASE_H
