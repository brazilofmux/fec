#ifndef RS_ENCODER_BASE_H
#define RS_ENCODER_BASE_H

#include "rs_common.h"
#include "rs_tables.h"

class RS_ENCODER_BASE {
public:
    virtual ~RS_ENCODER_BASE() = default;

    virtual void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) = 0;

    int get_tt() const { return tt_; }
    int get_kk() const { return kk_; }

protected:
    RS_ENCODER_BASE(int tt, int b0) : tt_(tt), kk_(nn - 2 * tt), b0_(b0) {}

    const int tt_;  // Number of correctable errors
    const int kk_;  // Number of data symbols
    const int b0_;  // Starting position for generator polynomial roots

    static const GF* get_pow2poly() { return RS_TABLES::instance().get_pow2poly(); }
    static const GF* get_poly2pow() { return RS_TABLES::instance().get_poly2pow(); }

    static GF mod_nn(int x) {
        return RS_TABLES::instance().mod_nn(x);
    }
};

#endif // RS_ENCODER_BASE_H
