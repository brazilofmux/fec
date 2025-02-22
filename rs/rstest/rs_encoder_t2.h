#ifndef RS_ENCODER_T2_H
#define RS_ENCODER_T2_H

#include <vector>
#include "rs_encoder_base.h"

class RS_ENCODER_T2 final : public RS_ENCODER_BASE {
public:
    explicit RS_ENCODER_T2(int b0);
    ~RS_ENCODER_T2() override;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;
    static constexpr int TABLE_SIZE = 256 * 4;  // 4 bytes per feedback value
    GF* ptable;

    // Utility function for bit rotation
    static inline uint32_t rotr32(uint32_t value, int shift) {
        return (value >> shift) | (value << (32 - shift));
    }
};

#endif // RS_ENCODER_T2_H
