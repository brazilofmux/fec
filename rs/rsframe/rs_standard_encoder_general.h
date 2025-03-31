#ifndef RS_STANDARD_ENCODER_GENERAL_H
#define RS_STANDARD_ENCODER_GENERAL_H

#include "rs_encoder_base.h"
#include <cstdint>
#include <cstring>
#include <vector>

class RS_STANDARD_ENCODER_GENERAL : public RS_ENCODER_BASE {
public:
    RS_STANDARD_ENCODER_GENERAL(int tt, int b0);
    ~RS_STANDARD_ENCODER_GENERAL();

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;     // Generator polynomial coefficients
    GF* ptable;             // Lookup table for encoding

    // Helpers for block operations
    static inline void process_block8(GF* bb, const GF* TableRow) {
        uint64_t next_block = *reinterpret_cast<uint64_t*>(bb - 1);
        uint64_t table_block = *reinterpret_cast<const uint64_t*>(TableRow);
        *reinterpret_cast<uint64_t*>(bb) = next_block ^ table_block;
    }

    static inline void process_block4(GF* bb, const GF* TableRow) {
        uint32_t next_block = *reinterpret_cast<uint32_t*>(bb - 1);
        uint32_t table_block = *reinterpret_cast<const uint32_t*>(TableRow);
        *reinterpret_cast<uint32_t*>(bb) = next_block ^ table_block;
    }

    static inline void process_byte1(GF* bb, const GF* TableRow, int pos) {
        bb[pos] = bb[pos - 1] ^ TableRow[pos];
    }
};

#endif // RS_STANDARD_ENCODER_GENERAL_H
