#ifndef RS_FLIPPED_ENCODER_T4_H
#define RS_FLIPPED_ENCODER_T4_H

#include <cstdint>
#include <cstring>
#include <vector>
#include "rs_encoder_base.h"

class RS_FLIPPED_ENCODER_T4 final : public RS_ENCODER_BASE {
public:
    explicit RS_FLIPPED_ENCODER_T4(int b0);
    ~RS_FLIPPED_ENCODER_T4() override;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;
    static constexpr int TABLE_SIZE = 256 * 8;  // 8 bytes per feedback value for tt=4
    GF* ptable;

    // Helper for 32-bit operations - processes 4 bytes at once
    static inline void process_block4(GF* bb, const GF* TableRow) {
        // Original DO4 macro converted to function:
        // #define DO4(x) {
        //     *((UINT32 *)(bb+x)) = *((UINT32 *)(bb+x+1))
        //                         ^ (*((UINT32 *)(TableRow+x))); }
        // memcpy compiles to single unaligned load/store instructions at -O1+
        // while staying alignment- and aliasing-clean.
        uint32_t next_block, table_block;
        memcpy(&next_block, bb + 1, sizeof(next_block));
        memcpy(&table_block, TableRow, sizeof(table_block));
        next_block ^= table_block;
        memcpy(bb, &next_block, sizeof(next_block));
    }

    // Helper for single byte operations - processes one byte
    static inline void process_byte1(GF* bb, const GF* TableRow, int pos) {
        // Original DO1 macro converted to function:
        // #define DO1(x) { bb[x] = bb[x+1]^TableRow[x]; }
        bb[pos] = bb[pos + 1] ^ TableRow[pos];
    }
};

#endif // RS_FLIPPED_ENCODER_T4_H
