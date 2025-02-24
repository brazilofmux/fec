#ifndef RS_ENCODER_T14_H
#define RS_ENCODER_T14_H

#include <cstdint>
#include <cstring>
#include <vector>
#include "rs_encoder_base.h"

class RS_ENCODER_T14 final : public RS_ENCODER_BASE {
public:
    explicit RS_ENCODER_T14(int b0);
    ~RS_ENCODER_T14() override;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;
    static constexpr int TABLE_SIZE = 256 * 28;  // 28 bytes per feedback value for tt=14
    GF* ptable;

    // Helper for 64-bit operations - processes 8 bytes at once
    static inline void process_block8(GF* bb, const GF* TableRow) {
        // Original DO8 macro:
        // #define DO8(x) {
        //     *((UINT64 *)(bb+x)) = *((UINT64 *)(bb+x+1))
        //                         ^ (*((UINT64 *)(TableRow+x))); }
        uint64_t next_block = *reinterpret_cast<uint64_t*>(bb + 1);
        uint64_t table_block = *reinterpret_cast<const uint64_t*>(TableRow);
        *reinterpret_cast<uint64_t*>(bb) = next_block ^ table_block;
    }

    // Helper for single byte operations
    static inline void process_byte1(GF* bb, const GF* TableRow, int pos) {
        // Original DO1 macro:
        // #define DO1(x) { bb[x] = bb[x+1]^TableRow[x]; }
        bb[pos] = bb[pos + 1] ^ TableRow[pos];
    }
};

#endif // RS_ENCODER_T14_H
