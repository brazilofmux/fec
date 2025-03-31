#ifndef RS_FLIPPED_ENCODER_T16_H
#define RS_FLIPPED_ENCODER_T16_H

#include <cstdint>
#include <cstring>
#include <vector>
#include "rs_encoder_base.h"

/**
 * Specialized Flipped CRC-style Reed-Solomon encoder for tt=16 (corrects 16 symbol errors)
 * 
 * This implementation uses a CRC-style lookup table approach for Reed-Solomon encoding
 * with the following key characteristics:
 * - Processes data bytes in reverse order (last data byte first)
 * - Uses CRC-style lookup tables for optimal performance
 * - Requires a symmetric generator polynomial (controlled by b0 = (kk+1)/2)
 * - Specialized for exactly tt=16 with hardcoded optimizations
 * 
 * The name "flipped" refers to the fact that bytes are processed in reverse
 * order compared to a traditional Reed-Solomon encoder. This approach works
 * correctly only when the generator polynomial is symmetric, which is achieved
 * by setting b0 = (kk+1)/2.
 * 
 * Historical note: For over 25 years, this implementation was misunderstood
 * as a traditional Reed-Solomon encoder. However, it actually processes bytes
 * in reverse order (CRC-style), which only produces correct results because
 * the generator polynomial is symmetric. This renaming makes the actual behavior
 * explicit.
 * 
 * The CRC-style approach offers significant performance advantages:
 * 1. It allows using lookup tables to avoid polynomial division
 * 2. It enables large-block (64-bit and 32-bit) operations for improved throughput
 * 3. It eliminates the need for polynomial multiplication on each byte
 * 
 * Technical details:
 * - A traditional RS encoder processes data in order from first byte to last
 * - This implementation processes data in CRC fashion, from last byte to first
 * - The algorithm depends on generator polynomial symmetry for correctness
 * - The b0 parameter must be set to (kk+1)/2 to ensure polynomial symmetry
 */
class RS_FLIPPED_ENCODER_T16 final : public RS_ENCODER_BASE {
public:
    explicit RS_FLIPPED_ENCODER_T16(int b0);
    ~RS_FLIPPED_ENCODER_T16() override;

    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) override;

private:
    void RSGenPoly();
    void RSGenTable();

    const GF* pow2poly_;
    const GF* poly2pow_;
    std::vector<GF> gg;
    static constexpr int TABLE_SIZE = 256 * 32;  // 32 bytes per feedback value for tt=16
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

    // Helper for 32-bit operations - processes 4 bytes at once
    static inline void process_block4(GF* bb, const GF* TableRow) {
        // Original DO4 macro:
        // #define DO4(x) {
        //     *((UINT32 *)(bb+x)) = *((UINT32 *)(bb+x+1))
        //                         ^ (*((UINT32 *)(TableRow+x))); }
        uint32_t next_block = *reinterpret_cast<uint32_t*>(bb + 1);
        uint32_t table_block = *reinterpret_cast<const uint32_t*>(TableRow);
        *reinterpret_cast<uint32_t*>(bb) = next_block ^ table_block;
    }

    // Helper for single byte operations
    static inline void process_byte1(GF* bb, const GF* TableRow, int pos) {
        // Original DO1 macro:
        // #define DO1(x) { bb[x] = bb[x+1]^TableRow[x]; }
        bb[pos] = bb[pos + 1] ^ TableRow[pos];
    }
};

#endif // RS_FLIPPED_ENCODER_T16_H