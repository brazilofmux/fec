#ifndef RS_TABLES_H
#define RS_TABLES_H

#include "rs_common.h"

class RS_TABLES {
public:

    static RS_TABLES& instance();

    const GF* get_pow2poly() const { return pow2poly_; }
    const GF* get_poly2pow() const { return poly2pow_; }

    // Split-nibble GF(256) multiply tables used by the SIMD Direct decoders
    // (NEON, AVX2, AVX-512). Both are 256 scalars × 16 bytes per row:
    //   get_split_lo()[s * 16 + k] = s * k        (poly form)   for k in 0..15
    //   get_split_hi()[s * 16 + k] = s * (k << 4) (poly form)   for k in 0..15
    // Depends only on the GF polynomial, so it's shared across all decoder
    // instances. SIMD code loads one 16-byte row with vld1q_u8 / _mm_loadu_si128
    // and broadcasts it before the shuffle-based multiply.
    const GF* get_split_lo() const { return split_lo_; }
    const GF* get_split_hi() const { return split_hi_; }

    void ensure_initialized();

    GF mod_nn(int x) const {
        return mod_table_[x + nn];
    }

private:
    RS_TABLES();
    ~RS_TABLES();

    // Prevent copying
    RS_TABLES(const RS_TABLES&) = delete;
    RS_TABLES& operator=(const RS_TABLES&) = delete;

    void initialize_gf();

    void initialize_split_tables();

    bool is_initialized_;
    GF* pow2poly_;
    GF* poly2pow_;
    GF* split_lo_;  // 256 * 16 bytes
    GF* split_hi_;  // 256 * 16 bytes

    void initialize_mod_table() {
        mod_table_ = new GF[4 * nn + 1];  // -nn to 3*nn range
        GF* pModTable = mod_table_ + nn;  // Center the table

        // Create Modulus Table
        for (int i = -nn; i <= 3 * nn; i++) {
            pModTable[i] = (nn + i) % nn;
        }
    }

    GF* mod_table_;
};

#endif // RS_TABLES_H
