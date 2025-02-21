#ifndef RS_TABLES_H
#define RS_TABLES_H

#include "rs_common.h"

// Class to manage GF tables
class RS_TABLES {
public:

    static RS_TABLES& instance();

    const GF* get_pow2poly() const { return pow2poly_; }
    const GF* get_poly2pow() const { return poly2pow_; }

    void ensure_initialized();

    GF mod_nn(int x) const {
        return mod_table_[x + nn]; // Access pre-computed table
    }

private:
    RS_TABLES();
    ~RS_TABLES();

    // Prevent copying
    RS_TABLES(const RS_TABLES&) = delete;
    RS_TABLES& operator=(const RS_TABLES&) = delete;

    void initialize_gf();

    bool is_initialized_;
    GF* pow2poly_;
    GF* poly2pow_;

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
