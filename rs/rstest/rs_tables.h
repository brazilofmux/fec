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
};

#endif // RS_TABLES_H
