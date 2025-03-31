#include "rs_flipped_encoder_t16.h"
#include <memory>

RS_FLIPPED_ENCODER_T16::RS_FLIPPED_ENCODER_T16(int b0)
    : RS_ENCODER_BASE(16, b0)  // tt is always 16
    , pow2poly_(get_pow2poly())
    , poly2pow_(get_poly2pow())
    , gg(2 * tt_ + 1)
    , ptable(new GF[TABLE_SIZE])
{
    RSGenPoly();
    RSGenTable();
}

RS_FLIPPED_ENCODER_T16::~RS_FLIPPED_ENCODER_T16() {
    delete[] ptable;
}

void RS_FLIPPED_ENCODER_T16::RSGenPoly() {
    // Generate the generator polynomial g(x) for Reed-Solomon encoding
    // For the flipped encoder to work correctly, this polynomial must be symmetric
    // This is achieved by setting b0 = (kk+1)/2 in the factory
    
    // Initially, g(x) = (X+@^b0)
    gg[0] = pow2poly_[b0_];
    gg[1] = 1;

    // Build the complete generator polynomial as a product of factors
    // g(x) = (x+a^b0)(x+a^(b0+1))...(x+a^(b0+2t-1))
    for (int i = 2; i <= 2 * tt_; i++) {
        gg[i] = 1;

        // Multiply (gg[0]+gg[1]*x + ... +gg[i]x^i) by (@^(b0+i-1) + x)
        for (int j = i - 1; j > 0; j--) {
            if (gg[j] != 0) {
                int iMod = mod_nn(poly2pow_[gg[j]] + b0_ + i - 1);
                gg[j] = gg[j - 1] ^ pow2poly_[iMod];
            }
            else {
                gg[j] = gg[j - 1];
            }
        }

        // gg[0] can never be zero
        int iMod = mod_nn(poly2pow_[gg[0]] + b0_ + i - 1);
        gg[0] = pow2poly_[iMod];
    }

    // Convert gg[] to power form for quicker encoding
    // Using power form simplifies multiplication operations
    for (int i = 0; i <= 2 * tt_; i++) {
        gg[i] = poly2pow_[gg[i]];
    }
    
    // The resulting generator polynomial has a critical property:
    // When b0 = (kk+1)/2, the polynomial is symmetric, meaning:
    // g_i = g_(2t-i) for all i from 0 to 2t
    // This symmetry is what allows the flipped CRC-style algorithm to work correctly
}

void RS_FLIPPED_ENCODER_T16::RSGenTable() {
    // Generate lookup tables for the CRC-style encoding approach
    // This is a critical optimization that eliminates the need for 
    // polynomial multiplication during encoding
    
    // For each possible feedback value (0-255), we pre-calculate
    // the contributions to each parity byte
    for (int ch = 0; ch < 256; ch++) {
        int feedback = poly2pow_[ch];

        // For tt=16, we generate 32 coefficients per feedback value (2*tt=32)
        // Each entry represents how this feedback value affects each parity byte
        for (int j = 0; j < 32; j++) {
            if (gg[j + 1] != GF_INFINITY && feedback != GF_INFINITY) {
                // Calculate the polynomial coefficient for this feedback/position pair
                int iMod = mod_nn(gg[j + 1] + feedback);
                ptable[ch * 32 + j] = pow2poly_[iMod];
            }
            else {
                ptable[ch * 32 + j] = 0;
            }
        }
    }
    
    // The table now contains pre-calculated values for all possible feedback bytes:
    // - Each row represents a possible feedback value (0-255)
    // - Each column represents the effect on a specific parity position
    // - During encoding, we simply XOR these pre-calculated values with the parity buffer
    // - This approach is similar to how CRC checksums are calculated
}

void RS_FLIPPED_ENCODER_T16::RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
    // Clear parity buffer
    memset(bb, 0, 2 * tt_);

    // Process all data bytes in CRC fashion
    // This approach differs from traditional Reed-Solomon encoding:
    // - Traditional RS builds polynomials by multiplying by x and adding terms
    // - This "flipped" implementation processes bytes like a CRC calculation
    // - It works correctly only because we have a symmetric generator polynomial
    // - This method offers significant performance optimizations through lookup tables
    for (int i = 0; i < kk_; i++) {
        // Calculate feedback term (data byte XOR highest parity byte)
        const GF feedback = bb[0] ^ data[i];
        
        // Lookup pre-calculated values for this feedback byte
        GF* TableRow = ptable + 32 * feedback;

        // Process parity buffer in large blocks for performance
        // This is the key optimization: instead of polynomial arithmetic,
        // we use table lookups and XOR operations
        
        // Process three 8-byte blocks using 64-bit operations
        process_block8(bb, TableRow);
        process_block8(bb + 8, TableRow + 8);
        process_block8(bb + 16, TableRow + 16);

        // Process one 4-byte block
        process_block4(bb + 24, TableRow + 24);

        // Process remaining bytes individually
        process_byte1(bb, TableRow, 28);
        process_byte1(bb, TableRow, 29);
        process_byte1(bb, TableRow, 30);

        // Final byte is feedback
        bb[31] = feedback;
    }
}