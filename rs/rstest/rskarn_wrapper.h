#ifndef RS_KARN_H
#define RS_KARN_H

#include "rs_common.h"
#include "rskarn.h"

// Forward declarations of Karn's original functions
void init_rs(int tt);
void generate_gf(void);
void gen_poly(int tt);
int encode_rs(unsigned char data[], unsigned char bb[], int tt);
int eras_dec_rs(unsigned char data[], int eras_pos[], int no_eras, int tt);

class RS_ENCODER_KARN {
public:
    RS_ENCODER_KARN(int CorrectableErrors) : tt(CorrectableErrors) {
        // Verify parameters match Karn's implementation
        static_assert(n == (1 << MM), "Symbol field size mismatch");
        static_assert(nn == ((1 << MM) - 1), "Block size mismatch");
        init_rs(tt);  // Initialize Karn's implementation
    }

    void RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]) {
        // Clear the parity buffer first
        memset(bb, 0, 2*tt);

        // Call Karn's encoder
        encode_rs(data, bb, tt);
    }

    // Main decode entry points - these use the staged functions internally
    int RSDecode(GF recd[nn]) {
        return eras_dec_rs(recd, nullptr, 0, tt);
    }

    int RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras) {
        return eras_dec_rs(recd, eras_pos, no_eras, tt);
    }

private:
    int tt;  // Number of correctable errors
};

#endif // RS_KARN_H
