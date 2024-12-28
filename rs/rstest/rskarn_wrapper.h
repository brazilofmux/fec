#ifndef RS_KARN_H
#define RS_KARN_H

#include "rs_common.h"  // For common definitions

// Forward declarations of Karn's original functions
extern "C" {
    void init_rs(void);
    void generate_gf(void);
    void gen_poly(void);
    int encode_rs(unsigned char data[], unsigned char bb[]);
    int eras_dec_rs(unsigned char data[], int eras_pos[], int no_eras);
}

class RS_ENCODER_KARN {
public:
    static void Init(void) {
        if (!initialized) {
            init_rs();  // Initialize Karn's implementation
            initialized = true;
        }
    }

    RS_ENCODER_KARN(int CorrectableErrors) : tt(CorrectableErrors) {
        // Verify parameters match Karn's implementation
        static_assert(n == (1 << MM), "Symbol field size mismatch");
        static_assert(nn == ((1 << MM) - 1), "Block size mismatch");
        
        if (!initialized) {
            Init();
        }
    }

    void RSEncode(GF data[MAX_KK], GF bb[2*MAX_TT]) {
        // Clear the parity buffer first
        memset(bb, 0, 2*tt);
        
        // Call Karn's encoder
        // Note: Karn's implementation uses unsigned char, which matches our GF type
        encode_rs(data, bb);
    }

    int RSDecode(GF recd[nn]) {
        return eras_dec_rs(recd, nullptr, 0);
    }

    int RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras) {
        return eras_dec_rs(recd, eras_pos, no_eras);
    }

private:
    static bool initialized;
    int tt;  // Number of correctable errors
};

bool RS_ENCODER_KARN::initialized = false;

#endif // RS_KARN_H
