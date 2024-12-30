#ifndef RS_KARN_H
#define RS_KARN_H

#include "rs_common.h"  // For common definitions
#include "rs_decoder_state.h"

// Forward declarations of Karn's original functions
extern "C" {
#include "rskarn.h"
    void init_rs(int tt);
    void generate_gf(void);
    void gen_poly(int tt);
    int encode_rs(unsigned char data[], unsigned char bb[], int tt);
    int eras_dec_rs(unsigned char data[], int eras_pos[], int no_eras, int tt);
}

class RS_ENCODER_KARN {
public:
    static void Init() {
        if (!initialized) {
            initialized = true;
        }
    }

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
#if 1
        return eras_dec_rs(recd, nullptr, 0, tt);
#else
        rs_decoder_state state;
        state.clear();

        // Run through decode stages
        int result = calculate_syndromes(recd, state);
        if (result != 0) return result;

        result = berlekamp_massey(state);
        if (result != 0) return result;

        result = chien_search(state);
        if (result != 0) return result;

        return forney_algorithm(recd, state);
#endif
    }

    int RSDecodeErasures(GF recd[nn], int eras_pos[2*MAX_TT], int no_eras) {
#if 1
        return eras_dec_rs(recd, eras_pos, no_eras, tt);
#else
        rs_decoder_state state;
        state.clear();

        // Calculate erasure locator polynomial first
        int result = calculate_phi(eras_pos, no_eras, state);
        if (result != 0) return result;

        // Then run through normal decode stages
        result = calculate_syndromes(recd, state);
        if (result != 0) return result;

        result = berlekamp_massey(state);
        if (result != 0) return result;

        result = chien_search(state);
        if (result != 0) return result;

        return forney_algorithm(recd, state);
#endif
    }

private:
    // Individual decode stages using Karn's code
    int calculate_syndromes(GF recd[nn], rs_decoder_state& state);
    int berlekamp_massey(rs_decoder_state& state);
    int chien_search(rs_decoder_state& state);
    int forney_algorithm(GF recd[nn], rs_decoder_state& state);
    int calculate_phi(int eras_pos[2*MAX_TT], int no_eras, rs_decoder_state& state);

    static bool initialized;
    int tt;  // Number of correctable errors
};

#endif // RS_KARN_H
