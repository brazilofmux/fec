#ifndef RS_CODEC_H
#define RS_CODEC_H

#include "rs_encoder_base.h"
#include "rs_decoder_base.h"
#include <memory>

class RS_CODEC {
public:
    RS_CODEC(std::unique_ptr<RS_ENCODER_BASE> enc,
        std::shared_ptr<RS_DECODER_BASE> dec)
        : encoder(std::move(enc))
        , decoder(dec)
        , tt_(encoder->get_tt())
        , kk_(encoder->get_kk()) {
    }

    // Core operations
    void RSEncode(GF data[MAX_KK], GF bb[2 * MAX_TT]) {
        encoder->RSEncode(data, bb);
    }

    int RSDecode(GF recd[nn]) {
        return decoder->RSDecode(recd);
    }

    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) {
        return decoder->RSDecodeErasures(recd, eras_pos, no_eras);
    }

    // Getters for parameters
    int get_tt() const { return tt_; }
    int get_kk() const { return kk_; }

    // Access to underlying implementations (for testing/comparison)
    const RS_ENCODER_BASE* get_encoder() const { return encoder.get(); }
    const RS_DECODER_BASE* get_decoder() const { return decoder.get(); }

private:
    std::unique_ptr<RS_ENCODER_BASE> encoder;
    std::shared_ptr<RS_DECODER_BASE> decoder;
    const int tt_;  // Cache these values from the encoder
    const int kk_;
};

#endif // RS_CODEC_H
