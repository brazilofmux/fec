#include "rs_factory.h"
#include "rs_tables.h"

RS_FACTORY& RS_FACTORY::instance() {
    static RS_FACTORY instance;
    RS_TABLES::instance().ensure_initialized();
    return instance;
}

std::unique_ptr<RS_CODEC> RS_FACTORY::create_codec(int tt, int b0) {
    auto decoder = get_shared_decoder(tt, b0);
    auto encoder = create_specialized_encoder(tt, b0);
    return std::make_unique<RS_CODEC>(std::move(encoder), decoder);
}

// Create a codec with specialized encoder AND decoder
std::unique_ptr<RS_CODEC> RS_FACTORY::create_specialized_codec(int tt, int b0) {
    RS_TABLES::instance().ensure_initialized();

    // Create specialized encoder and decoder
    auto encoder = create_specialized_encoder(tt, b0);
    auto decoder = create_specialized_decoder(tt, b0);

    return std::make_unique<RS_CODEC>(std::move(encoder), decoder);
}

std::shared_ptr<RS_DECODER_BASE> RS_FACTORY::get_shared_decoder(int tt, int b0) {
    if (!shared_decoder || shared_decoder->get_tt() != tt) {
        shared_decoder = std::shared_ptr<RS_DECODER_BASE>(new RS_DECODER_GENERAL(tt, b0));
    }
    return shared_decoder;
}

std::unique_ptr<RS_ENCODER_BASE> RS_FACTORY::create_specialized_encoder(int tt, int b0) {
    switch (tt) {
    case 1:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T1(b0));
    case 2:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T2(b0));
    case 4:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T4(b0));
    case 8:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T8(b0));
    case 14:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T14(b0));
    case 16:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T16(b0));
    case 32:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T32(b0));
    case 64:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_T64(b0));
    default:
        return std::unique_ptr<RS_ENCODER_BASE>(new RS_ENCODER_GENERAL(tt, b0));
    }
}

std::shared_ptr<RS_DECODER_BASE> RS_FACTORY::create_specialized_decoder(int tt, int b0) {
    switch (tt) {
#if 0
    case 1:
        return std::shared_ptr<RS_DECODER_BASE>(new RS_DECODER_T1(tt, b0));
#endif
    default:
        return std::shared_ptr<RS_DECODER_BASE>(new RS_DECODER_GENERAL(tt, b0));
    }
}
