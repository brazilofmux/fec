#include "rs_factory.h"
#include "rs_tables.h"

RS_FACTORY& RS_FACTORY::instance() {
    static RS_FACTORY instance;
    RS_TABLES::instance().ensure_initialized();
    return instance;
}

RS_FACTORY::RS_FACTORY()
{
    register_implementations();
}

std::unique_ptr<RS_CODEC> RS_FACTORY::create_codec(int tt, int b0) {
    auto encoder = create_best_encoder(tt, b0);
    auto decoder = create_best_decoder(tt, b0);
    return std::make_unique<RS_CODEC>(std::move(encoder), decoder);
}

void RS_FACTORY::register_implementations() {
    encoder_registry[1] = [](int b0) { return std::make_unique<RS_ENCODER_T1>(b0); };
    encoder_registry[2] = [](int b0) { return std::make_unique<RS_ENCODER_T2>(b0); };
    encoder_registry[4] = [](int b0) { return std::make_unique<RS_ENCODER_T4>(b0); };
    encoder_registry[8] = [](int b0) { return std::make_unique<RS_ENCODER_T8>(b0); };
    encoder_registry[14] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T14>(b0); };
    encoder_registry[16] = [](int b0) { return std::make_unique<RS_ENCODER_T16>(b0); };
    encoder_registry[32] = [](int b0) { return std::make_unique<RS_ENCODER_T32>(b0); };
    encoder_registry[64] = [](int b0) { return std::make_unique<RS_ENCODER_T64>(b0); };

    decoder_registry[1] = [](int b0) { return std::make_shared<RS_DECODER_T1>(b0); };
}

std::unique_ptr<RS_ENCODER_BASE> RS_FACTORY::create_best_encoder(int tt, int b0) {
    int kk = nn - 2 * tt;
    if (b0 == (kk + 1) / 2) {
        auto it = encoder_registry.find(tt);
        if (it != encoder_registry.end()) {
            return it->second(b0);
        }
        return std::make_unique<RS_ENCODER_GENERAL>(tt, b0);
    }
    else {
        return std::make_unique<RS_ENCODER_GENERAL_DESCENDING>(tt, b0);
    }
}

std::shared_ptr<RS_DECODER_BASE> RS_FACTORY::create_best_decoder(int tt, int b0) {
    auto it = decoder_registry.find(tt);
    if (it != decoder_registry.end()) {
        return it->second(b0);
    }
    return std::make_unique<RS_DECODER_GENERAL>(tt, b0);
}
