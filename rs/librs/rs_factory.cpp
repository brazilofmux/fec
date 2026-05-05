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

std::unique_ptr<RS_CODEC> RS_FACTORY::create_codec(int tt, int b0, DecoderKind decoder_kind) {
    auto encoder = create_best_encoder(tt, b0);
    auto decoder = create_decoder(tt, b0, decoder_kind);
    return std::make_unique<RS_CODEC>(std::move(encoder), decoder);
}

std::shared_ptr<RS_DECODER_BASE> RS_FACTORY::create_decoder(int tt, int b0, DecoderKind kind) {
    switch (kind) {
        case DecoderKind::General:
            return std::make_shared<RS_DECODER_GENERAL>(tt, b0);
        case DecoderKind::Direct:
            return std::make_shared<RS_DECODER_DIRECT>(tt, b0);
        case DecoderKind::DirectNeon:
            return std::make_shared<RS_DECODER_DIRECT_NEON>(tt, b0);
        case DecoderKind::Auto:
        default:
            return create_best_decoder(tt, b0);
    }
}

void RS_FACTORY::register_implementations() {
    flipped_encoder_registry[1] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T1>(b0); };
    flipped_encoder_registry[2] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T2>(b0); };
    flipped_encoder_registry[4] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T4>(b0); };
    flipped_encoder_registry[8] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T8>(b0); };
    flipped_encoder_registry[14] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T14>(b0); };
    flipped_encoder_registry[16] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T16>(b0); };
    flipped_encoder_registry[32] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T32>(b0); };
    flipped_encoder_registry[64] = [](int b0) { return std::make_unique<RS_FLIPPED_ENCODER_T64>(b0); };

    decoder_registry[1] = [](int b0) { return std::make_shared<RS_DECODER_T1>(b0); };
}

std::unique_ptr<RS_ENCODER_BASE> RS_FACTORY::create_best_encoder(int tt, int b0) {
    int kk = nn - 2 * tt;
    // Symmetrical generator polynomial check: b0 == (kk + 1) / 2
    // When the generator polynomial is symmetrical, we can use the more efficient
    // "flipped" encoder implementation which processes feedback in ascending order
    // rather than descending order. This is an optimization that produces identical
    // results to the standard encoder, but only works with symmetrical polynomials.
    if (b0 == (kk + 1) / 2) {
        auto it = flipped_encoder_registry.find(tt);
        if (it != flipped_encoder_registry.end()) {
            return it->second(b0);
        }
        return std::make_unique<RS_FLIPPED_ENCODER_GENERAL>(tt, b0);
    }
    else {
        // For non-symmetrical generator polynomials, we must use the standard encoder
        auto it = standard_encoder_registry.find(tt);
        if (it != standard_encoder_registry.end()) {
            return it->second(b0);
        }
        return std::make_unique<RS_STANDARD_ENCODER_GENERAL>(tt, b0);
    }
}

std::shared_ptr<RS_DECODER_BASE> RS_FACTORY::create_best_decoder(int tt, int b0) {
    auto it = decoder_registry.find(tt);
    if (it != decoder_registry.end()) {
        return it->second(b0);
    }
    return std::make_unique<RS_DECODER_GENERAL>(tt, b0);
}
