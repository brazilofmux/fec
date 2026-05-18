#ifndef RS_FACTORY_H
#define RS_FACTORY_H

#include <memory>
#include <map>
#include <functional>
#include "rs_codec.h"
#include "rs_flipped_encoder_t1.h"
#include "rs_flipped_encoder_t2.h"
#include "rs_flipped_encoder_t4.h"
#include "rs_flipped_encoder_t8.h"
#include "rs_flipped_encoder_t14.h"
#include "rs_flipped_encoder_t16.h"
#include "rs_flipped_encoder_t32.h"
#include "rs_flipped_encoder_t64.h"
#include "rs_flipped_encoder_general.h"
#include "rs_standard_encoder_general.h"
#include "rs_decoder_t1.h"
#include "rs_decoder_general.h"
#include "rs_decoder_direct.h"
#include "rs_decoder_direct_neon.h"
#include "rs_decoder_direct_avx2.h"
#include "rs_decoder_direct_avx512.h"

// Per-platform µs/decode benchmark numbers (scalar / NEON / AVX2 / AVX-512)
// live in rs/BENCHMARKS.md. Re-run with `./rstest -b` and append a new
// dated section there rather than editing this header.

class RS_FACTORY {
public:
    enum class DecoderKind {
        Auto,        // factory's default selection (specialized if available, else GENERAL)
        General,     // force RS_DECODER_GENERAL      (LFSR / Horner form)
        Direct,      // force RS_DECODER_DIRECT       (direct evaluation, scalar)
        DirectNeon,   // force RS_DECODER_DIRECT_NEON   (direct evaluation, NEON SIMD)
        DirectAvx2,   // force RS_DECODER_DIRECT_AVX2   (direct evaluation, AVX2 SIMD)
        DirectAvx512, // force RS_DECODER_DIRECT_AVX512 (direct evaluation, AVX-512 SIMD)
    };

    static RS_FACTORY& instance();

    std::unique_ptr<RS_CODEC> create_codec(int tt, int b0,
                                           DecoderKind decoder_kind = DecoderKind::Auto);

    // Direct decoder construction (used by benchmarks / experiments).
    std::shared_ptr<RS_DECODER_BASE> create_decoder(int tt, int b0,
                                                    DecoderKind kind = DecoderKind::Auto);

private:
    RS_FACTORY();

    using EncoderCreator = std::function<std::unique_ptr<RS_ENCODER_BASE>(int b0)>;
    std::map<int, EncoderCreator> flipped_encoder_registry;
    std::map<int, EncoderCreator> standard_encoder_registry;

    using DecoderCreator = std::function<std::shared_ptr<RS_DECODER_BASE>(int b0)>;
    std::map<int, DecoderCreator> decoder_registry;

    void register_implementations();

    std::unique_ptr<RS_ENCODER_BASE> create_best_encoder(int tt, int b0);
    std::shared_ptr<RS_DECODER_BASE> create_best_decoder(int tt, int b0);
};

#endif // RS_FACTORY_H
