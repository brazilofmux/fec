#ifndef RS_FACTORY_H
#define RS_FACTORY_H

#include <memory>
#include <map>
#include <functional>
#include "rs_codec.h"
#include "rs_encoder_t1.h"
#include "rs_encoder_t2.h"
#include "rs_encoder_t4.h"
#include "rs_encoder_t8.h"
#include "rs_flipped_encoder_t14.h"
#include "rs_flipped_encoder_t16.h"
#include "rs_encoder_t32.h"
#include "rs_encoder_t64.h"
#include "rs_encoder_general.h"
#include "rs_encoder_general_descending.h"
#include "rs_decoder_t1.h"
#include "rs_decoder_general.h"

//             ---Errors Only---  -Errors+Erasures-
//             ----Decoding-----  ----Decoding-----
// tt  Encode  w/errors  without  w/errors  without
//  1    0.51      0.69      0.70    -----    -----
//  2    0.54      1.10      1.11    -----    -----
//  3    1.29      1.34      1.46    -----    -----
//  4    1.30      1.56      1.56    -----    -----
//  5    1.28      1.82      1.83    -----    -----
//  6    1.28      2.22      2.24    -----    -----
//  7    1.28      2.61      2.39    -----    -----
//  8    1.22      2.81      2.80    -----    -----
// 16    1.19      5.08      5.24    -----    -----
// 24    1.26      7.61      7.69    -----    -----
// 32    1.05     10.34     10.32    -----    -----
// 40    1.08     13.03     13.13    -----    -----
// 48    1.00     17.32     17.30    -----    -----
// 56    0.93     19.79     19.81    -----    -----
// 64    0.83     23.10     23.10    -----    -----

class RS_FACTORY {
public:
    static RS_FACTORY& instance();

    std::unique_ptr<RS_CODEC> create_codec(int tt, int b0);

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
