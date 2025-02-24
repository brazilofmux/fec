#ifndef RS_FACTORY_H
#define RS_FACTORY_H

#include <memory>
#include "rs_codec.h"
#include "rs_encoder_t1.h"
#include "rs_encoder_t2.h"
#include "rs_encoder_t4.h"
#include "rs_encoder_t8.h"
#include "rs_encoder_t14.h"
#include "rs_encoder_t16.h"
#include "rs_encoder_t32.h"
#include "rs_encoder_t64.h"
#include "rs_encoder_general.h"
#include "rs_decoder_t1.h"
#include "rs_decoder_general.h"

//             ---Errors Only---  -Errors+Erasures-
//             ----Decoding-----  ----Decoding-----
// tt  Encode  w/errors  without  w/errors  without
//  1    0.52      1.25      1.24    -----    -----
//  2    0.55      1.76      1.75    -----    -----
//  3    1.31      2.30      2.30    -----    -----
//  4    1.32      2.98      3.08    -----    -----
//  5    1.30      3.73      3.73    -----    -----
//  6    1.31      4.21      4.17    -----    -----
//  7    1.30      4.85      4.85    -----    -----
//  8    1.24      5.38      5.40    -----    -----
// 16    1.15     10.37     10.44    -----    -----
// 24    1.29     15.42     15.51    -----    -----
// 32    1.08     20.77     20.87    -----    -----
// 40    1.10     26.45     26.35    -----    -----
// 48    1.04     32.78     32.59    -----    -----
// 56    0.94     38.13     38.02    -----    -----
// 64    0.83     43.46     43.38    -----    -----

class RS_FACTORY {
public:
    static RS_FACTORY& instance();

    std::unique_ptr<RS_CODEC> create_codec(int tt, int b0);
    std::unique_ptr<RS_CODEC> create_specialized_codec(int tt, int b0);

private:
    RS_FACTORY() = default;

    std::shared_ptr<RS_DECODER_BASE> get_shared_decoder(int tt, int b0);
    std::unique_ptr<RS_ENCODER_BASE> create_specialized_encoder(int tt, int b0);
    std::shared_ptr<RS_DECODER_BASE> create_specialized_decoder(int tt, int b0);
    std::shared_ptr<RS_DECODER_BASE> shared_decoder;
};

#endif // RS_FACTORY_H
