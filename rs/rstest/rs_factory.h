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
//  1    0.51      2.60      2.60    -----    -----
//  2    0.54      3.72      3.73    -----    -----
//  3    1.29      4.75      4.74    -----    -----
//  4    1.30      5.81      5.81    -----    -----
//  5    1.28      6.91      6.90    -----    -----
//  6    1.28      8.02      7.99    -----    -----
//  7    1.28      9.10      9.08    -----    -----
//  8    1.22     10.20     10.19    -----    -----
// 16    1.14     20.38     20.43    -----    -----
// 24    1.27     29.34     29.28    -----    -----
// 32    1.05     38.58     38.65    -----    -----
// 40    1.10     47.91     47.97    -----    -----
// 48    1.02     57.21     57.08    -----    -----
// 56    0.95     66.30     66.53    -----    -----
// 64    0.83     75.55     75.52    -----    -----

class RS_FACTORY {
public:
    static RS_FACTORY& instance();

    std::unique_ptr<RS_CODEC> create_codec(int tt, int b0);
    std::unique_ptr<RS_CODEC> create_specialized_codec(int tt, int b0);

private:
    RS_FACTORY() = default;

    std::shared_ptr<RS_DECODER_BASE> get_shared_decoder(int tt);
    std::unique_ptr<RS_ENCODER_BASE> create_specialized_encoder(int tt, int b0);
    std::shared_ptr<RS_DECODER_BASE> create_specialized_decoder(int tt);
    std::shared_ptr<RS_DECODER_BASE> shared_decoder;
};

#endif // RS_FACTORY_H
