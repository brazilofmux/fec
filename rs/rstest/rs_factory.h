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
//  1    0.52      3.56      3.64    -----    -----
//  2    0.62      5.04      5.04    -----    -----
//  3    1.31      6.55      6.56    -----    -----
//  4    1.32      8.08      8.10    -----    -----
//  5    1.46      9.64      9.64    -----    -----
//  6    1.65     11.30     11.20    -----    -----
//  7    1.46     12.76     12.88    -----    -----
//  8    1.45     14.42     14.39    -----    -----
// 16    1.51     28.49     28.61    -----    -----
// 24    1.61     41.61     41.02    -----    -----
// 32    1.54     54.01     54.67    -----    -----
// 40    1.71     66.64     66.88    -----    -----
// 48    1.72     80.12     80.21    -----    -----
// 56    1.71     93.90     93.12    -----    -----
// 64    1.59    106.62    106.42    -----    -----

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
