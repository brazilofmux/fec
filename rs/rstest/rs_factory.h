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
//  1    0.51      4.22      4.21    ---- - ---- -
//  2    0.63      5.61      5.73    ---- - ---- -
//  3    1.22      9.21      9.14    ---- - ---- -
//  4    1.33      7.94      8.01    ---- - ---- -
//  5    1.61     13.47     13.34    ---- - ---- -
//  6    1.82     16.19     16.26    ---- - ---- -
//  7    2.00     19.29     18.68    ---- - ---- -
//  8    1.51     15.04     14.79    ---- - ---- -
// 16    1.53     28.90     28.83    ---- - ---- -
// 24    4.57     69.72     69.82    ---- - ---- -
// 32    1.54     55.01     54.73    ---- - ---- -
// 40    6.13    121.36    121.69    ---- - ---- -
// 48    6.57    152.24    150.15    ---- - ---- -
// 56    6.87    185.50    185.03    ---- - ---- -
// 64    1.58    107.17    108.80    ---- - ---- -

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
