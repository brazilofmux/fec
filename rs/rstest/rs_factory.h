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
#include "rs_encoder_general.h"
#include "rs_decoder_t1.h"
#include "rs_decoder_general.h"

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
