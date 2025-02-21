#ifndef RS_FACTORY_H
#define RS_FACTORY_H

#include <memory>
#include "rs_encoder_base.h"

class RS_FACTORY {
public:
    static RS_FACTORY& instance();

    std::unique_ptr<RS_ENCODER_BASE> create_encoder(int tt, int b0);

private:
    RS_FACTORY() {}

    // Prevent copying
    RS_FACTORY(const RS_FACTORY&) = delete;
    RS_FACTORY& operator=(const RS_FACTORY&) = delete;
};

#endif // RS_FACTORY_H
