#ifndef RS_DECODER_T1_H
#define RS_DECODER_T1_H

#include <vector>
#include "rs_decoder_base.h"

class RS_DECODER_T1 : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_T1(int b0);
    ~RS_DECODER_T1() = default;

    int RSDecode(GF recd[nn]) override;
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) override;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;
    GF precomp_power1[256];
    GF precomp_power2[256];
};

#endif // RS_DECODER_T1_H
