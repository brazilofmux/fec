#ifndef RS_DECODER_GENERAL_H
#define RS_DECODER_GENERAL_H

#include "rs_decoder_base.h"
#include <vector>

class RS_DECODER_GENERAL : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_GENERAL(int tt, int b0);
    ~RS_DECODER_GENERAL() = default;

    int RSDecode(GF recd[nn]) override;
    int RSDecodeErasures(GF recd[nn], int eras_pos[2 * MAX_TT], int no_eras) override;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;
    // Precomputed power tables as class members
    std::vector<std::vector<GF>> precomp_powers;
};

#endif // RS_DECODER_GENERAL_H
