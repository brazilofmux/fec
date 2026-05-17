#ifndef RS_DECODER_GENERAL_H
#define RS_DECODER_GENERAL_H

#include "rs_decoder_base.h"
#include <vector>

class RS_DECODER_GENERAL : public RS_DECODER_BASE {
public:
    explicit RS_DECODER_GENERAL(int tt, int b0);
    ~RS_DECODER_GENERAL() = default;

private:
    void calculate_syndromes(const GF recd[nn], std::vector<GF>& syndromes) override;
    // Precomputed power tables as class members
    std::vector<std::vector<GF>> precomp_powers;
};

#endif // RS_DECODER_GENERAL_H
