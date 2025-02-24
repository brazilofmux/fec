#ifndef RS_ENCODING_TEST_H
#define RS_ENCODING_TEST_H

#include <random>
#include "rs_factory.h"
#include "RsTestConfig.h"
#include "RsVerification.h"

class RsEncodingTest {
public:
    explicit RsEncodingTest(const RsTestConfig& config);
    bool run();

private:
    const RsTestConfig& config;
    int tt;
    int kk;
    std::unique_ptr<RS_CODEC> codec;
    std::mt19937 rng;

    void initialize();
    bool encode_and_save_codewords(const char* filename);
    bool validate_shortened_params(int nn_short, int kk_short) const;
};

#endif // RS_ENCODING_TEST_H
