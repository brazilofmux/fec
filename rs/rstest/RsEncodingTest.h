#ifndef RS_ENCODING_TEST_H
#define RS_ENCODING_TEST_H

#include <random>
#include "rs.h"
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
    RS_ENCODER* encoder;
    std::mt19937 rng;
    
    void initialize();
    bool encode_and_save_codewords(const char* filename);
    bool validate_shortened_params(int nn_short, int kk_short) const;
};

#endif // RS_ENCODING_TEST_H
