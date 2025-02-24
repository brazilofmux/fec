#ifndef RS_DECODING_TEST_H
#define RS_DECODING_TEST_H

#include <random>
#include "rs_factory.h"
#include "RsTestConfig.h"
#include "RsVerification.h"

class RsDecodingTest {
public:
    explicit RsDecodingTest(const RsTestConfig& config);
    bool run();

private:
    struct DecodingStats {
        int total_codewords;
        int total_errors_injected;
        int errors_corrected;
        int decode_failures;
        int errors_by_count[MAX_TT + 1];  // Distribution of error counts
        int uncorrectable_errors;         // When errors > tt

        DecodingStats() : total_codewords(0), total_errors_injected(0),
                         errors_corrected(0), decode_failures(0),
                         uncorrectable_errors(0) {
            memset(errors_by_count, 0, sizeof(errors_by_count));
        }
    };

    const RsTestConfig& config;
    int tt;
    int kk;
    DecodingStats stats;
    std::mt19937 rng;
    std::uniform_real_distribution<double> uniform_dist;
    std::uniform_int_distribution<int> byte_dist;
    std::unique_ptr<RS_CODEC> codec;

    bool process_codeword(const GF* codeword, const int nn_short, const int kk_short);
    void print_stats() const;
    void initialize();
    bool open_and_process_file(const char* filename);
};

#endif // RS_DECODING_TEST_H
