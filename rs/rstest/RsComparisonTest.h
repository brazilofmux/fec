#ifndef RS_COMPARISON_TEST_H
#define RS_COMPARISON_TEST_H

#include "rs.h"
#include "rsref.h"
#include "rskarn_wrapper.h"
#include "RsTestConfig.h"
#include "RsVerification.h"

enum class Implementation {
    Original,
    Reference,
    Karn
};

class RsComparisonTest {
public:
    explicit RsComparisonTest(const RsTestConfig& config);
    bool run();

private:
    struct TestContext {
        RsVerification::VerificationResults results;
        GF recd[nn];
        int decode_result;
        int eras_pos[2*MAX_TT];
        int no_eras;
    };

    const RsTestConfig& config;
    int tt;
    GF test_data[nn];

    void prepare_test_data();
    TestContext run_encoder(Implementation impl, bool use_erasures) const;
    bool compare_results(const TestContext& left, const TestContext& right,
                        const char* test_name) const;

    // New test methods
    bool test_implementations_match() ;
    struct EncodingTestContext {
        GF data[MAX_KK];
        GF bb[2*MAX_TT];
        RsVerification::VerificationResults results;
    };
    EncodingTestContext run_encoder_test(bool use_reference) const;
    bool compare_encoding_results(const EncodingTestContext& left,
                                const EncodingTestContext& right,
                                const char* test_name) const;

    // Decoder execution helpers
    int run_decoder(Implementation impl, GF recd[nn]);
    int run_erasure_decoder(Implementation impl, TestContext& ctx);
    std::string get_impl_name(Implementation impl);

    // Comparison runners
    bool run_decoder_comparison(Implementation impl_a, Implementation impl_b);
    bool run_erasure_decoder_comparison(Implementation impl_a, Implementation impl_b);
};

#endif // RS_COMPARISON_TEST_H
