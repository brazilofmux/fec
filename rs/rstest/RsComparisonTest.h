#ifndef RS_COMPARISON_TEST_H
#define RS_COMPARISON_TEST_H

#include "rs.h"
#include "rsref.h"
#include "RsTestConfig.h"
#include "RsVerification.h"

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
    TestContext run_encoder(bool use_reference, bool use_erasures = false) const;
    bool compare_results(const TestContext& left, const TestContext& right,
                        const char* test_name) const;

    // New test methods
    bool test_reference_implementation() const;
    bool test_original_implementation() const;
    bool test_implementations_match() const;
};

#endif // RS_COMPARISON_TEST_H
