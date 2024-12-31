#include <cstring>
#include <iostream>
#include "RsComparisonTest.h"

RsComparisonTest::RsComparisonTest(const RsTestConfig& cfg)
    : config(cfg), tt(cfg.getTT()) {
    prepare_test_data();
}

void RsComparisonTest::prepare_test_data() {
    GF data[MAX_KK];
    GF bb[2*MAX_TT];

    // Fill data with test pattern
    for (int i = 0; i < MAX_KK; i++) {
        data[i] = i & 0xFF;
    }

    // Initial encoding
    RS_ENCODER rs0(tt);
    rs0.RSEncode(data, bb);

    // Copy parity and data into test data
    const int kk = nn - 2*tt;
    memcpy(test_data, bb, 2*tt);        // Copy parity
    memcpy(test_data + nn - kk, data, kk);  // Copy data

    // Inject errors
    for (int i = 0; i < tt/2; i++) {
        test_data[i*2] ^= 0xFF;  // Flip some bits
    }
}

RsComparisonTest::TestContext RsComparisonTest::run_encoder(Implementation impl, bool use_erasures) const {
    TestContext ctx;

    // Initialize verification framework with debug enabled
    RsVerification::init(config.getVerboseLevel() == Verbosity::Debug);
    RsVerification::setResults(&ctx.results);

    // Copy test data
    memcpy(ctx.recd, test_data, nn);

    // Verify initial state
    std::cout << "\nInitial received word state:\n";
    RsVerification::verify_received_word(ctx.recd, nn);
    uint32_t initial_hash = ctx.results.received_word_hash;
    std::cout << "Initial received word hash: 0x" << std::hex << initial_hash << std::dec << "\n";

    // Initialize erasure data
    ctx.no_eras = 0;
    if (use_erasures) {
        ctx.no_eras = tt/2;
        std::cout << "\nSetting up " << ctx.no_eras << " erasures at positions: ";
        for (int i = 0; i < ctx.no_eras; i++) {
            ctx.eras_pos[i] = i*2;
            std::cout << ctx.eras_pos[i] << " ";
        }
        std::cout << "\n";
    }

    // Clear verification results before decoder run
    ctx.results.clear();
    RsVerification::setResults(&ctx.results);

    switch(impl) {
        case Implementation::Original:
            {
                std::cout << "\nUsing original implementation with "
                          << (use_erasures ? "erasures" : "normal decoding") << "\n";
                RS_ENCODER encoder(tt);
                if (use_erasures) {
                    std::cout << "Before RSDecodeErasures, hash: ";
                    RsVerification::verify_received_word(ctx.recd, nn);
                    std::cout << "0x" << std::hex << ctx.results.received_word_hash << std::dec << "\n";

                    ctx.decode_result = encoder.RSDecodeErasures(ctx.recd,
                        ctx.eras_pos, ctx.no_eras);

                    std::cout << "After RSDecodeErasures, hash: ";
                    RsVerification::verify_received_word(ctx.recd, nn);
                    std::cout << "0x" << std::hex << ctx.results.received_word_hash << std::dec << "\n";
                    std::cout << "Decode result: " << ctx.decode_result << "\n";
                } else {
                    ctx.decode_result = encoder.RSDecode(ctx.recd);
                }
            }
            break;

        case Implementation::Reference:
	    {
                std::cout << "\nUsing reference implementation with "
                          << (use_erasures ? "erasures" : "normal decoding") << "\n";
                RS_ENCODER_REF encoder(tt);
                if (use_erasures) {
                    std::cout << "Before RSDecodeErasures, hash: ";
                    RsVerification::verify_received_word(ctx.recd, nn);
                    std::cout << "0x" << std::hex << ctx.results.received_word_hash << std::dec << "\n";

                    ctx.decode_result = encoder.RSDecodeErasures(ctx.recd, ctx.eras_pos, ctx.no_eras);

                    std::cout << "After RSDecodeErasures, hash: ";
                    RsVerification::verify_received_word(ctx.recd, nn);
                    std::cout << "0x" << std::hex << ctx.results.received_word_hash << std::dec << "\n";
                    std::cout << "Decode result: " << ctx.decode_result << "\n";
                } else {
                    ctx.decode_result = encoder.RSDecode(ctx.recd);
                }
            }
            break;

        case Implementation::Karn:
            {
                std::cout << "\nUsing Karn implementation with "
                          << (use_erasures ? "erasures" : "normal decoding") << "\n";
                RS_ENCODER_KARN karn(tt);
                if (use_erasures) {
                    std::cout << "Before RSDecodeErasures, hash: ";
                    RsVerification::verify_received_word(ctx.recd, nn);
                    std::cout << "0x" << std::hex << ctx.results.received_word_hash << std::dec << "\n";

                    ctx.decode_result = karn.RSDecodeErasures(ctx.recd,
                        ctx.eras_pos, ctx.no_eras);

                    std::cout << "After RSDecodeErasures, hash: ";
                    RsVerification::verify_received_word(ctx.recd, nn);
                    std::cout << "0x" << std::hex << ctx.results.received_word_hash << std::dec << "\n";
                    std::cout << "Decode result: " << ctx.decode_result << "\n";
                } else {
                    ctx.decode_result = karn.RSDecode(ctx.recd);
                }
            }
            break;
    }
    return ctx;
}

bool RsComparisonTest::compare_results(const TestContext& left,
                                     const TestContext& right,
                                     const char* test_name) const {
#if 0
    // For same-operation comparisons (between implementations)
    if (test_name && strstr(test_name, "vs original")) {
#endif
        // Results should match exactly
        bool results_match = RsVerification::compare_results(
            left.results,
            right.results,
            config.getVerboseLevel() == Verbosity::Debug
        );
        bool decode_match = (left.decode_result == right.decode_result);

        if (!results_match || !decode_match) {
            std::cout << "FAIL: " << test_name << " mismatch\n";
            if (!decode_match) {
                std::cout << "  Decode results differ: " << left.decode_result
                         << " vs " << right.decode_result << "\n";
            }
            return false;
        }
#if 0
    }
    // For RSDecode vs RSDecodeErasures comparisons
    else {
        // If both decoders agree (whether finding errors or not), that's success
        if (left.decode_result == right.decode_result) {
            if (config.getVerboseLevel() == Verbosity::Debug) {
                std::cout << "PASS: " << test_name << "\n";
                std::cout << "  Normal decode found " << left.decode_result << " errors\n";
                std::cout << "  Erasure decode found " << right.decode_result << " errors\n";
            }
        }
	else
	{
            std::cout << "FAIL: " << test_name << " - decoding failed\n";
            std::cout << "  Decode results: " << left.decode_result
                     << " vs " << right.decode_result << "\n";
            return false;
        }
    }
#endif

    std::cout << "PASS: " << test_name << "\n";
    return true;
}

RsComparisonTest::EncodingTestContext RsComparisonTest::run_encoder_test(bool use_reference) const {
    EncodingTestContext ctx;

    // Initialize test data
    const int kk = nn - 2*tt;
    std::cout << "\n" << (use_reference ? "Reference" : "Original") << " implementation:\n";
    for (int i = 0; i < kk; i++) {
        ctx.data[i] = i & 0xFF;  // Simple test pattern
    }
    std::cout << "\n";

    // Initialize verification
    RsVerification::setResults(&ctx.results);

    // Show generator polynomial
    if (use_reference) {
        RS_ENCODER_REF encoder(tt);
        encoder.RSEncode(ctx.data, ctx.bb);
    } else {
        RS_ENCODER encoder(tt);
        encoder.RSEncode(ctx.data, ctx.bb);
    }

    return ctx;
}

bool RsComparisonTest::compare_encoding_results(const EncodingTestContext& left,
                                              const EncodingTestContext& right,
                                              const char* test_name) const {
    // Compare parity bytes with detailed output
    bool parity_matches = true;
    std::cout << "\nParity byte comparison:\n";
    for (int i = 0; i < 2*tt; i++) {
        if (left.bb[i] != right.bb[i]) {
            std::cout << "Position " << i << ": Reference=0x" << std::hex
                     << (int)left.bb[i] << " Original=0x" << (int)right.bb[i]
                     << std::dec << "\n";
            parity_matches = false;
        }
    }

    bool success = true;
    if (!parity_matches) {
        std::cout << "\nFAIL: " << test_name << " - parity bytes differ (shown above)\n";
        success = false;
    }

    // Verify encoding integrity
    bool results_match = RsVerification::compare_results(
        left.results,
        right.results,
        config.getVerboseLevel() == Verbosity::Debug
    );

    if (!results_match) {
        std::cout << "FAIL: " << test_name << " - verification results differ\n";
        success = false;
    }

    if (success) std::cout << "PASS: " << test_name << "\n";
    return success;
}

bool RsComparisonTest::run_decoder_comparison(Implementation impl_a, Implementation impl_b) {
    // Maintain separate verification contexts
    TestContext ctx_a;
    TestContext ctx_b;

    // Set up initial state
    memcpy(ctx_a.recd, test_data, nn);
    memcpy(ctx_b.recd, test_data, nn);

    // Run first decoder with its own verification context
    RsVerification::setResults(&ctx_a.results);
    ctx_a.decode_result = run_decoder(impl_a, ctx_a.recd);

    // Run second decoder with its own context
    RsVerification::setResults(&ctx_b.results);
    ctx_b.decode_result = run_decoder(impl_b, ctx_b.recd);

    // Now we can compare the complete states
    std::string comparison_name = get_impl_name(impl_a) + " vs " + get_impl_name(impl_b);
    return compare_results(ctx_a, ctx_b, comparison_name.c_str());
}

// Similar structure for erasure decoder comparison
bool RsComparisonTest::run_erasure_decoder_comparison(Implementation impl_a, Implementation impl_b) {
    TestContext ctx_a;
    TestContext ctx_b;

    // Set up erasures identically for both contexts
    ctx_a.no_eras = tt/2;
    ctx_b.no_eras = tt/2;
    for (int i = 0; i < ctx_a.no_eras; i++) {
        ctx_a.eras_pos[i] = i*2;
        ctx_b.eras_pos[i] = i*2;
    }

    // Copy initial state
    memcpy(ctx_a.recd, test_data, nn);
    memcpy(ctx_b.recd, test_data, nn);

    // Run decoders with separate verification contexts
    RsVerification::setResults(&ctx_a.results);
    ctx_a.decode_result = run_erasure_decoder(impl_a, ctx_a);

    RsVerification::setResults(&ctx_b.results);
    ctx_b.decode_result = run_erasure_decoder(impl_b, ctx_b);

    std::string comparison_name = get_impl_name(impl_a) + " vs " + get_impl_name(impl_b);
    return compare_results(ctx_a, ctx_b, comparison_name.c_str());
}

bool RsComparisonTest::test_implementations_match() {
    std::cout << "\nComparing All Decoder Implementations:\n";

    bool success = true;

    // Compare Reference vs Original
    success &= run_erasure_decoder_comparison(Implementation::Reference, Implementation::Original);

    // Compare Karn vs Original
    success &= run_erasure_decoder_comparison(Implementation::Karn, Implementation::Original);

    // Compare Karn vs Reference
    success &= run_erasure_decoder_comparison(Implementation::Karn, Implementation::Reference);

    // Test Karn's internal consistency
    success &= run_erasure_decoder_comparison(Implementation::Karn, Implementation::Karn);

    return success;
}

bool RsComparisonTest::run() {
    bool success = true;

    success &= test_implementations_match();

    return success;
}
// Pull out decoder execution into its own function
int RsComparisonTest::run_decoder(Implementation impl, GF recd[nn]) {
    switch(impl) {
        case Implementation::Original:
            {
                RS_ENCODER encoder(tt);
                return encoder.RSDecode(recd);
            }
            break;

        case Implementation::Reference:
            {
                RS_ENCODER_REF encoder(tt);
                return encoder.RSDecode(recd);
            }
            break;

        case Implementation::Karn:
            {
                RS_ENCODER_KARN encoder(tt);
                return encoder.RSDecode(recd);
            }
            break;
    }
    return -1;
}

// And same for erasure decoding
int RsComparisonTest::run_erasure_decoder(Implementation impl, TestContext& ctx) {
    switch(impl) {
        case Implementation::Original:
            {
                RS_ENCODER encoder(tt);
                return encoder.RSDecodeErasures(ctx.recd, ctx.eras_pos, ctx.no_eras);
            }
            break;

        case Implementation::Reference:
            {
                RS_ENCODER_REF encoder(tt);
                return encoder.RSDecodeErasures(ctx.recd, ctx.eras_pos, ctx.no_eras);
            }
            break;

        case Implementation::Karn:
            {
                RS_ENCODER_KARN encoder(tt);
                return encoder.RSDecodeErasures(ctx.recd, ctx.eras_pos, ctx.no_eras);
            }
            break;
    }
    return -1;
}

// Helper to get implementation name for messages
std::string RsComparisonTest::get_impl_name(Implementation impl) {
    switch(impl) {
        case Implementation::Original: return "Original";
        case Implementation::Reference: return "Reference";
        case Implementation::Karn: return "Karn";
        default: return "Unknown";
    }
}
