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

RsComparisonTest::TestContext RsComparisonTest::run_encoder(bool use_reference, bool use_erasures) const {
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

    // Run appropriate encoder version
    if (use_reference) {
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
    } else {
        std::cout << "\nUsing original implementation with "
                  << (use_erasures ? "erasures" : "normal decoding") << "\n";
        RS_ENCODER encoder(tt);
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

    return ctx;
}

bool RsComparisonTest::compare_results(const TestContext& left,
                                     const TestContext& right,
                                     const char* test_name) const {
    // For same-operation comparisons (between implementations)
    if (test_name && strstr(test_name, "vs original")) {
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
    }
    // For RSDecode vs RSDecodeErasures comparisons
    else {
        // Results should show successful correction but hashes will differ
        if (left.decode_result <= 0 || right.decode_result <= 0) {
            std::cout << "FAIL: " << test_name << " - decoding failed\n";
            std::cout << "  Decode results: " << left.decode_result
                     << " vs " << right.decode_result << "\n";
            return false;
        }

        // For debugging
        if (config.getVerboseLevel() == Verbosity::Debug) {
            std::cout << "SUCCESS: " << test_name << "\n";
            std::cout << "  Normal decode found " << left.decode_result << " errors\n";
            std::cout << "  Erasure decode corrected " << right.decode_result
                     << " positions\n";
        }
    }

    std::cout << "PASS: " << test_name << "\n";
    return true;
}

bool RsComparisonTest::test_reference_implementation() const {
    std::cout << "\nTesting Reference Implementation Internal Consistency:\n";

    TestContext ref_decode = run_encoder(true, false);
    TestContext ref_decode_eras = run_encoder(true, true);

    return compare_results(ref_decode, ref_decode_eras,
                         "Reference implementation RSDecode vs RSDecodeErasures");
}

bool RsComparisonTest::test_original_implementation() const {
    std::cout << "\nTesting Original Implementation Internal Consistency:\n";

    TestContext orig_decode = run_encoder(false, false);
    TestContext orig_decode_eras = run_encoder(false, true);

    return compare_results(orig_decode, orig_decode_eras,
                         "Original implementation RSDecode vs RSDecodeErasures");
}

bool RsComparisonTest::test_implementations_match() const {
    std::cout << "\nComparing Reference vs Original Implementation:\n";

    TestContext ref_decode = run_encoder(true, false);
    TestContext orig_decode = run_encoder(false, false);

    return compare_results(ref_decode, orig_decode,
                         "Reference implementation vs original");
}

bool RsComparisonTest::run() {
    bool success = true;

    success &= test_reference_implementation();
    success &= test_original_implementation();
    success &= test_implementations_match();

    return success;
}
