// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All rights reserved.
// Stephen V. Dennis <brazilofmux@gmail.com>
//

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rs.h"
#include "rsref.h"
#include "RsComparisonTest.h"
#include "RsDecodingTest.h"
#include "RsEncodingTest.h"
#include "RsTestConfig.h"
#include "RsVerification.h"
#include "rs_factory.h"

struct DecodingStats {
    int total_codewords;
    int total_errors_injected;
    int errors_corrected;
    int decode_failures;
    int errors_by_count[MAX_TT + 1];  // Distribution of error counts
    int uncorrectable_errors;         // When errors > tt
};

// Replace old test context with new one using RsVerification
struct EncoderTestContext {
    RsVerification::VerificationResults results;
    GF recd[nn];
    int decode_result;
    int eras_pos[2*MAX_TT];
    int no_eras;
};

void print_progress(const RsTestConfig& config, const char* msg, ...) {
    if (config.getVerboseLevel() == Verbosity::Quiet) {
        return;
    }

    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
}

void print_debug(const RsTestConfig& config, const char* msg, ...) {
    if (config.getVerboseLevel() != Verbosity::Debug) {
        return;
    }

    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
}

long int nrand48()
{
    return rand();
}

bool rs_fopen(FILE** pFile, const char* filename, const char* mode)
{
    if (pFile)
    {
        *pFile = nullptr;
        if (nullptr != filename
            && nullptr != mode)
        {
#if !defined(__INTEL_COMPILER) && (_MSC_VER >= 1400)
            // 1400 is Visual C++ 2005
            //
            return (fopen_s(pFile, filename, mode) == 0);
#else
            * pFile = fopen((char*)filename, (char*)mode);
            if (nullptr != *pFile)
            {
                return true;
            }
#endif // WINDOWS_FILES
        }
    }
    return false;
}

void run_benchmarks() {
    const int tt_values[] = {1,2,3,4,5,6,7,8,16,24,32,40,48,56,64};
    const int iterations = 100000;  // Adjustable based on speed needed

    // Pre-allocate our buffers
    GF data[MAX_KK];
    GF recd[nn];
    GF bb[2*MAX_TT];

    // Fill data buffer with some pattern
    for (int i = 0; i < MAX_KK; i++) {
        data[i] = i & 0xFF;
    }

    printf("//             ---Errors Only---  -Errors+Erasures-\n");
    printf("//             ----Decoding-----  ----Decoding-----\n");
    printf("// tt  Encode  w/errors  without  w/errors  without\n");

    for (const int tt_val : tt_values) {
        RS_ENCODER rs(tt_val);
        const int kk = nn - 2*tt_val;
        double encode_time = 0.0;
        double decode_time_clean = 0.0;
        double decode_time_errors = 0.0;

        // Time encoding
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            rs.RSEncode(data, bb);
        }
        auto end = std::chrono::high_resolution_clock::now();
        encode_time = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // Time decoding without errors
        memcpy(recd, bb, 2*tt_val);        // Copy parity
        memcpy(recd + nn - kk, data, kk);  // Copy data

        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            rs.RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_time_clean = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // Time decoding with errors
        memcpy(recd, bb, 2*tt_val);        // Copy parity
        memcpy(recd + nn - kk, data, kk);  // Copy data
        // Inject tt_val/2 errors at fixed positions
        for (int i = 0; i < tt_val/2; i++) {
            recd[i*2] ^= 0xFF;
        }

        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            rs.RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_time_errors = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        printf("%3d  %6.2f  %8.2f  %8.2f    -----    -----\n",
               tt_val, encode_time, decode_time_errors, decode_time_clean);
    }
}

int main(int argc, char *argv[]) {
    try {
        RsTestConfig config = RsTestConfig::parse_args(argc, argv);
        if (config.getVerboseLevel() != Verbosity::Quiet) {
            config.print();
        }

        RS_Init();
        RS_ENCODER::Init();

        // Initialize verification framework
        RsVerification::init(config.getVerboseLevel() == Verbosity::Debug);

        switch (config.getMode()) {
            case TestMode::Encoding:
                {
                    RsEncodingTest test(config);
                    bool success = test.run();
                    return success ? 0 : 1;
                }
                break;

            case TestMode::Decoding:
                {
                    RsDecodingTest test(config);
                    bool success = test.run();
                    return success ? 0 : 1;
                }
                break;

            case TestMode::Benchmark:
                run_benchmarks();
                break;

            case TestMode::Comparison:
                {
                    RsComparisonTest test(config);
                    bool success = test.run();
                    return success ? 0 : 1;
                }
                break;

            case TestMode::Unknown:
                std::cerr << "Unknown test mode." << std::endl;
                break;
        }

    } catch (const ConfigurationError& e) {
        std::cerr << "Configuration error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

void example_usage() {
    // Initialize test data
    GF data[MAX_KK] = { 0 };
    for (int i = 0; i < MAX_KK; i++) {
        data[i] = i & 0xFF;
    }
    GF parity[2 * MAX_TT] = { 0 };

    // Get our two test codecs
    auto codec1 = RS_FACTORY::instance().create_codec(1, 1);
    auto codec2 = RS_FACTORY::instance().create_specialized_codec(1, 1);

    // Make full test blocks
    GF block1[nn] = { 0 };
    GF block2[nn] = { 0 };
    memcpy(block1, data, MAX_KK);
    memcpy(block2, data, MAX_KK);

    // Encode both blocks
    codec1->RSEncode(data, parity);
    memcpy(block1 + MAX_KK, parity, 2 * MAX_TT);

    codec2->RSEncode(data, parity);
    memcpy(block2 + MAX_KK, parity, 2 * MAX_TT);

    // Compare results
    if (memcmp(block1, block2, nn) == 0) {
        printf("Both encoders give same result\n");
    }

    // Now test decoders
    // Introduce an error in both blocks
    block1[10] ^= 0x01;  // Flip a bit
    block2[10] ^= 0x01;  // Same error

    int result1 = codec1->RSDecode(block1);
    int result2 = codec2->RSDecode(block2);

    if (result1 == result2 && memcmp(block1, block2, nn) == 0) {
        printf("Both decoders give same result\n");
    }
}

