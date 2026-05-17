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
    const int tt_values[] = { 1,2,3,4,5,6,7,8,16,24,32,40,48,56,64 };
    const int iterations = 100000;  // Adjustable based on speed needed

    // Pre-allocate our buffers
    GF data[MAX_KK];
    GF recd[nn];
    GF bb[2 * MAX_TT];

    // Fill data buffer with some pattern
    for (int i = 0; i < MAX_KK; i++) {
        data[i] = i & 0xFF;
    }

    printf("//                Auto             Direct (scalar)    Direct (NEON)        Direct (AVX2)       Direct (AVX512)\n");
    printf("// tt  Encode  w/err   clean      w/err   clean      w/err   clean      w/err   clean      w/err   clean\n");

    for (const int tt_val : tt_values) {
        const int b0 = (nn - 2 * tt_val + 1) / 2;
        auto codec_auto    = RS_FACTORY::instance().create_codec(tt_val, b0, RS_FACTORY::DecoderKind::Auto);
        auto codec_direct  = RS_FACTORY::instance().create_codec(tt_val, b0, RS_FACTORY::DecoderKind::Direct);
        auto codec_neon    = RS_FACTORY::instance().create_codec(tt_val, b0, RS_FACTORY::DecoderKind::DirectNeon);
        auto codec_avx2    = RS_FACTORY::instance().create_codec(tt_val, b0, RS_FACTORY::DecoderKind::DirectAvx2);
        auto codec_avx512  = RS_FACTORY::instance().create_codec(tt_val, b0, RS_FACTORY::DecoderKind::DirectAvx512);

        // Correctness check: all three decoders must agree.
        {
            const int kk_check = nn - 2 * tt_val;
            GF golden[nn], probe[nn];
            codec_auto->RSEncode(data, bb);
            memcpy(golden, bb, 2 * tt_val);
            memcpy(golden + nn - kk_check, data, kk_check);
            memcpy(probe, golden, nn);
            for (int i = 0; i < tt_val / 2; i++) {
                probe[i * 2] ^= 0xFF;
            }
            GF probe_a[nn], probe_d[nn], probe_n[nn], probe_x[nn], probe_5[nn];
            memcpy(probe_a, probe, nn);
            memcpy(probe_d, probe, nn);
            memcpy(probe_n, probe, nn);
            memcpy(probe_x, probe, nn);
            memcpy(probe_5, probe, nn);
            int ra = codec_auto->RSDecode(probe_a);
            int rd = codec_direct->RSDecode(probe_d);
            int rn = codec_neon->RSDecode(probe_n);
            int rx = codec_avx2->RSDecode(probe_x);
            int r5 = codec_avx512->RSDecode(probe_5);
            const bool buf_ok = (memcmp(probe_a, golden, nn) == 0)
                             && (memcmp(probe_d, golden, nn) == 0)
                             && (memcmp(probe_n, golden, nn) == 0)
                             && (memcmp(probe_x, golden, nn) == 0)
                             && (memcmp(probe_5, golden, nn) == 0);
            if (ra != rd || ra != rn || ra != rx || ra != r5 || !buf_ok) {
                fprintf(stderr, "ERROR: decoder disagreement at tt=%d (ra=%d rd=%d rn=%d rx=%d r5=%d, bufs=%d)\n",
                    tt_val, ra, rd, rn, rx, r5, buf_ok);
                return;
            }
        }
        const int kk = nn - 2 * tt_val;
        double encode_time = 0.0;
        double decode_auto_clean = 0.0,   decode_auto_errors = 0.0;
        double decode_direct_clean = 0.0, decode_direct_errors = 0.0;
        double decode_neon_clean = 0.0,    decode_neon_errors = 0.0;
        double decode_avx2_clean = 0.0,    decode_avx2_errors = 0.0;
        double decode_avx512_clean = 0.0,  decode_avx512_errors = 0.0;

        // Time encoding (uses Auto codec; the encoder is identical either way)
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_auto->RSEncode(data, bb);
        }
        auto end = std::chrono::high_resolution_clock::now();
        encode_time = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // --- Auto decoder ---
        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_auto->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_auto_clean = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        for (int i = 0; i < tt_val / 2; i++) {
            recd[i * 2] ^= 0xFF;
        }
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_auto->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_auto_errors = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // --- Direct decoder ---
        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_direct->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_direct_clean = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        for (int i = 0; i < tt_val / 2; i++) {
            recd[i * 2] ^= 0xFF;
        }
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_direct->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_direct_errors = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // --- NEON decoder ---
        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_neon->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_neon_clean = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        for (int i = 0; i < tt_val / 2; i++) {
            recd[i * 2] ^= 0xFF;
        }
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_neon->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_neon_errors = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // --- AVX2 decoder ---
        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_avx2->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_avx2_clean = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        for (int i = 0; i < tt_val / 2; i++) {
            recd[i * 2] ^= 0xFF;
        }
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_avx2->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_avx2_errors = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        // --- AVX-512 decoder ---
        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_avx512->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_avx512_clean = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        memcpy(recd, bb, 2 * tt_val);
        memcpy(recd + nn - kk, data, kk);
        for (int i = 0; i < tt_val / 2; i++) {
            recd[i * 2] ^= 0xFF;
        }
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            codec_avx512->RSDecode(recd);
        }
        end = std::chrono::high_resolution_clock::now();
        decode_avx512_errors = std::chrono::duration<double, std::micro>(end - start).count() / iterations;

        printf("//%3d  %6.2f  %7.2f %7.2f    %7.2f %7.2f    %7.2f %7.2f    %7.2f %7.2f    %7.2f %7.2f\n",
            tt_val, encode_time,
            decode_auto_errors,    decode_auto_clean,
            decode_direct_errors,  decode_direct_clean,
            decode_neon_errors,    decode_neon_clean,
            decode_avx2_errors,    decode_avx2_clean,
            decode_avx512_errors,  decode_avx512_clean);

        // --- Phase breakdown on the fastest ("Auto") decoder for the errored case ---
        // This is the key data for deciding the next target while squeezing.
        {
            // Re-inject the same error pattern we used for timing
            memcpy(recd, bb, 2 * tt_val);
            memcpy(recd + nn - kk, data, kk);
            for (int i = 0; i < tt_val / 2; i++) {
                recd[i * 2] ^= 0xFF;
            }

            auto* dec = const_cast<RS_DECODER_BASE*>(codec_auto->get_decoder());
            auto prof = dec->profile_decode(recd);
            if (prof.success || prof.errors_found > 0) {
                printf("//     [Auto phases @ tt=%d]  synd=%.2f  BM=%.2f  Chien=%.2f  Forney=%.2f  (total=%.2f)  errs=%d\n",
                       tt_val,
                       prof.syndrome_us,
                       prof.berlekamp_massey_us,
                       prof.chien_search_us,
                       prof.forney_us,
                       prof.total_us,
                       prof.errors_found);
            }
        }
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
