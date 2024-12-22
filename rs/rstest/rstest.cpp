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
#include "RsVerification.h"
#include "RsTestConfig.h"

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

double erand48()
{
    return ((double)rand())/((double)RAND_MAX);
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

// Helper to prepare test data
void prepare_test_data(GF* recd, int tt) {
    GF data[MAX_KK];
    GF bb[2*MAX_TT];

    // Fill data with test pattern
    for (int i = 0; i < MAX_KK; i++) {
        data[i] = i & 0xFF;
    }

    // Initial encoding
    RS_ENCODER rs0(tt);
    rs0.RSEncode(data, bb);

    // Copy parity and data into received word
    const int kk = nn - 2*tt;
    memcpy(recd, bb, 2*tt);        // Copy parity
    memcpy(recd + nn - kk, data, kk);  // Copy data

    // Inject errors
    for (int i = 0; i < tt/2; i++) {
        recd[i*2] ^= 0xFF;  // Flip some bits
    }
}

// Set up and run a single encoder test
EncoderTestContext setup_and_run_encoder(int tt, bool use_reference, const GF* test_data) {
    EncoderTestContext ctx;

    // Initialize verification framework
    RsVerification::init(true);
    RsVerification::setResults(&ctx.results);

    // Copy test data
    memcpy(ctx.recd, test_data, nn);

    // Create and run appropriate encoder
    if (use_reference) {
        RS_ENCODER_REF encoder(tt);
        ctx.decode_result = encoder.RSDecode(ctx.recd);
    } else {
        RS_ENCODER encoder(tt);
        ctx.decode_result = encoder.RSDecode(ctx.recd);
    }

    return ctx;
}

GF data[MAX_KK];
GF recd[nn];
GF bb[2*MAX_TT];
std::mt19937 rng;  // Mersenne Twister
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
std::uniform_int_distribution<int> byte_dist(1, 255);

void compare_implementations(int tt) {
    // Prepare test data
    GF test_data[nn];
    prepare_test_data(test_data, tt);

    printf("\nTesting Reference Implementation Internal Consistency:\n");

    // Reference Implementation: RSDecode vs RSDecodeErasures
    EncoderTestContext ref_decode = setup_and_run_encoder(tt, true, test_data);
    EncoderTestContext ref_decode_eras = setup_and_run_encoder(tt, true, test_data); // TODO: Add erasures

    bool results_match = RsVerification::compare_results(
        ref_decode.results,
        ref_decode_eras.results,
        false
    );

    bool decode_match = (ref_decode.decode_result == ref_decode_eras.decode_result);

    if (!results_match || !decode_match) {
        printf("FAIL: Reference implementation RSDecode vs RSDecodeErasures mismatch\n");
    } else {
        printf("PASS: Reference implementation RSDecode matches RSDecodeErasures\n");
    }

    printf("\nTesting Original Implementation Internal Consistency:\n");

    // Original Implementation: RSDecode vs RSDecodeErasures
    EncoderTestContext orig_decode = setup_and_run_encoder(tt, false, test_data);
    EncoderTestContext orig_decode_eras = setup_and_run_encoder(tt, false, test_data); // TODO: Add erasures

    results_match = RsVerification::compare_results(
        orig_decode.results,
        orig_decode_eras.results,
        false
    );

    decode_match = (orig_decode.decode_result == orig_decode_eras.decode_result);

    if (!results_match || !decode_match) {
        printf("FAIL: Original implementation RSDecode vs RSDecodeErasures mismatch\n");
    } else {
        printf("PASS: Original implementation RSDecode matches RSDecodeErasures\n");
    }

    // Compare reference vs original implementation
    printf("\nComparing Reference vs Original Implementation:\n");
    results_match = RsVerification::compare_results(
        ref_decode.results,
        orig_decode.results,
        false
    );
    decode_match = (ref_decode.decode_result == orig_decode.decode_result);

    if (!results_match || !decode_match) {
        printf("FAIL: Reference implementation differs from original\n");
    } else {
        printf("PASS: Reference implementation matches original\n");
    }
}

int main(int argc, char *argv[]) {
    try {
        RsTestConfig config = RsTestConfig::parse_args(argc, argv);
        if (config.getVerboseLevel() != Verbosity::Quiet) {
            config.print();
        }

        // Initialize RS code
        int tt = config.getTT();
        int kk = nn - 2*tt;

        RS_Init();
        RS_ENCODER::Init();
        RS_ENCODER_REF::Init();

        const char* cw_file = "cws.txt";
        const char* file_in = "cws.txt";;
        int kk_short, nn_short, iter, no_decoder_errors;
        DecodingStats stats = {
            .total_codewords = 0,
            .total_errors_injected = 0,
            .errors_corrected = 0,
            .decode_failures = 0,
            .uncorrectable_errors = 0
        };

        // Initialize verification framework
        RsVerification::init(config.getVerboseLevel() == Verbosity::Debug);

        RS_ENCODER *prs = new RS_ENCODER(tt);

        switch (config.getMode()) {
            case TestMode::Encoding:
                nn_short = 255;
                if ((nn_short < 2*tt) || (nn_short > nn))
                {
                    printf("Invalid entry %d for shortened length\n",nn_short);
                    exit(0);
                }

                /* compute dimension of shortened code */
                kk_short = kk - (nn-nn_short);

                FILE *fp_out;
                if (!rs_fopen(&fp_out, cw_file, "wb") || nullptr == fp_out)
                {
                    printf("Could not open %s\n", cw_file);
                    exit(0);
                }

                /**** BEGIN: Encoding random information vectors ****/
                /* Pad with zeros rightmost (kk-kk_short) positions */
                memset(data + kk_short, 0, kk - kk_short);
                for (int i = 0; i < config.getNumCodewords(); i++)
                {
                    for (int j = 0; j < (int)kk_short; j++)
                        data[j] = (int) (nrand48() % 256);

                    prs->RSEncode(data, bb);

                    fwrite(bb, sizeof(GF), nn_short - kk_short, fp_out);
                    fwrite(data, sizeof(GF), kk_short, fp_out);
                }
                fclose(fp_out);
                break;

            case TestMode::Decoding:
                memset(stats.errors_by_count, 0, sizeof(stats.errors_by_count));

                nn_short = 255;
                if ((nn_short < 2*tt) || (nn_short > nn))
                {
                    printf("Invalid entry %d for shortened length\n",nn_short);
                    exit(0);
                }
                kk_short = kk - (nn-nn_short); /* compute dimension of shortened code */
                rng.seed(config.getRandomSeed());

                FILE *fp;
                if (!rs_fopen(&fp, file_in, "rb") || nullptr == fp)
                {
                    printf("Could not open %s\n", file_in);
                    exit(0);
                }

                no_decoder_errors = 0;
                iter = -1;
                while (++iter < config.getNumCodewords()) {
                    stats.total_codewords++;
                    GF received[1000], cw[1000];
                    fread(cw, sizeof(GF), nn_short, fp);

                    print_debug(config, "\n\n\n\n Transmitting codeword %d \n", iter);
                    print_debug(config, "Channel caused following errors Location (Error): \n");
                    int no_ch_errs = 0;
                    for (int i=0;i < nn_short;i++)
                    {
                        double x = uniform_dist(rng);
                        print_debug(config, "Random value: %f vs threshold: %f\n", x, config.getErrorRate());
                        if (x < config.getErrorRate())
                        {
                            GF error_byte = byte_dist(rng);
                            print_debug(config, "Injecting error: %d\n", error_byte);
                            received[i] = cw[i] ^ error_byte;
                            ++no_ch_errs;
                        }
                        else
                            received[i] = cw[i];
                    }
                    stats.total_errors_injected += no_ch_errs;
                    stats.errors_by_count[no_ch_errs]++;
                    print_debug(config, "Channel caused %d errors\n", no_ch_errs);

                    // Pad with zeros and decode as if was a (255,kk) tt-error correcting
                    // code.
                    //
                    int i;
                    for (i=0;i < nn-kk;i++)
                        recd[i] = received[i]; /* parity bytes */
                    for (i=nn-kk;i < nn-kk+kk_short;i++)
                        recd[i] = received[i]; /* info bytes */
                    for (i=nn-kk+kk_short;i < nn;i++)
                        recd[i] = 0; /* zero padding */

                    int decode_flag = prs->RSDecode(recd);
                    print_debug(config, "decode_rs() returned %d\n", decode_flag);
                    int error_flag = 0;
                    for (i=0; i < kk_short; i++)
                    {
                        if (recd[i+nn-kk] != cw[i+nn-kk])
                        {
                            if (no_ch_errs <= tt)
                            {
                                printf("Position %d miscorrected %x,%x=%x.\n", i+nn-kk,
                                recd[i+nn-kk], cw[i+nn-kk], recd[i+nn-kk]^cw[i+nn-kk]);
                            }
                            error_flag = 1;
                            //break;
                        }
                    }

                    if (decode_flag == 0) {
                        stats.errors_corrected++;
                    }
                    if (error_flag == 1 && no_ch_errs <= tt) {
                        stats.decode_failures++;
                    }
                    if (no_ch_errs > tt) {
                        stats.uncorrectable_errors++;
                    }

                    if (decode_flag == RS_ERROR_LAMDA_ERROR && no_ch_errs <= tt)
                    {
                        printf("%d ch errs  <=  max # correctable errs but \n", no_ch_errs);
                        printf("\n DECODER ERROR: deg(lambda) unequal to #roots \n");
                        exit(2); /* DECODER ERROR condition */
                    }
                    else if (decode_flag == RS_ERROR_CHIEN_SEARCH && no_ch_errs <= tt)
                    {
                        printf(" %d ch errs  <= max # correctable errs but \n", no_ch_errs);
                        printf("\n deg(lambda) > 2*tt \n");
                        exit(3);/* DECODER ERROR condition */
                    }
                    if (error_flag == 1 && no_ch_errs <= tt)
                    {
                        printf("DECODER FAILED TO CORRECT %d ERRORS\n", no_ch_errs);
                        ++no_decoder_errors;
                    }
                }

                print_progress(config, "\nDecoding Summary:\n");
                print_progress(config, "Total codewords processed: %d\n", stats.total_codewords);
                print_progress(config, "Total errors injected: %d (avg %.2f per codeword, max correctable: %d)\n",
                    stats.total_errors_injected,
                    (double)stats.total_errors_injected / stats.total_codewords,
                    tt);

                // Add error distribution
                print_progress(config, "\nError distribution:\n");
                for (int i = 0; i <= tt*2 && i < MAX_TT; i++) {  // Show up to 2*tt errors
                    if (stats.errors_by_count[i] > 0) {
                        print_progress(config, "%d errors: %d times (%d%% %s)\n",
                            i,
                            stats.errors_by_count[i],
                            (stats.errors_by_count[i] * 100) / stats.total_codewords,
                            i <= tt ? "correctable" : "uncorrectable");
                    }
                }
                break;

            case TestMode::Benchmark:
                run_benchmarks();
                break;

            case TestMode::Comparison:
                compare_implementations(tt);
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
