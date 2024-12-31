#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstring>
#include "RsDecodingTest.h"

RsDecodingTest::RsDecodingTest(const RsTestConfig& cfg)
    : config(cfg),
      tt(cfg.getTT()),
      kk(nn - 2*tt),
      uniform_dist(0.0, 1.0),
      byte_dist(1, 255) {
    initialize();
}

void RsDecodingTest::initialize() {
    // Initialize random number generator
    rng.seed(config.getRandomSeed());

    // Create encoder
    encoder = new RS_ENCODER(tt);

    // Initialize verification
    RsVerification::init(config.getVerboseLevel() == Verbosity::Debug);
}

bool RsDecodingTest::process_codeword(const GF* codeword, const int nn_short, const int kk_short) {
    stats.total_codewords++;

    // Create received word with errors
    GF received[nn];
    int no_ch_errs = 0;

    // Inject errors based on error rate
    for (int i = 0; i < nn_short; i++) {
        double x = uniform_dist(rng);
        if (x < config.getErrorRate()) {
            GF error_byte = byte_dist(rng);
            received[i] = codeword[i] ^ error_byte;
            ++no_ch_errs;
        } else {
            received[i] = codeword[i];
        }
    }

    stats.total_errors_injected += no_ch_errs;
    stats.errors_by_count[no_ch_errs]++;

    if (config.getVerboseLevel() == Verbosity::Debug) {
        std::cout << "Channel caused " << no_ch_errs << " errors\n";
    }

    // Pad with zeros for decoding
    GF recd[nn];
    for (int i = 0; i < nn-kk; i++) {
        recd[i] = received[i];  // parity bytes
    }
    for (int i = nn-kk; i < nn-kk+kk_short; i++) {
        recd[i] = received[i];  // info bytes
    }
    for (int i = nn-kk+kk_short; i < nn; i++) {
        recd[i] = 0;  // zero padding
    }

    // Decode
    RsVerification::VerificationResults verify_results;
    RsVerification::setResults(&verify_results);

    int decode_flag = encoder->RSDecode(recd);

    if (config.getVerboseLevel() == Verbosity::Debug) {
        std::cout << "decode_rs() returned " << decode_flag << "\n";
    }

    // Check for errors
    bool error_flag = false;
    for (int i = 0; i < kk_short; i++) {
        if (recd[i+nn-kk] != codeword[i+nn-kk]) {
            if (no_ch_errs <= tt) {
                std::cout << "Position " << i+nn-kk
                         << " miscorrected " << std::hex
                         << (int)recd[i+nn-kk] << ","
                         << (int)codeword[i+nn-kk] << "="
                         << (int)(recd[i+nn-kk]^codeword[i+nn-kk])
                         << std::dec << ".\n";
            }
            error_flag = true;
        }
    }

    // Update statistics
    if (decode_flag == 0) {
        stats.errors_corrected++;
    }
    if (error_flag && no_ch_errs <= tt) {
        stats.decode_failures++;
    }
    if (no_ch_errs > tt) {
        stats.uncorrectable_errors++;
    }

    // Check for specific decode failure modes
    if (decode_flag == RS_ERROR_LAMDA_ERROR && no_ch_errs <= tt) {
        std::cout << no_ch_errs << " ch errs  <=  max # correctable errs but\n"
                  << "DECODER ERROR: deg(lambda) unequal to #roots\n";
        return false;
    }
    if (decode_flag == RS_ERROR_CHIEN_SEARCH && no_ch_errs <= tt) {
        std::cout << no_ch_errs << " ch errs  <= max # correctable errs but\n"
                  << "deg(lambda) > 2*tt\n";
        return false;
    }
    if (error_flag && no_ch_errs <= tt) {
        std::cout << "DECODER FAILED TO CORRECT " << no_ch_errs << " ERRORS\n";
        return false;
    }

    return true;
}

void RsDecodingTest::print_stats() const {
    if (config.getVerboseLevel() == Verbosity::Quiet) {
        return;
    }

    std::cout << "\nDecoding Summary:\n"
              << "Total codewords processed: " << stats.total_codewords << "\n"
              << "Total errors injected: " << stats.total_errors_injected
              << " (avg " << (double)stats.total_errors_injected / stats.total_codewords
              << " per codeword, max correctable: " << tt << ")\n\n"
              << "Error distribution:\n";

    for (int i = 0; i <= tt*2 && i < MAX_TT; i++) {
        if (stats.errors_by_count[i] > 0) {
            std::cout << i << " errors: " << stats.errors_by_count[i]
                     << " times ("
                     << (stats.errors_by_count[i] * 100) / stats.total_codewords
                     << "% " << (i <= tt ? "correctable" : "uncorrectable")
                     << ")\n";
        }
    }
}

bool RsDecodingTest::open_and_process_file(const char* filename) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        std::cerr << "Could not open " << filename << "\n";
        return false;
    }

    const int nn_short = 255;
    if ((nn_short < 2*tt) || (nn_short > nn)) {
        std::cerr << "Invalid entry " << nn_short << " for shortened length\n";
        fclose(fp);
        return false;
    }

    const int kk_short = kk - (nn-nn_short);

    bool success = true;
    for (int i = 0; i < config.getNumCodewords(); i++) {
        GF codeword[nn_short];
        if (fread(codeword, sizeof(GF), nn_short, fp) != nn_short) {
            std::cerr << "Error reading codeword\n";
            success = false;
            break;
        }

        if (!process_codeword(codeword, nn_short, kk_short)) {
            success = false;
            break;
        }
    }

    fclose(fp);
    return success;
}

bool RsDecodingTest::run() {
    bool success = open_and_process_file("cws.txt");
    print_stats();
    return success;
}
