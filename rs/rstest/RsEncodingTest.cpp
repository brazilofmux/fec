#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstring>
#include "RsEncodingTest.h"

RsEncodingTest::RsEncodingTest(const RsTestConfig& cfg)
    : config(cfg),
      tt(cfg.getTT()),
      kk(nn - 2*tt),
      codec(nullptr) {
    initialize();
}

void RsEncodingTest::initialize() {
    // Initialize random number generator
    rng.seed(config.getRandomSeed());

    // Create encoder
    codec = RS_FACTORY::instance().create_codec(tt, (kk+1)/2);

    // Initialize verification framework
    RsVerification::init(config.getVerboseLevel() == Verbosity::Debug);
}

bool RsEncodingTest::validate_shortened_params(int nn_short, int kk_short) const {
    if ((nn_short < 2*tt) || (nn_short > nn)) {
        std::cerr << "Invalid entry " << nn_short << " for shortened length\n";
        return false;
    }

    if (kk_short <= 0 || kk_short > kk) {
        std::cerr << "Invalid shortened message length " << kk_short << "\n";
        return false;
    }

    return true;
}

bool RsEncodingTest::encode_and_save_codewords(const char* filename) {
    FILE* fp_out = fopen(filename, "wb");
    if (!fp_out) {
        std::cerr << "Could not open " << filename << " for writing\n";
        return false;
    }

    const int nn_short = 255;  // Fixed shortened length
    const int kk_short = kk - (nn-nn_short);

    if (!validate_shortened_params(nn_short, kk_short)) {
        fclose(fp_out);
        return false;
    }

    // Prepare data buffer
    GF data[MAX_KK];
    GF bb[2*MAX_TT];

    // Pad with zeros rightmost (kk-kk_short) positions
    memset(data + kk_short, 0, kk - kk_short);

    RsVerification::VerificationResults verify_results;
    RsVerification::setResults(&verify_results);

    // Generate and write codewords
    for (int i = 0; i < config.getNumCodewords(); i++) {
        // Generate random data
        for (int j = 0; j < kk_short; j++) {
            data[j] = rng() % 256;
        }

        // Encode
        codec->RSEncode(data, bb);

        // Write parity and data
        if (fwrite(bb, sizeof(GF), nn_short - kk_short, fp_out) != static_cast<size_t>(nn_short - kk_short) ||
            fwrite(data, sizeof(GF), kk_short, fp_out) != static_cast<size_t>(kk_short)) {
            std::cerr << "Error writing codeword " << i << "\n";
            fclose(fp_out);
            return false;
        }

        if (config.getVerboseLevel() == Verbosity::Debug) {
            std::cout << "Generated codeword " << i + 1 << " of "
                     << config.getNumCodewords() << "\n";
            std::cout.flush();
        }
    }

    if (config.getVerboseLevel() == Verbosity::Debug) {
        std::cout << "\nEncoding complete\n";
    }

    fclose(fp_out);
    return true;
}

bool RsEncodingTest::run() {
    return encode_and_save_codewords("cws.txt");
}
