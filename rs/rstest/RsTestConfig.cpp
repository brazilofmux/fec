#include <iostream>
#include <cstdlib>
#include <cstring>
#include "RsTestConfig.h"
#include "rs.h"  // For MIN_TT and MAX_TT

void RsTestConfig::print_usage() {
    std::cout << "Usage: rstest [-e|-d|-c|-b] [tt] [-v level] [-n count] [-s seed] [-r rate]\n"
              << "  -e tt : encode mode with error correction capacity tt\n"
              << "  -d tt : decode mode with error correction capacity tt\n"
              << "  -c tt : compare mode with error correction capacity tt\n"
              << "  -b    : benchmarking mode\n"
              << "  -v level    : verbosity (0=quiet, 1=normal, 2=debug)\n"
              << "  -n count    : number of codewords (default: 10000)\n"
              << "  -s seed     : random seed (default: 1093)\n"
              << "  -r rate     : error injection rate (default: 0.01)\n";
}

RsTestConfig RsTestConfig::parse_args(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        throw ConfigurationError("Insufficient arguments");
    }

    RsTestConfig config;
    bool needs_tt = false;

    // Parse mode first
    if (argv[1][0] != '-') {
        throw ConfigurationError("First argument must be mode (-e, -d, -c, or -b)");
    }

    switch (argv[1][1]) {
        case 'e':
        case 'E':
            config.mode = TestMode::Encoding;
            needs_tt = true;
            break;

        case 'd':
        case 'D':
            config.mode = TestMode::Decoding;
            needs_tt = true;
            break;

        case 'c':
        case 'C':
            config.mode = TestMode::Comparison;
            needs_tt = true;
            break;

        case 'b':
        case 'B':
            config.mode = TestMode::Benchmark;
            break;

        default:
            throw ConfigurationError("Invalid mode specified");
    }

    // Parse tt if needed
    if (needs_tt) {
        if (argc < 3) {
            throw ConfigurationError("tt parameter required for selected mode");
        }
        config.tt = std::atoi(argv[2]);
        if (config.tt < MIN_TT || config.tt > MAX_TT) {
            throw ConfigurationError("tt must be between " + 
                                   std::to_string(MIN_TT) + " and " + 
                                   std::to_string(MAX_TT));
        }
    }

    // Parse optional arguments
    for (int i = needs_tt ? 3 : 2; i < argc; i++) {
        if (argv[i][0] != '-' || strlen(argv[i]) != 2) {
            throw ConfigurationError("Invalid option: " + std::string(argv[i]));
        }

        if (i + 1 >= argc) {
            throw ConfigurationError("Missing value for option: " + std::string(argv[i]));
        }

        switch (argv[i][1]) {
            case 'v': {
                int v = std::atoi(argv[++i]);
                if (v < 0 || v > 2) {
                    throw ConfigurationError("Verbosity must be between 0 and 2");
                }
                config.verbose_level = static_cast<Verbosity>(v);
                break;
            }

            case 'n': {
                int nCodeWords = std::atoi(argv[++i]);
                if (nCodeWords <= 0) {
                    throw ConfigurationError("Number of codewords must be positive");
                }
                config.num_codewords = nCodeWords;
                break;
            }

            case 's':
                config.random_seed = static_cast<unsigned int>(std::atoi(argv[++i]));
                break;

            case 'r': {
                double r = std::atof(argv[++i]);
                if (r <= 0.0 || r >= 1.0) {
                    throw ConfigurationError("Error rate must be between 0 and 1");
                }
                config.error_rate = r;
                break;
            }

            default:
                throw ConfigurationError("Unknown option: -" + 
                                       std::string(1, argv[i][1]));
        }
    }

    // Set default error rate if not specified
    if (config.error_rate == 0.0 && config.mode == TestMode::Decoding) {
        config.error_rate = 0.01;  // Default 1% error rate
    }

    return config;
}

void RsTestConfig::print() const {
    std::cout << "RS Test Configuration:\n"
              << "  Mode: ";
    
    switch (mode) {
        case TestMode::Encoding:   std::cout << "Encoding\n"; break;
        case TestMode::Decoding:   std::cout << "Decoding\n"; break;
        case TestMode::Comparison: std::cout << "Comparison\n"; break;
        case TestMode::Benchmark:  std::cout << "Benchmark\n"; break;
        default: std::cout << "Unknown\n"; break;
    }

    if (mode != TestMode::Benchmark) {
        std::cout << "  Error correction capacity (tt): " << tt << "\n";
    }

    std::cout << "  Number of codewords: " << num_codewords << "\n"
              << "  Verbosity level: ";
    
    switch (verbose_level) {
        case Verbosity::Quiet:  std::cout << "Quiet\n"; break;
        case Verbosity::Normal: std::cout << "Normal\n"; break;
        case Verbosity::Debug:  std::cout << "Debug\n"; break;
    }

    std::cout << "  Random seed: " << random_seed << "\n";
    
    if (mode == TestMode::Decoding) {
        std::cout << "  Error rate: " << error_rate << "\n";
    }
    
    std::cout << "  Verify correction: " << (verify_correction ? "Yes" : "No") << "\n";
}
