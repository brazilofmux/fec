#ifndef RS_TEST_CONFIG_H
#define RS_TEST_CONFIG_H

#include <string>
#include <stdexcept>

enum class TestMode {
    Unknown,
    Encoding,
    Decoding,
    Comparison,
    Benchmark
};

enum class Verbosity {
    Quiet,      // Summary only
    Normal,     // Basic progress
    Debug       // Full details
};

class RsTestConfig {
public:
    // Parse command line arguments and create config
    static RsTestConfig parse_args(int argc, char* argv[]);

    // Access methods
    TestMode getMode() const { return mode; }
    int getTT() const { return tt; }
    int getNumCodewords() const { return num_codewords; }
    Verbosity getVerboseLevel() const { return verbose_level; }
    unsigned int getRandomSeed() const { return random_seed; }
    double getErrorRate() const { return error_rate > 0.0 ? error_rate : 0.01; }
    bool getVerifyCorrection() const { return verify_correction; }

    // Print configuration summary
    void print() const;

private:
    RsTestConfig() :
        mode(TestMode::Unknown),
        tt(0),
        num_codewords(10000),
        verbose_level(Verbosity::Normal),
        random_seed(1093),
        error_rate(0.0),
        verify_correction(true) {}

    TestMode mode;
    int tt;                     // Error correction capacity
    int num_codewords;          // Default 10000
    Verbosity verbose_level;
    unsigned int random_seed;   // For reproducibility
    double error_rate;          // Override default tt/(8*nn_short)
    bool verify_correction;     // Extra validation

    static void print_usage();
};

// Custom exception for configuration errors
class ConfigurationError : public std::runtime_error {
public:
    explicit ConfigurationError(const std::string& msg) : std::runtime_error(msg) {}
};

#endif // RS_TEST_CONFIG_H
