#include <iostream>
#include <iomanip>
#include <cstring>
#include <random>
#include <vector>
#include <functional>
#include "rs.h"          // Original implementation
#include "rsref.h"       // Reference implementation 
#include "rskarn_wrapper.h" // Wrapper for Karn's implementation
#include "RsVerification.h"

void print_hex(const char* label, const GF* data, int len) {
    std::cout << label << ": ";
    for (int i = 0; i < len; i++) {
        std::cout << std::hex << std::setfill('0') << std::setw(2) 
                  << static_cast<int>(data[i]) << " ";
        if ((i + 1) % 16 == 0 && i + 1 < len) {
            std::cout << "\n" << std::string(strlen(label) + 2, ' ');
        }
    }
    std::cout << std::dec << std::endl;
}

struct EncoderResults {
    GF bb[2*MAX_TT];
    void clear() { memset(bb, 0, sizeof(bb)); }
};

void compare_results(const EncoderResults& karn, 
                    const EncoderResults& ref,
                    const EncoderResults& orig,
                    int tt) {
    bool all_match = true;
    std::cout << "\nComparing parity bytes:\n";
    print_hex("Karn     ", karn.bb, 2*tt);
    print_hex("Reference", ref.bb, 2*tt);
    print_hex("Original ", orig.bb, 2*tt);
    
    for (int i = 0; i < 2*tt; i++) {
        if (karn.bb[i] != ref.bb[i] || karn.bb[i] != orig.bb[i]) {
            all_match = false;
            std::cout << "Mismatch at position " << i << ":\n"
                      << "  Karn: 0x" << std::hex << static_cast<int>(karn.bb[i])
                      << "  Ref: 0x" << static_cast<int>(ref.bb[i])
                      << "  Orig: 0x" << static_cast<int>(orig.bb[i])
                      << std::dec << std::endl;
        }
    }
    
    if (all_match) {
        std::cout << "✓ All implementations match\n";
    } else {
        std::cout << "✗ Implementations differ\n";
    }
}

void run_standard_tests(int tt) {
    RS_ENCODER_KARN karn(tt);
    RS_ENCODER_REF ref(tt);
    RS_ENCODER orig(tt);

    const int kk = nn - 2*tt;  // Data length
    
    // Define test cases
    struct TestCase {
        const char* name;
        void (*fill_data)(GF*);  // Function pointer instead of std::function
    };
    
    // Create test cases array
    TestCase test_cases[] = {
        {"Zero data", +[](GF* data) { memset(data, 0, MAX_KK); }},
        {"All ones", +[](GF* data) { memset(data, 1, MAX_KK); }},
        {"Sequential", +[](GF* data) { 
            for (int i = 0; i < MAX_KK; i++) 
                data[i] = i & 0xFF; 
        }},
        {"Random data", +[](GF* data) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(0, 255);
            for (int i = 0; i < MAX_KK; i++) {
                data[i] = dis(gen);
            }
        }}
    };

    for (const auto& test : test_cases) {
        std::cout << "\nRunning test: " << test.name << std::endl;
        
        GF data[MAX_KK];
        test.fill_data(data);
        print_hex("Input", data, kk);

        EncoderResults karn_result, ref_result, orig_result;
        karn_result.clear();
        ref_result.clear();
        orig_result.clear();

        karn.RSEncode(data, karn_result.bb);
        ref.RSEncode(data, ref_result.bb);
        orig.RSEncode(data, orig_result.bb);

        compare_results(karn_result, ref_result, orig_result, tt);
    }
}

int main(int argc, char* argv[]) {
    // Initialize all implementations
    RS_Init();  // Common initialization
    RS_ENCODER::Init();
    
    std::cout << "Reed-Solomon Encoder Comparison Test\n";
    std::cout << "===================================\n";

    // Test different error correction capacities
    int tt_values[] = {1, 2, 4, 8, 16, 32};
    
    for (int tt : tt_values) {
        std::cout << "\nTesting with tt=" << tt << " (correctable errors)\n";
        std::cout << "--------------------------------------------\n";
        run_standard_tests(tt);
    }

    return 0;
}
