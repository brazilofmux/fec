// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
// Copyright (C) 2025 Stephen Dennis. All rights reserved.
//

#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <iostream>
#include <fcntl.h>

// Include necessary headers
#include "projdefs.h"
#include "rs_factory.h"
#include "rsframe.h"
#include "rsframe_processor.h"

// Helper function for opening files
bool rs_fopen(FILE** pFile, const char* filename, const char* mode) {
    if (pFile) {
        *pFile = nullptr;
        if (nullptr != filename && nullptr != mode) {
#if !defined(__INTEL_COMPILER) && (_MSC_VER >= 1400)
            // 1400 is Visual C++ 2005
            return (fopen_s(pFile, filename, mode) == 0);
#else
            *pFile = fopen((char*)filename, (char*)mode);
            if (nullptr != *pFile) {
                return true;
            }
#endif // WINDOWS_FILES
        }
    }
    return false;
}

// Command-line utility function
void print_usage(const char* programName) {
    std::cerr << "Usage: " << programName << " [-e|-d] infile outfile" << std::endl;
    std::cerr << "  -e            Encode mode" << std::endl;
    std::cerr << "  -d            Decode mode" << std::endl;
}

// ------------------- Main Function -------------------

int main(int argc, char* argv[]) {
    try {
        if (argc != 4) {
            print_usage(argv[0]);
            return 1;
        }

        std::string mode = argv[1];
        std::string inputFile = argv[2];
        std::string outputFile = argv[3];

        if (mode != "-e" && mode != "-d") {
            print_usage(argv[0]);
            return 1;
        }

        RSFrameProcessor processor;

        if (mode == "-e") {
            if (!processor.encode(inputFile, outputFile)) {
                return 1;
            }
            std::cout << "File encoded successfully." << std::endl;
        } else {
            if (!processor.decode(inputFile, outputFile)) {
                return 1;
            }
            std::cout << "File decoded successfully." << std::endl;
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
