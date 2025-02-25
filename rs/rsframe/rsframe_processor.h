// Copyright (C) 1996, 1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
// Copyright (C) 2025 Stephen Dennis. All rights reserved.
//

#ifndef RSFRAME_PROCESSOR_H
#define RSFRAME_PROCESSOR_H

#include <memory>
#include "projdefs.h"
#include "crc.h"

// Maintain exact memory layout compatibility with original implementation
class RSFrameProcessor {
public:
    RSFrameProcessor(ULONG frameLength = 64 * 1024 - sizeof(UINT32), ULONG maxCorrections = 16)
        : frameLength_(frameLength), maxCorrections_(maxCorrections) {

        // Design the Reed-Solomon code
        if (!rsframe_.Design1(frameLength_, maxCorrections_, &designInfo_) ||
            designInfo_.cbPayloadLength <= sizeof(UINT64) + sizeof(UINT32) + sizeof(USHORT)) {
            throw std::runtime_error("Infeasible code design");
        }

        fillBytes_ = frameLength_ - designInfo_.cbFrameLength;

        // Allocate memory for payload and check data - exactly as in the original
        pPayload_ = new UBYTE[designInfo_.cbPayloadLength];
        pCheck_ = new UBYTE[designInfo_.cbCheckLength];

        if (!pPayload_ || !pCheck_) {
            if (pPayload_) delete[] pPayload_;
            if (pCheck_) delete[] pCheck_;
            throw std::runtime_error("Out of memory");
        }

        // Set up pointers exactly as in the original implementation
        pBlockCount_ = (INT64 *)pPayload_;
        pLength_ = (USHORT *)(pPayload_ + sizeof(INT64));
        pBuffer_ = pPayload_ + sizeof(INT64) + sizeof(USHORT);
        pCRC_ = (UINT32 *)(pPayload_ + designInfo_.cbPayloadLength - sizeof(UINT32));

        bufferSize_ = designInfo_.cbPayloadLength - sizeof(INT64) - sizeof(UINT32) - sizeof(USHORT);
        crcLength_ = designInfo_.cbPayloadLength - sizeof(UINT32);
    }

    ~RSFrameProcessor() {
        if (pPayload_) delete[] pPayload_;
        if (pCheck_) delete[] pCheck_;
    }

    // Disable copy and move to prevent double-free
    RSFrameProcessor(const RSFrameProcessor&) = delete;
    RSFrameProcessor& operator=(const RSFrameProcessor&) = delete;
    RSFrameProcessor(RSFrameProcessor&&) = delete;
    RSFrameProcessor& operator=(RSFrameProcessor&&) = delete;

    bool encode(const std::string& inputFile, const std::string& outputFile) {
        FILE* fpIn = fopen(inputFile.c_str(), "rb");
        if (!fpIn) {
            std::cerr << "Error: Cannot open input file: " << inputFile << std::endl;
            return false;
        }

        FILE* fpOut = fopen(outputFile.c_str(), "wb");
        if (!fpOut) {
            fclose(fpIn);
            std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
            return false;
        }

        // Set up file buffers exactly as in the original
        char inFileBuffer[16 * 1024];
        char outFileBuffer[16 * 1024];
        setvbuf(fpIn, inFileBuffer, _IOFBF, sizeof(inFileBuffer));
        setvbuf(fpOut, outFileBuffer, _IOFBF, sizeof(outFileBuffer));

        // Initialize exactly as in the original
        INT64 blockCount = -1;
        memset(pPayload_, 0, designInfo_.cbPayloadLength);

        size_t bytesRead;
        while ((bytesRead = fread(pBuffer_, 1, bufferSize_, fpIn)) > 0) {
            // Process block exactly as in the original
            *pBlockCount_ = blockCount++;
            *pLength_ = (USHORT)bytesRead;
            *pCRC_ = CRC32_ProcessBuffer(0UL, pPayload_, crcLength_);

            if (!rsframe_.Encode(pPayload_, pCheck_)) {
                std::cerr << "Fatal: Error encoding." << std::endl;
                fclose(fpOut);
                fclose(fpIn);
                return false;
            }

            fwrite("SVDL", 4, 1, fpOut);
            fwrite(pPayload_, designInfo_.cbPayloadLength, 1, fpOut);
            fwrite(pCheck_, designInfo_.cbCheckLength, 1, fpOut);

            memset(pPayload_, 0, designInfo_.cbPayloadLength);
            if (fillBytes_) {
                fwrite(pPayload_, fillBytes_, 1, fpOut);
            }
        }

        fclose(fpOut);
        fclose(fpIn);
        return true;
    }

    bool decode(const std::string& inputFile, const std::string& outputFile) {
        FILE* fpIn = fopen(inputFile.c_str(), "rb");
        if (!fpIn) {
            std::cerr << "Error: Cannot open input file: " << inputFile << std::endl;
            return false;
        }

        FILE* fpOut = fopen(outputFile.c_str(), "wb");
        if (!fpOut) {
            fclose(fpIn);
            std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
            return false;
        }

        // Set up file buffers exactly as in the original
        char inFileBuffer[16 * 1024];
        char outFileBuffer[16 * 1024];
        setvbuf(fpIn, inFileBuffer, _IOFBF, sizeof(inFileBuffer));
        setvbuf(fpOut, outFileBuffer, _IOFBF, sizeof(outFileBuffer));

        // Initialize exactly as in the original
        INT64 blockCount = -1;

        for (;;) {
            memset(pPayload_, 0, designInfo_.cbPayloadLength);

            // Read frame header
            char header[4];
            size_t headerRead = fread(header, 1, 4, fpIn);
            if (headerRead != 4) {
                break;  // End of file or error
            }

            // Check frame signature - exactly as in the original
            if (memcmp(header, "SVDL", 4) != 0) {
                std::cerr << "Warning: Frame flags don't match." << std::endl;
            }

            // Read payload and check blocks - exactly as in the original
            fread(pPayload_, designInfo_.cbPayloadLength, 1, fpIn);
            fread(pCheck_, designInfo_.cbCheckLength, 1, fpIn);

            // Decode block - exactly as in the original
            ULONG iError = 0;
            if (!rsframe_.Decode(pPayload_, pCheck_, &iError)) {
                std::cerr << "Error: Hard decoding error." << std::endl;
            }

            // Check block sequence - exactly as in the original
            if (*pBlockCount_ == -1) {
                blockCount = *pBlockCount_;
            } else if (blockCount != *pBlockCount_) {
                std::cerr << "Warning: blocks out of sequence." << std::endl;
            }
            blockCount++;

            // Check CRC - exactly as in the original
            if (CRC32_ProcessBuffer(0UL, pPayload_, crcLength_ + sizeof(UINT32)) != 0UL) {
                std::cerr << "Fatal: CRC-32 error." << std::endl;
                fclose(fpOut);
                fclose(fpIn);
                return false;
            }

            // Write decoded data
            size_t bytesToWrite = *pLength_;
            fwrite(pBuffer_, bytesToWrite, 1, fpOut);

            // Skip fill bytes if present - exactly as in the original
            if (fillBytes_) {
                fread(pPayload_, fillBytes_, 1, fpIn);
            }
        }

        fclose(fpOut);
        fclose(fpIn);
        return true;
    }

    // Get design information for debugging
    const RSFRAME_DESIGN_INFO& getDesignInfo() const {
        return designInfo_;
    }

private:
    ULONG frameLength_;
    ULONG maxCorrections_;
    ULONG fillBytes_;
    RSFRAME rsframe_;
    RSFRAME_DESIGN_INFO designInfo_;

    // Raw pointers to match original implementation exactly
    UBYTE* pPayload_;
    UBYTE* pCheck_;

    // Pointers into the payload buffer, exactly as in the original
    INT64* pBlockCount_;
    USHORT* pLength_;
    UBYTE* pBuffer_;
    UINT32* pCRC_;

    int bufferSize_;
    int crcLength_;
};

#endif // RSFRAME_PROCESSOR_H
