// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <fcntl.h>

#include "projdefs.h"
#include "rs_factory.h"
#include "rsframe.h"
#include "crc.h"

bool rs_fopen(FILE** pFile, const char* filename, const char* mode)
{
    if (pFile)
    {
        *pFile = nullptr;
        if (  nullptr != filename
           && nullptr != mode)
        {
#if !defined(__INTEL_COMPILER) && (_MSC_VER >= 1400)
            // 1400 is Visual C++ 2005
            //
            return (fopen_s(pFile, filename, mode) == 0);
#else
            *pFile = fopen((char*)filename, (char*)mode);
            if (nullptr != *pFile)
            {
                return true;
            }
#endif // WINDOWS_FILES
        }
    }
    return false;
}

RSFRAME::RSFRAME(void)
{
    memset(work, 0, sizeof(work));
    bDesigned = FALSE;
    memset(&DesignInfo, 0, sizeof(DesignInfo));
}

BOOL RSFRAME::Design1
(
    ULONG req_cbFrameLength,
    ULONG req_nMaxCorrectionsPerCodeword,
    struct RSFRAME_DESIGN_INFO *pDesignInfo
)
{
    ULONG T = req_cbFrameLength;
    ULONG t = req_nMaxCorrectionsPerCodeword;

TryAgain:
    if (T == 0) return FALSE;
    if (t == 0 || MAX_TT < t) return FALSE;

    ULONG CW_P = 2*t;  // Check bytes per codeword.

    // Calculate the minimum number of codewords we need.
    //
    ULONG N = (T+255-1)/255;

TryAgain2:
    // Calculate the minimum check length
    //
    ULONG P = CW_P*N; // Total number of check bytes.

    // Calculate the maximum payload length.
    //
    if (T <= P)
    {
        N--;
        P = CW_P*N;
        if (T <= P) return FALSE;
    }
    ULONG D = T-P;   // Total number of payload bytes.

    // Calculate maximum codeword length
    //
    ULONG CW = T/N;

    // Calculate maximum payload per codeword.
    //
    if (CW <= CW_P) return FALSE;
    ULONG CW_D = CW - CW_P;

    // Calculate length of partial codeword.
    //
    if (T < N*CW) return FALSE;
    ULONG CWPartial = T - N*CW;
    if (CWPartial >= CW)
    {
        N++;
        goto TryAgain2;
    }

    if (0 < CWPartial && CWPartial <= CW_P)
    {
        // The partial codeword isn't big enough to cover it's own checking
        // so, we need to reduce the size of the frame to get rid of it.
        //
        T -= CWPartial;
        goto TryAgain;
    }

    // Calculate payload of partial keyword.
    //
    ULONG CW_D_Partial = 0;
    if (CWPartial != 0)
    {
        CW_D_Partial = CWPartial - CW_P;

        // Re-Calculate the check length including partial codeword.
        //
        P = CW_P*(N+1); // Total number of check bytes.

        // Re-Calculate the maximum payload length.
        //
        if (T <= P) return FALSE;
        D = T-P;   // Total number of payload bytes.
    }

    DesignInfo.cbFrameLength = T;
    DesignInfo.cbPayloadLength = D;
    DesignInfo.cbCheckLength = P;

    DesignInfo.cCWFull = D/CW_D;
    DesignInfo.cbCWPayload = CW_D;
    DesignInfo.cbCWPayloadPartial = CW_D_Partial;
    DesignInfo.cbCWCheck = 2*t;
    DesignInfo.cbCWMaxCorrections = t;

    bDesigned = TRUE;
    *pDesignInfo = DesignInfo;
    codec = RS_FACTORY::instance().create_specialized_codec(t, (nn-2*t+1)/2);
    return TRUE;
}

RSFRAME::~RSFRAME(void)
{
}

BOOL RSFRAME::Encode(UBYTE *pPayload, UBYTE *pCheck)
{
    if (!bDesigned) return FALSE;
    int mm = DesignInfo.cbCWCheck;
    int kk = DesignInfo.cbCWPayload;
    int kk_Partial = DesignInfo.cbCWPayloadPartial;
    ULONG i;
    for (i = 0; i < DesignInfo.cCWFull; i++)
    {
        memset(work, 0, MAX_KK);
        memcpy(work, pPayload + i*kk, kk);
        codec->RSEncode(work, pCheck + i*mm);
    }
    if (DesignInfo.cbCWPayloadPartial)
    {
        memset(work, 0, MAX_KK);
        memcpy(work, pPayload + i*kk, kk_Partial);
        codec->RSEncode(work, pCheck + i*mm);
    }
    return TRUE;
}

BOOL RSFRAME::Decode(UBYTE *pPayload, UBYTE *pCheck, ULONG *iError)
{
    if (!bDesigned) return FALSE;
    int mm = DesignInfo.cbCWCheck;
    int kk = DesignInfo.cbCWPayload;
    int kk_Partial = DesignInfo.cbCWPayloadPartial;
    ULONG i;
    for (i = 0; i < DesignInfo.cCWFull; i++)
    {
        memset(work, 0, nn);
        memcpy(work, pCheck + i*mm, mm);
        memcpy(work + mm, pPayload + i*kk, kk);
        int cc = codec->RSDecode(work);
        if (cc < 0 || (ULONG)cc == DesignInfo.cbCWMaxCorrections)
        {
            *iError = i*kk;
            return FALSE;
        }
        memcpy(pPayload + i*kk, work + mm, kk);
        memcpy(pCheck + i*mm, work, mm);
    }
    if (DesignInfo.cbCWPayloadPartial)
    {
        memset(work, 0, nn);
        memcpy(work, pCheck + i*mm, mm);
        memcpy(work + mm, pPayload + i*kk, kk_Partial);
        int cc= codec->RSDecode(work);
        if (cc < 0 || (ULONG)cc == DesignInfo.cbCWMaxCorrections)
        {
            *iError = i*kk;
            return FALSE;
        }
        memcpy(pPayload + i*kk, work + mm, kk_Partial);
        memcpy(pCheck + i*mm, work, mm);
    }
    return TRUE;
}

void usage(char *ProgramName)
{
    fprintf(stderr, "Usage: %s [-e|-d] infile outfile\r\n", ProgramName);
    exit(0);
}

int Mode = 0;
#define ENCODING 1
#define DECODING 2

char out_file_buffer[16*1024];
char InFileBuffer[16*1024];

int main(int argc, char *argv[])
{
    if (argc != 4) usage(argv[0]);
    if (argv[1][0] != '-') usage(argv[0]);
    if (argv[1][2] != '\0') usage(argv[0]);
    if (argv[1][1] == 'e' || argv[1][1] == 'E')
    {
        Mode = ENCODING;
    } else if (argv[1][1] == 'd' || argv[1][1] == 'D')
    {
        Mode = DECODING;
    }
    else usage(argv[0]);

    RS_Init();
    RSFRAME rsframe;
    RSFRAME_DESIGN_INFO InfoStruct;
#define FRAME_LENGTH (64*1024-sizeof(UINT32))
    if (  !rsframe.Design1(FRAME_LENGTH, 16UL, &InfoStruct)
       || InfoStruct.cbPayloadLength <= sizeof(UINT64) + sizeof(UINT32) + sizeof(USHORT))
    {
        fprintf(stderr, "Fatal: Infeasible code design.\n");
        exit(0);
    }
    int FillBytes = FRAME_LENGTH - InfoStruct.cbFrameLength;
    UBYTE *pPayload = new UBYTE[InfoStruct.cbPayloadLength];
    UBYTE *pCheck   = new UBYTE[InfoStruct.cbCheckLength];
    if (!pPayload || !pCheck)
    {
        fprintf(stderr, "Fatal: Out of memory.\r\n");
        if (pPayload) delete[] pPayload;
        if (pCheck) delete[] pCheck;
        return 0;
    }

    FILE* fpIn = nullptr;
    if (!rs_fopen(&fpIn, argv[2], "rb") || nullptr == fpIn)
    {
        delete[] pCheck;
        delete[] pPayload;
        return 0;
    }
    setvbuf(fpIn, InFileBuffer, _IOFBF, sizeof(InFileBuffer));
    FILE* fpOut = nullptr;
    if (!rs_fopen(&fpOut, argv[3], "wb") || nullptr == fpOut)
    {
        fclose(fpIn);
        delete[] pCheck;
        delete[] pPayload;
        return 0;
    }
    setvbuf(fpOut, out_file_buffer, _IOFBF, sizeof(out_file_buffer));

    INT64 *pBlockCount     = (INT64 *)pPayload;
    USHORT  *pLength       = (USHORT *)(pPayload + sizeof(INT64));
    UBYTE   *pBuffer       = pPayload + sizeof(INT64) + sizeof(USHORT);
    UINT32  *pCRC          = (UINT32 *)(pPayload + InfoStruct.cbPayloadLength - sizeof(UINT32));
    int      nBuffer       = InfoStruct.cbPayloadLength - sizeof(INT64) - sizeof(UINT32) - sizeof(USHORT);
    int      nCRCLength    = InfoStruct.cbPayloadLength - sizeof(UINT32);

    if (ENCODING == Mode)
    {
        // Initialize
        //
        INT64 BlockCount = -1;
        memset(pPayload, 0, InfoStruct.cbPayloadLength);
        size_t cc = fread(pBuffer, 1, nBuffer, fpIn);
        while (cc > 0)
        {
            // Do a block
            //
            *pBlockCount = BlockCount++;
            *pLength = (USHORT)cc;
            *pCRC = CRC32_ProcessBuffer(0UL, pPayload, nCRCLength);
            if (!rsframe.Encode(pPayload, pCheck))
            {
                fprintf(stderr, "Fatal: Error encoding.\r\n");
                break;
            }
            fwrite("SVDL", 4, 1, fpOut);
            fwrite(pPayload, InfoStruct.cbPayloadLength, 1, fpOut);
            fwrite(pCheck, InfoStruct.cbCheckLength, 1, fpOut);
            memset(pPayload, 0, InfoStruct.cbPayloadLength);
            if (FillBytes) fwrite(pPayload, FillBytes, 1, fpOut);
            cc = fread(pBuffer, 1, nBuffer, fpIn);
        }
    }
    else
    {
        // Initialize
        //
        INT64 BlockCount = -1;

        for (;;)
        {
            memset(pPayload, 0, InfoStruct.cbPayloadLength);
            size_t cc = fread(pBuffer, 4, 1, fpIn);
            if (cc != 1)
            {
                break;
            }
            if (memcmp(pBuffer, "SVDL", 4))
            {
                fprintf(stderr, "Warning: Frame flags don't match.\r\n");
            }
            fread(pPayload, InfoStruct.cbPayloadLength, 1, fpIn);
            fread(pCheck, InfoStruct.cbCheckLength, 1, fpIn);

            // Do a block
            //
            ULONG iError = 0;
            if (!rsframe.Decode(pPayload, pCheck, &iError))
            {
                fprintf(stderr, "Error: Hard decoding error.\r\n");
            }

            if (*pBlockCount == -1)
            {
                BlockCount = *pBlockCount;
            }
            else if (BlockCount != *pBlockCount)
            {
                fprintf(stderr, "Warning: blocks out of sequence.\r\n");
            }
            BlockCount++;

            cc = *pLength;
            if (CRC32_ProcessBuffer(0UL, pPayload, nCRCLength+sizeof(UINT32)) != 0UL)
            {
                fprintf(stderr, "Fatal: CRC-32 error.\r\n");
                break;
            }
            fwrite(pBuffer, cc, 1, fpOut);
            if (FillBytes) fread(pPayload, FillBytes, 1, fpIn);
        }
    }
    fclose(fpOut);
    fclose(fpIn);
    delete[] pPayload;
    delete[] pCheck;
}
