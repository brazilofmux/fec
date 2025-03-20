// Copyright (C) 1996, 1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
// Copyright (C) 2025 Stephen Dennis. All rights reserved.
//

#ifndef RSFRAME_H
#define RSFRAME_H

#include <memory>
#include "projdefs.h"

struct RSFRAME_DESIGN_INFO
{
    // cbFrameLength = cbPayloadLength + cbCheckLength
    // These length dictate the size of buffers passed into
    // the encoding and decoding routines.
    //
    ULONG cbFrameLength;
    ULONG cbPayloadLength;
    ULONG cbCheckLength;

    // Codeword array design.
    //
    ULONG cCWFull;                // Must be non-zero.
    ULONG cbCWPayload;            // Payload bytes in each codeword.
    ULONG cbCWPayloadPartial;     // Payload bytes in partial codeword.
    ULONG cbCWCheck;              // Check bytes in each codeword.
    ULONG cbCWMaxCorrections;     // Maximum number of correctable errors
                                  // per codeword.
};

class RSFRAME
{
public:
    RSFRAME(void)
    {
        memset(work, 0, sizeof(work));
        bDesigned = FALSE;
        memset(&DesignInfo, 0, sizeof(DesignInfo));
    }

    ~RSFRAME(void)
    {
    }

    BOOL Design1(ULONG req_cbFrameLength, ULONG req_nMaxCorrectionsPerCodeword,
                 RSFRAME_DESIGN_INFO *pDesignInfo)
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
        codec = RS_FACTORY::instance().create_codec(t, (nn-2*t+1)/2);
        return TRUE;
    }

    BOOL Encode(UBYTE *pPayload, UBYTE *pCheck)
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

    BOOL Decode(UBYTE *pPayload, UBYTE *pCheck, ULONG *iError)
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

private:
    BOOL bDesigned;
    RSFRAME_DESIGN_INFO DesignInfo;
    std::unique_ptr<RS_CODEC> codec;
    GF work[nn];
};

#endif // RSFRAME_H
