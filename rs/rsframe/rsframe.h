// Copyright (C) 1996, 1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

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
    ULONG cbCWPayloadPartial;    // Payload bytes in partial codeword.
    ULONG cbCWCheck;            // Check bytes in each codeword.
    ULONG cbCWMaxCorrections;    // Maximum number of correctable errors
                                // per codeword.
};

class RSFRAME
{
public:
    RSFRAME(void);
    ~RSFRAME(void);
    BOOL Design1(ULONG cbFrameLength, ULONG nMaxCorrectionsPerCodeword,
        struct RSFRAME_DESIGN_INFO *pDesignInfo);
    BOOL Encode(UBYTE *pPayload, UBYTE *pCheck);
    BOOL Decode(UBYTE *pPayload, UBYTE *pCheck, ULONG *iError);

private:
    BOOL bDesigned;
    struct RSFRAME_DESIGN_INFO DesignInfo;
    std::unique_ptr<RS_CODEC>  codec;
    GF work[nn];
};
