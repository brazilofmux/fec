// Copyright (C) 1996, 1997 Solid Vertical Domains, Ltd. All rights reserved.
// Copyright (C) 2008 Stephen Dennis. All rights reserved.
//

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include "projdefs.h"
#include "rs.h"
#include "crc.h"

long int nrand48(unsigned short seed[3])
{
    return rand();
}

double erand48(unsigned short seed[3])
{
    return ((double)rand())/0xFFFFUL;
}

// Each segment is made up of 16 consecutively numbered sections. The first
// section number is a multiple of 16 (i.e., 16n, 16n+1, 16n+2, ..., 16n+15,
// where n is the segment number).
//
// Each section is 4KB long and consists of a 4 byte flag ("SVDL"), a 4 byte
// section number, 3468 bytes of C2 data+parity, 4 byte CRC-32, and 616 bytes
// of C1 parity bytes. C1 Data consists of the section number, C2 data+parity,
// and the 4 byte CRC-32 in an interlaced format to be defined later.
//
// C1 is a Reed-Solomon (186,158) 14-error correcting code. There are 22
// C1 codewords per section.
//
// C2 is a Reed-Solomon (192,164) 14-error correcting code. There are 289
// C2 codewords per segment. Each codeword occupies 12 bytes of each
// section -- below the 14-error correcting capability of it's companion C1
// code.
//
// The C2 data consists of 47392 bytes of user data and a 4 byte CRC-32 to
// guard against miscorrections. 
//
// The overall space utilization is 72.31%.
//

#define SIZEOF_C1_CODE       186
#define SIZEOF_C1_DATA       158
#define SIZEOF_C1_PARITY   (SIZEOF_C1_CODE-SIZEOF_C1_DATA)
#define NUM_C1_CODEWORDS      22
#define SIZEOF_C2_CODE       192
#define SIZEOF_C2_DATA       164
#define SIZEOF_C2_PARITY   (SIZEOF_C2_CODE-SIZEOF_C2_DATA)
#define NUM_C2_CODEWORDS     289
#define SIZEOF_CLIENT_DATA 47392
#define NUM_SECTIONS          16
#define SIZEOF_FLAG            4
#define SIZEOF_NUMBER          4
#define SIZEOF_SEC_DATA        (SIZEOF_C2_CODE*NUM_C2_CODEWORDS/NUM_SECTIONS)
#define SIZEOF_CRC32           4
#define SIZEOF_SEC_PARITY      (SIZEOF_C1_PARITY*NUM_C1_CODEWORDS)
#define SIZEOF_SECTION         4096
#define SIZEOF_SECTION2        (SIZEOF_FLAG+SIZEOF_NUMBER+SIZEOF_SEC_DATA+SIZEOF_CRC32+SIZEOF_SEC_PARITY)
#if (SIZEOF_SECTION == SIZEOF_SECTION2)
#else
#error Invalid Size Assumptions
#endif

typedef struct
{
    UBYTE Flag[4];
    ULONG Number;
    UBYTE C2[SIZEOF_C2_CODE*NUM_C2_CODEWORDS/NUM_SECTIONS];
    ULONG CRC32;
    UBYTE C1Parity[SIZEOF_C1_PARITY*NUM_C1_CODEWORDS];
} SECTION, *PSECTION;

SECTION Section[NUM_SECTIONS];

struct
{
    UBYTE C2Parity[SIZEOF_C2_PARITY];
    UBYTE C2Data[SIZEOF_C2_DATA];
} C2Codewords[NUM_C2_CODEWORDS];


// Copy the C1 Data and Parity bytes of the specified C1 codeword from
// the specified section.
//
// The first data byte is located at SIZEOF_C1_DATA*iCodeword+SIZEOF_FLAG from
// the beginning of the section. The first parity byte is located at
// SIZEOF_C1_PARITY*iCodeword within C1Parity.
//
// If pC1Data is NULL, then data bytes are not copied. Likewise, if pC1Parity
// is NULL, then parity bytes are not copied.
//
// Codewords are numbered starting at 0.
//
void GetC1Codeword
(   PSECTION pSection,
    int      iCodeword,
    UBYTE *  pC1Data,
    UBYTE *  pC1Parity)
{
    if (pC1Data != NULL)
    {
        UBYTE *p = (UBYTE *)(&(pSection->Number));
        p += SIZEOF_C1_DATA*iCodeword;
        memcpy(pC1Data, p, SIZEOF_C1_DATA);
    }
    if (pC1Parity != NULL)
    {
        UBYTE *p = pSection->C1Parity;
        p += SIZEOF_C1_PARITY*iCodeword;
        memcpy(pC1Data, p, SIZEOF_C1_PARITY);
    }
}

// Copy the C1 Data and Parity bytes of the specified C1 codeword to
// the specified section.
//
// The first data byte is located at SIZEOF_C1_DATA*iCodeword+SIZEOF_FLAG from
// the beginning of the section. The first parity byte is located at
// SIZEOF_C1_PARITY*iCodeword within C1Parity.
//
// If pC1Data is NULL, then data bytes are not copied. Likewise, if pC1Parity
// is NULL, then parity bytes are not copied.
//
// Codewords are numbered starting at 0.
//
void SetC1Codeword
(   PSECTION pSection,
    int      iCodeword,
    UBYTE *  pC1Data,
    UBYTE *  pC1Parity)
{
    if (pC1Data != NULL)
    {
        UBYTE *p = (UBYTE *)(&(pSection->Number));
        p += SIZEOF_C1_DATA*iCodeword;
        memcpy(p, pC1Data, SIZEOF_C1_DATA);
    }
    if (pC1Parity != NULL)
    {
        UBYTE *p = pSection->C1Parity;
        p += SIZEOF_C1_PARITY*iCodeword;
        memcpy(p, pC1Data, SIZEOF_C1_PARITY);
    }
}

// Set Client Data into Segment in the right places.
//
void SetClientDataAndCRC32(UBYTE *pClientData, ULONG *pCRC32)
{
    int cData = SIZEOF_CLIENT_DATA;
    int cC2Parity = SIZEOF_C2_PARITY;
    int cC2Data = 0;
    UBYTE *p = pClientData;
    for (int i = 0; i < SIZEOF_C2_DATA; i++)
    {
        for (int j = 0; j < NUM_SECTIONS; j++)
        {
            // Note: The (cC2Parity == 0) AND (cC2Data == 0) state is illegal.
            // cC2Parity != 0 is the "Skip Parity Symbols" State, and
            // cC2Parity != 0 is the "Do Data Symbols" State.
            // The (cC2Parity != 0) AND (cC2Data != 0) state is illegal.
            //
            if (cC2Parity)
            {
                // We are in the "Skip Parity Symbols" State
                if (--cC2Parity == 0)
                {
                    cC2Data = SIZEOF_C2_DATA;
                }
            }
            else
            {
                // We are in the "Do Data Symbols" State.
                if (!cData)
                {
                    p = (UBYTE *)(pCRC32);
                    cData = SIZEOF_CRC32;
                }
                Section[j].C2[i] = *(p++);
                if (--cC2Data == 0)
                {
                    cC2Parity = SIZEOF_C2_PARITY;
                }
            }
        }
    }
}

// Get Client Data from the right places in the Segment.
//
void GetClientDataAndCRC32(UBYTE *pClientData, ULONG *pCRC32)
{
    HEARTBEAT;
    int cData = SIZEOF_CLIENT_DATA;
    int cC2Parity = SIZEOF_C2_PARITY;
    int cC2Data = 0;
    UBYTE *p = pClientData;
    for (int i = 0; i < SIZEOF_C2_DATA; i++)
    {
        for (int j = 0; j < NUM_SECTIONS; j++)
        {
            // Note: The (cC2Parity == 0) AND (cC2Data == 0) state is illegal.
            // cC2Parity != 0 is the "Skip Parity Symbols" State, and
            // cC2Parity != 0 is the "Do Data Symbols" State.
            // The (cC2Parity != 0) AND (cC2Data != 0) state is illegal.
            //
            if (cC2Parity)
            {
                // We are in the "Skip Parity Symbols" State
                if (--cC2Parity == 0)
                {
                    cC2Data = SIZEOF_C2_DATA;
                }
            }
            else
            {
                // We are in the "Do Data Symbols" State.
                if (!cData)
                {
                    p = (UBYTE *)(pCRC32);
                    cData = SIZEOF_CRC32;
                }
                *(p++) = Section[j].C2[i];
                if (--cC2Data == 0)
                {
                    cC2Parity = SIZEOF_C2_PARITY;
                }
            }
        }
    }
}

// Calculate CRC32 on ClientData and organize the result into C2 codewords.
// Calculate C2 parity and place result into Segment. Add flags, counts, etc.
// Calculate C1 parity and finish constructing Segment.
//
void EncodeSegment(ULONG SegmentNumber, UBYTE ClientData[SIZEOF_CLIENT_DATA])
{
    HEARTBEAT;
    ULONG CRC32;
    CRC32 = lpUpdateCRC(0, SIZEOF_CLIENT_DATA, ClientData,
                        pulCRCTable);

    SetClientDataAndCRC32(ClientData, &CRC32);

    // Zero-out Reed-Solomon data
    //
    memset(data, 0, 255);

    // Transfer client data and it's CRC-32 to C2Codewords.
    //
    HEARTBEAT;
    UBYTE *p = ClientData;
    int cnt = SIZEOF_CLIENT_DATA;
    for (int j = 0; j < NUM_C2_CODEWORDS; j++)
    {
        HEARTBEAT;
        for (int i = 0; i < SIZEOF_C2_DATA; i++)
        {
            data[i] = *p;
            C2Codewords[j].C2Data[i] = *(p++);
            if (--cnt == 0)
            {
                p = ubCRC;
                cnt = sizeof ULONG;
            }
        }

        // Generate C2 Parity Symbols.
        //
        HEARTBEAT;
        encode_rs();
        for (i = 0; i < SIZEOF_C2_PARITY; i++)
        {
            C2Codewords[j].C2Parity[i] = bb[i];
        }
    }

    // Zero-out Reed-Solomon data
    //
    HEARTBEAT;
    memset(data, 0, 255);

    // Transfer C2Codewords into Segment.
    //
    HEARTBEAT;
    p = (UBYTE *)C2Codewords;
    ULONG SectionNumber = NUM_SECTIONS*SegmentNumber;
    for (j = 0; j < NUM_SECTIONS; j++)
    {
        HEARTBEAT;
        memcpy(Section[j].Flag, "SVDL", sizeof Section[j].Flag);
        Section[j].Number = SectionNumber++;
        memcpy(Section[j].C2, p, sizeof Section[j].C2);
        p += sizeof Section[j].C2;

        HEARTBEAT;
        UBYTE *pC1Data = (UBYTE *)(&Section[j].Number);
        Section[j].CRC32 = lpUpdateCRC(0,
                  NUM_C1_CODEWORDS*SIZEOF_C1_DATA-sizeof ULONG,
                  pC1Data, pulCRCTable);

        HEARTBEAT;
        UBYTE *pC1Parity = Section[j].C1Parity;
        for (int i = 0; i < NUM_C1_CODEWORDS; i++)
        {
            memcpy(data, pC1Data, SIZEOF_C1_DATA);
            pC1Data += SIZEOF_C1_DATA;
            encode_rs();
            memcpy(pC1Parity, bb, SIZEOF_C1_PARITY);
            pC1Parity += SIZEOF_C1_PARITY;
        }
    }
}

int DecodeSegment
(
    UBYTE UserData[SIZEOF_CLIENT_DATA],
    SECTION Section[NUM_SECTIONS],
    int Erasures[NUM_SECTIONS]
)
{
    // Check CRC-32, first.
    //
    for (int j = 0; j < NUM_SECTIONS; j++)
    {
        if (Erasures[j]) continue;
        if (memcmp(Section[j].Flag, "SVDL", sizeof(Section[j].Flag)) != 0)
        {
            Erasures[j] = TRUE;
            continue;
        }
        UBYTE *pC1Data = (UBYTE *)(&Section[j].Number);
        ULONG ulCRC32 = lpUpdateCRC(0, NUM_C1_CODEWORDS*SIZEOF_C1_DATA,
                  pC1Data, pulCRCTable);
        if (ulCRC32 != 0UL)
        {
            printf("CRC32(0x%08lx) failed for section %d.\n", ulCRC32, j);
            Erasures[j] = TRUE;
        }
    }

#if 0
    // Check the section number base number first. All sections should
    // share the same the same base number and be offsets from that base
    // number.
    //
    int BestBase = 0;
    int BestMatches = 0;
    for (int i = 0; i < NUM_SECTIONS; i++)
    {
        if (Erasures[i]) continue;
        int Matches = 1;
        int Base = Section[i].Number/NUM_SECTIONS;
        for (j = i+1; j < NUM_SECTIONS; j++)
        {
            if (Base == Section[j].Number/NUM_SECTIONS)
            {
                Matches++;
            }
        }
        if (Matches > BestMatches)
        {
            BestMatches = Matches;
            BestBase = Base;
        }
    }
    if (BestMatches <= 0)
    {
        printf("No clear base section number.\n");
        return FALSE;
    }

    // Mark all of the non-matches as errors
    //
    for (i = 0; i < NUM_SECTIONS; i++)
    {
        if (Section[i].Number/NUM_SECTIONS != BestBase)
        {
            printf("Mismatched base section number for section %d.\n", i);
            Erasures[i] = TRUE;
        }
    }

    // Verify that all the offsets are unique
    //
    int OffsetPresent[NUM_SECTIONS];
    for (i = 0; i < NUM_SECTIONS; i++)
    {
        OffsetPresent[i] = FALSE;
    }
    for (i = 0; i < NUM_SECTIONS; i++)
    {
        if (Erasures[i]) continue;
        int iOffset = Section[i].Number % NUM_SECTIONS;
        if (OffsetPresent[iOffset])
        {
            printf("Section %d is not unique.\n", i);
            Erasures[i] = TRUE;
        }
        else
        {
            OffsetPresent[iOffset] = TRUE;
        }
    }

    // Verify that the offset corresponds to the section position
    // Move the section into the correct position, if it's not already
    // there.
    //
    for (i = 0; i < NUM_SECTIONS; i++)
    {
        if (Erasures[i]) continue;
        int iOffset = Section[i].Number % NUM_SECTIONS;
        if (iOffset == i ) continue;

        // At this point, the section number isn't where it says that it
        // should be. Let's ask a few questions about the spot it says
        // that it should go.
        //
        if (Erasures[iOffset]) continue;
        if (iOffset != i)
        {
            printf("Swapping Sections %d and %d.\n", iOffset, i);
            SECTION tmp;
            tmp = Section[i];
            Section[i] = Section[iOffset];
            Section[iOffset] = tmp;
            int tmp1;
            tmp1 = Erasures[i];
            Erasures[i] = Erasures[iOffset];
            Erasures[iOffset] = tmp1;
        }
    }

    // Ok, now every section is correctly position and/or marked for erasure
    // Do, we have any erasures?
    //
    int TotalErasures = 0;
    for (i = 0; i < NUM_SECTIONS; i++)
    {
        if (Erasures[i]) TotalErasures++;
    }

    if (TotalErasures)
    {
        //printf("TotalErasures = %d\n", TotalErasures);
        // Employ RS code to correct errors
        //
        for (j = 0; j < 292; j++)
        {
            int erased_symbols[nn];
            int *pes = erased_symbols;
            int how_many_erased_symbols = 0;
            int m = 0;
            int Position = j;
            for (int i = 0; i < NUM_SECTIONS; i++)
            {
                for (int k = 0; k < 14; k++)
                {
                    if (Erasures[i] && how_many_erased_symbols < nn)
                    {
                        how_many_erased_symbols++;
                        *(pes++) = m;
                    }
                    recd[m++] = Section[i].Data[Position];
                    Position += 292;
                    if (Position >= SIZEOF_SECTION_DATA)
                        Position -= SIZEOF_SECTION_DATA;
                }
                Position += 18;
                if (Position >= SIZEOF_SECTION_DATA)
                    Position -= SIZEOF_SECTION_DATA;
            }

            // Pad with zeros.
            memset(recd+224, 0, 255-224);

            // Decode codeword in recd[]
            //
            if (how_many_erased_symbols > 28)
            {
                printf("Too many erasures %d.\n", how_many_erased_symbols);
                return FALSE;
            }
#ifndef NO_PRINT
            printf("Symbol erasures = %d.\n", how_many_erased_symbols);
            for (i = 0; i < how_many_erased_symbols; i++)
            {
                printf("%d ", erased_symbols[i]);
            }
#endif
            int cc;
            //cc = decode_rs();
            //if (cc != 0 && cc != 1)
                cc = eras_dec_rs(erased_symbols, how_many_erased_symbols);
            if (cc == 0 || cc == 1)
            {
                // Corrected successfully
                //
                UBYTE *p = recd;
                int Position = j;
                for (int i = 0; i < NUM_SECTIONS; i++)
                {
                    for (int k = 0; k < 14; k++)
                    {
                        Section[i].Data[Position] = *(p++);
                        Position += 292;
                        if (Position >= SIZEOF_SECTION_DATA)
                            Position -= SIZEOF_SECTION_DATA;
                    }
                    Position += 18;
                    if (Position >= SIZEOF_SECTION_DATA)
                        Position -= SIZEOF_SECTION_DATA;
                }
            }
            else
            {
                // Could not correct successfully
                //
                printf("Codeword (%d) decoder failure %d.\n", j, cc);
                return FALSE;
            }
        }

        // Force the section numbers to be correct.
        //
        for (i = 0; i < NUM_SECTIONS; i++)
        {
            Section[i].Number = BestBase*NUM_SECTIONS+i;
        }
    
        // Check CRC-32 again to make sure everything is wonderful.
        // BUG!: ECC doesn't protect the section numbers or the CRC-32, so
        // when errors occur, the CRC-32 may still not be correct.
        //
        for (int j = 0; j < NUM_SECTIONS; j++)
        {
            if (!DecodeCRC32(Section+j))
            {
                printf("Unexpected CRC32 failure for section %d.\n", j);
                //return FALSE;
            }
        }
    }

    // Now that all the errors that can be corrected are corrected, let's
    // transfer all the users data out of the segment.
    //
    UBYTE *pUserData = UserData;
    for (j = 0; j < 292; j++)
    {
        int Position = j+36;
        for (int i = 2; i < NUM_SECTIONS; i++)
        {
            for (int k = 0; k < 14; k++)
            {
                *(pUserData++) = Section[i].Data[Position];
                Position += 292;
                if (Position >= SIZEOF_SECTION_DATA)
                    Position -= SIZEOF_SECTION_DATA;
            }
            Position += 18;
            if (Position >= SIZEOF_SECTION_DATA)
                Position -= SIZEOF_SECTION_DATA;
        }
    }
#endif
    ULONG CRC32;
    GetClientDataAndCRC32(UserData, &CRC32);
    
  
    ULONG CRC32Received = lpUpdateCRC(0, SIZEOF_CLIENT_DATA, UserData,
                        pulCRCTable);
    if (CRC32Received != CRC32)
    {
        printf("Client Data failed CRC check.\n");
        exit(-1);
    }
    return TRUE;
}

#define MAX_WRITE_SIZE 16384
// Returns 1 if successfull, 0 if disk is full, and -1 is error.
//
BOOL WriteToFile(int hFile, UBYTE *pBuffer, USHORT cAsk, USHORT &cDone)
{
    cDone = 0;
    while (cAsk)
    {
        int cRequest = cAsk;
        if (cRequest > MAX_WRITE_SIZE)
        {
            cRequest = MAX_WRITE_SIZE;
        }
        int cc = _write(hFile, pBuffer, cRequest);
        if (cc > 0)
        {
            cAsk -= cc;
            pBuffer += cc;
            cDone += cc;
        }
        else return cc;
    }
    return 1;
}

#define MAX_READ_SIZE 16384
// Return 1 if successfull, 0 if EOF, and -1 if error.
//
int ReadFromFile(int hFile, UBYTE *pBuffer, USHORT cAsk, USHORT &cDone)
{
    cDone = 0;
    while (cAsk)
    {
        int cRequest = cAsk;
        if (cRequest > MAX_READ_SIZE)
        {
            cRequest = MAX_READ_SIZE;
        }
        int cc = _read(hFile, pBuffer, cRequest);
        if (cc > cRequest)
        {
            printf("Internal Read Error.\n");
            exit(-1);
        }
        else if (cc > 0)
        {
            cAsk -= cc;
            pBuffer += cc;
            cDone += cc;
        }
        else // EOF(cc==0) and ERROR(cc<0)
        {
            return cc;
        }
    }
    return 1;
}

class FileWithErrors
{
    int   hFile;
    LONG  Length;
    LONG  Offset;
public:
    void OpenForRead(char *filename);
    void OpenForWrite(char *filename);
    void Read(UBYTE *pBuffer, USHORT cAsk);
    void Write(UBYTE *pBuffer, USHORT cAsk);
    void Close(void);
    FileWithErrors(void);
    ~FileWithErrors(void);
};

#define MIN_READ_SIZE 32
void FileWithErrors::Read(UBYTE *pBuffer, USHORT cAsk)
{
    // Put buffer in known state, and calculate
    // ending offset.
    //
    memset(pBuffer, 0, cAsk);
    LONG NewOffset = Offset + cAsk;

    // Verify file handle
    //
    if (hFile == -1)
 { Offset = NewOffset; return; }

    // Verify file offset and length of request
    //
    if (Offset >= Length)
 { Offset = NewOffset; return; }
    if (NewOffset > Length)
    {
        cAsk = Length - Offset;
    }

    // Attempt progressively smaller read requests. Start by asking for
    // cAsk bytes at Offset. If that doesn't work, then ask in smaller
    // requests.
    //
    USHORT cRequest = cAsk;
    USHORT cTotal = 0;
    while (cTotal != cAsk)
    {
        // Attempt a seek
        //
        LONG lcc = _lseek(hFile, Offset, SEEK_SET);
        if (lcc != Offset) { Offset = NewOffset; return; }
    
        // Attempt a read. If EOF, then we got as many bytes as we could
        // anyway. Sucess speaks for itself. Only errors need further
        // processing.
        //
        USHORT cGiven;
        if (cTotal + cRequest > cAsk) cRequest = cTotal - cAsk;
        int cc =  ReadFromFile(hFile, pBuffer, cRequest, cGiven);
        if (cc == -1)
        {
            // Read Error
            //
            cRequest /= 2;
            if (cRequest < MIN_READ_SIZE)  { Offset = NewOffset; return; }
        }
        else
        {
            cTotal += cRequest;
            pBuffer += cRequest;
            Offset += cRequest;
        }
    }
}

void FileWithErrors::Write(UBYTE *pBuffer, USHORT cAsk)
{
    // Calculate ending offset.
    //
    long NewOffset = Offset + cAsk;

    // Verify file handle
    //
    if (hFile == -1) { Offset = NewOffset; return; }

    // Attempt a seek
    //
    LONG lcc = _lseek(hFile, Offset, SEEK_SET);
    if (lcc != Offset) { Offset = NewOffset; return; }

    // Attempt a write. If FULL DISK, then we wrote as many bytes as we
    // could anyway. Sucess speaks for itself. Only errors need further
    // processing.
    //
    USHORT cTaken;
    WriteToFile(hFile, pBuffer, cAsk, cTaken);

BumpOffset:
    Offset = NewOffset;
}

FileWithErrors::FileWithErrors(void)
{
    hFile = -1; // File not opened
    Length = 0;
    Offset = 0;
}

FileWithErrors::~FileWithErrors(void)
{
    if (hFile != -1) Close();
}

void FileWithErrors::OpenForRead(char *filename)
{
    Length = 0;
    Offset = 0;
    hFile = _sopen(filename, _O_BINARY+_O_RDONLY, _SH_DENYRW);
    if (hFile == -1) return;
    Length = _lseek(hFile, 0, SEEK_END);
    if (Length < -1L) Length = 0;
    _lseek(hFile, 0L, SEEK_SET);
}

void FileWithErrors::OpenForWrite(char *filename)
{
    Length = 0;
    Offset = 0;
    hFile = _sopen(filename, _O_BINARY+_O_WRONLY+_O_TRUNC, _SH_DENYRW);
}

void FileWithErrors::Close(void)
{
    if (hFile == -1) return;
    _close(hFile);
    hFile = -1;
}

#pragma pack(push,1)
typedef struct
{
   UBYTE HeaderSize;
   _int64  Length;
} Header;
#pragma pack(pop)

#if 1
void main(void)
{
    CalculateTable(); // CRC-32 table
    generate_gf();
    b0 = (kk+1)/2;
    gen_poly();
    gen_table();

    UBYTE UserData[57232];
    SECTION Sections[NUM_SECTIONS];
    int Erasures[NUM_SECTIONS];
    for (int i = 0; i < 57232; i++)
    {
        UserData[i] = (UBYTE)(i % 196);
    }
    for (i = 0; i < 100; i++)
    {
        EncodeSegment(Sections, UserData);
        Sections[i%NUM_SECTIONS].Data[i%SIZEOF_SECTION_DATA] = i;
        for (int i = 0; i < NUM_SECTIONS; i++)
        {
            Erasures[i] = FALSE;
        }
        if (!DecodeSegment(UserData, Sections, Erasures))
        {
            printf("Decode Failure for pass %d.\n", i);
        }
        for (int j = 0; j < 57232; j++)
        {
            if (UserData[j] != (UBYTE)(j % 196))
            {
                printf("Didn't decode segment correctly. %x != %x at position %d\n", UserData[j], (UBYTE)(j % 196), j);
                exit(-1);
            }
        }
    }
    //DisplaySegment(Sections);
}
#elif 0
// Encode
void main(void)
{
    CalculateTable(); // CRC-32 table
    generate_gf();
    b0 = (kk+1)/2;
    gen_poly();
    gen_table();

    UBYTE UserData[SIZEOF_CLIENT_DATA];
    int Erasures[NUM_SECTIONS];
    int hInput = _sopen("input.dat", _O_BINARY+_O_RDONLY, _SH_DENYRW);
    if (hInput == -1) { printf("Couldn't open input.dat\n"); exit(-1); }
    LONG Length = _lseek(hInput, 0, SEEK_END);
    _lseek(hInput, 0, SEEK_SET);

    int hOutputs[NUM_SECTIONS];
    for (int i = 0; i < NUM_SECTIONS; i++)
    {
        char filename[100];
        sprintf(filename, "output_%x.dat", i);
        hOutputs[i] = _sopen(filename, _O_BINARY+_O_WRONLY+_O_TRUNC+_O_CREAT, _SH_DENYRW, 0755);
        if (hOutputs[i] == -1)
        {
            printf("Couldn't open %s.\n", filename);
            exit(-1);
        }
    }
    
    ULONG SegmentNumber = 0;
    USHORT length;
    BOOL isFirstBlock = TRUE;
    do
    {
        int cRequest;
        UBYTE *pUserData;
        memset(UserData, 0, SIZEOF_CLIENT_DATA);
        if (isFirstBlock)
        {
            Header Hdr;
            Hdr.HeaderSize = sizeof(Header);
            Hdr.Length = Length;
            cRequest = SIZEOF_CLIENT_DATA - sizeof(Hdr);
            pUserData = UserData + sizeof(Hdr);
            *(Header *)(UserData) = Hdr;
            isFirstBlock = FALSE;
        }
        else
        {
            cRequest = SIZEOF_CLIENT_DATA;
            pUserData = UserData;
        }
        ReadFromFile(hInput, pUserData, cRequest, length);
        Length -= length;
        EncodeSegment(SegmentNumber++, UserData);
        for (i = 0; i < NUM_SECTIONS; i++)
        {
            USHORT cWritten;
            WriteToFile(hOutputs[i], (UBYTE *)(Section+i), sizeof SECTION, cWritten);
            if (cWritten != sizeof SECTION)
            {
                printf("Couldn't write.\n");
                exit(-1);
            }
        }
    } while (Length);
}
#else
// Decode
void main(void)
{
    CalculateTable(); // CRC-32 table
    generate_gf();
    b0 = (kk+1)/2;
    gen_poly();
    gen_table();

    UBYTE UserData[SIZEOF_CLIENT_DATA];
    int Erasures[NUM_SECTIONS];
    int hOutput = _sopen("output.dat", _O_BINARY+_O_WRONLY+_O_TRUNC, _SH_DENYRW);
    if (hOutput == -1)
    {
        printf("Couldn't open output.dat\n");
        exit(-1);
    }

    int hInputs[NUM_SECTIONS];
    for (int i = 0; i < NUM_SECTIONS; i++)
    {
        char filename[100];
        sprintf(filename, "input_%x.dat", i);
        hInputs[i] = _sopen(filename, _O_BINARY+_O_RDONLY, _SH_DENYRW);
        if (hInputs[i] == -1)
        {
            printf("Couldn't open %s.\n", filename);
            exit(-1);
        }
    }
    
    USHORT length;
    ULONG Length;
    BOOL isFirstBlock = TRUE;
    UBYTE *pUserData;
    do
    {
        for (i = 0; i < NUM_SECTIONS; i++)
        {
            USHORT cRead;
            ReadFromFile(hInputs[i], (UBYTE *)(Section+i), sizeof SECTION, cRead);
            if (cRead != sizeof SECTION)
            {
                printf("Problems reading....\n");
            }
            Erasures[i] = FALSE;
        }
        //Erasures[0] = TRUE;
        memset(UserData, 0, SIZEOF_CLIENT_DATA);
        if (!DecodeSegment(UserData, Section, Erasures))
        {
            printf("Decode Failure.\n");
            exit(-1);
        }
        if (isFirstBlock)
        {
            Header Hdr = *(Header *)UserData;
            if (Hdr.HeaderSize != sizeof(Header))
            {
                printf("Unknown Header format.\n");
                exit(-1);
            }
            if (Hdr.Length > 0xFFFFFFFFUL)
            {
                printf("64bit file is too large for me.\n");
                exit(-1);
            }
            Length = Hdr.Length;
            length = SIZEOF_CLIENT_DATA - sizeof(Header);
            pUserData = UserData + sizeof(Header);
            isFirstBlock = FALSE;
        }
        else
        {
            length = SIZEOF_CLIENT_DATA;
            pUserData = UserData;
        }
        USHORT Written;
        printf("Writting %d bytes.\n", length);
        WriteToFile(hOutput, pUserData, length, Written);
        if (Written != length)
        {
            printf("Problems writting.\n");
            exit(-1);
        }
        Length -= length;
    } while (Length);
}

#endif
