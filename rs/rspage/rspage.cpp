#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <memory.h>
#include <fcntl.h>

#include "projdefs.h"
#include "rs_factory.h"
#include "crc.h"
#include "base64.h"

size_t getline(char** line, size_t* length, FILE* stream);

bool rs_fopen(FILE** pFile, const char* filename, const char* mode)
{
    if (pFile)
    {
        *pFile = nullptr;
        if (nullptr != filename
            && nullptr != mode)
        {
#if !defined(__INTEL_COMPILER) && (_MSC_VER >= 1400)
            // 1400 is Visual C++ 2005
            //
            return (fopen_s(pFile, filename, mode) == 0);
#else
            * pFile = fopen((char*)filename, (char*)mode);
            if (nullptr != *pFile)
            {
                return true;
            }
#endif // WINDOWS_FILES
        }
    }
    return false;
}

void usage(char* program_name)
{
    fprintf(stderr, "Usage: %s [-e|-d] infile outfile\r\n", program_name);
    exit(0);
}

int mode = 0;
#define ENCODING 1
#define DECODING 2

char out_file_buffer[16 * 1024];
char in_file_buffer[16 * 1024];

int main(int argc, char* argv[])
{
    if (argc != 4) usage(argv[0]);
    if (argv[1][0] != '-') usage(argv[0]);
    if (argv[1][2] != '\0') usage(argv[0]);
    if (argv[1][1] == 'e' || argv[1][1] == 'E')
    {
        mode = ENCODING;
    }
    else if (argv[1][1] == 'd' || argv[1][1] == 'D')
    {
        mode = DECODING;
    }
    else usage(argv[0]);

    FILE* fp_in = nullptr;
    if (!rs_fopen(&fp_in, argv[2], "rb") || nullptr == fp_in)
    {
        return 0;
    }
    setvbuf(fp_in, in_file_buffer, _IOFBF, sizeof(in_file_buffer));
    FILE* fp_out = nullptr;
    if (!rs_fopen(&fp_out, argv[3], "wb") || nullptr == fp_out)
    {
        fclose(fp_in);
        return 0;
    }
    setvbuf(fp_out, out_file_buffer, _IOFBF, sizeof(out_file_buffer));

    const auto p_encoder = RS_FACTORY::instance().create_specialized_codec(16, 0);
    if (mode == ENCODING)
    {
        // Initialize
        //
        int row_number = 0;
        UBYTE input_block[39 * 27] = {};
        UBYTE check[39 * 32] = {};
        size_t cc = fread(input_block, 1, sizeof input_block, fp_in);
        while (cc > 0)
        {
            // Do a block
            //
            // 1. Encode 39 rows of 27 bytes with a RS(59,27) code producing 39 rows of 59 bytes.
            //
            UBYTE step1_block[39 * 59] = {};
            for (int i = 0; i < 39; i++)
            {
                UBYTE codeword[255] = {};
                memcpy(codeword + 255 - 32 - 27, input_block + i * 27, 27);
                p_encoder->RSEncode(codeword, codeword + 255 - 32);
                memcpy(step1_block + i * 59, codeword + 255 - 59, 59);
            }

            // 2. Rotate 39 rows of 59 bytes (59x39) to produce 59 rows of 39 bytes (39x59).
            //
            UBYTE step2_block[59 * 39] = {};
            for (int i = 0; i < 59; i++)
            {
                for (int j = 0; j < 39; j++)
                {
                    step2_block[i * 39 + j] = step1_block[j * 59 + i];
                }
            }

            // 3. Append one (1) byte to each row to produce 59 rows of 40 bytes (40x59).
            //
            UBYTE step3_block[59 * 40] = {};
            for (int i = 0; i < 59; i++)
            {
                step3_block[i * 40 + 39] = static_cast<UBYTE>(row_number++);
                memcpy(step3_block + i * 40, step2_block + i * 39, 39);
            }

            // 4. Encode 59 rows of 40 bytes (40x59) with a RS(72,40) into 59 rows of 72 bytes (72x59).
            //
            UBYTE step4_block[59 * 72] = {};
            for (int i = 0; i < 59; i++)
            {
                UBYTE codeword[255] = {};
                memcpy(codeword + 255 - 32 - 40, step3_block + i * 40, 40);
                p_encoder->RSEncode(codeword, codeword + 255 - 32);
                memcpy(step4_block + i * 72, codeword + 255 - 72, 72);
            }                        

            // 5. Encode 59 rows of 72 bytes (72x59) using Base-64 to produce 59 rows of 96 characters (96x59).
            //
            for (int i = 0; i < 59; i++)
            {
                auto s = base64_encode(step4_block + i * 72, 72);
                fwrite(s.c_str(), 1, s.length(), fp_out);
                fwrite("\r\n", 1, 2, fp_out);
            }

            memset(input_block, 0, sizeof input_block);
            cc = fread(input_block, 1, sizeof input_block, fp_in);
        }
    }
    else if (DECODING == mode)
    {
        char* line = nullptr;
        size_t len = 0;
        while ((getline(&line, &len, fp_in)) != -1)
        {
            auto s = std::string(line);
            if (s.length() == 0) continue;

            if (s.length() != 96)
            {
                printf("Not 96 characters long.\n");
            }

            // Decode 96-character line to 72 bytes.
            //
            auto t = base64_decode(s, true);
            if (t.length() != 72)
            {
                printf("Not 72 characters long.\n");
            }

            UBYTE codeword[255] = {};
            memcpy(codeword + 255 - 72, t.c_str(), 72);
            if (p_encoder->RSDecode(codeword) < 0)
            {
                printf("Too many errors. Rewrite code to handle erasures.");
            }
        }
        delete[] line;
    }

    fclose(fp_out);
    fclose(fp_in);
}
