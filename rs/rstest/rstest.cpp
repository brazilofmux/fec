// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All rights reserved.
// Stephen V. Dennis <brazilofmux@gmail.com>
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <cstdint>
#include "rs.h"

enum class Verbosity {
    Quiet,      // Summary only
    Normal,     // Basic progress
    Debug       // Full details (old behavior)
};

struct RsTestConfig {
    int tt;                     // Error correction capacity
    int num_codewords;          // Default 10000
    Verbosity verbose_level;
    unsigned int random_seed;   // For reproducibility
    double error_rate;          // Override default tt/(8*nn_short)
    bool verify_correction;     // Extra validation
};

struct RsTestFileHeader {
    uint32_t magic;            // Identifier
    uint32_t version;          // Format version
    uint32_t tt;               // Error correction capacity
    uint32_t num_codewords;    // How many vectors follow
    uint32_t codeword_size;    // Total size (parity + data)
};

struct DecodingStats {
    int total_codewords;
    int total_errors_injected;
    int errors_corrected;
    int decode_failures;
    int errors_by_count[MAX_TT + 1];  // Distribution of error counts
    int uncorrectable_errors;         // When errors > tt
};

void print_progress(const RsTestConfig& config, const char* msg, ...) {
    if (config.verbose_level == Verbosity::Quiet) {
        return;
    }

    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
}

void print_debug(const RsTestConfig& config, const char* msg, ...) {
    if (config.verbose_level != Verbosity::Debug) {
        return;
    }

    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
}

long int nrand48()
{
    return rand();
}

double erand48()
{
    return ((double)rand())/((double)RAND_MAX);
}

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

int tt;
int kk;
GF data[MAX_KK];
GF recd[nn];
GF bb[2*MAX_TT];
RS_ENCODER *prs;

int mode = 0;
#define ENCODING 1
#define DECODING 2

RsTestConfig parse_args(int argc, char* argv[]) {
    RsTestConfig config = {
        .tt = 0,
        .num_codewords = 10000,
        .verbose_level = Verbosity::Normal,
        .random_seed = 1093,
        .error_rate = 0.0,
        .verify_correction = true
    };

    if (argc < 3) {
        printf("Usage: rstest [-e|-d] tt [-v level] [-n count] [-s seed] [-r rate]\n");
        printf("  -e/-d tt    : encode/decode mode and error correction capacity\n");
        printf("  -v level    : verbosity (0=quiet, 1=normal, 2=debug)\n");
        printf("  -n count    : number of codewords (default: 10000)\n");
        printf("  -s seed     : random seed (default: 1093)\n");
        printf("  -r rate     : error injection rate (default: 0.01)\n");
        exit(1);
    }

    // Parse mode first
    if (argv[1][0] != '-' || (argv[1][1] != 'e' && argv[1][1] != 'E' && 
                             argv[1][1] != 'd' && argv[1][1] != 'D')) {
        printf("First argument must be -e or -d\n");
        exit(1);
    }
    mode = (argv[1][1] == 'e' || argv[1][1] == 'E') ? ENCODING : DECODING;

    // Parse tt
    config.tt = tt = atoi(argv[2]);
    if (config.tt < MIN_TT || config.tt > MAX_TT) {
        printf("tt must be between %d and %d\n", MIN_TT, MAX_TT);
        exit(1);
    }

    // Parse optional arguments
    int v;
    for (int i = 3; i < argc; i++) {
        if (argv[i][0] != '-' || strlen(argv[i]) != 2) {
            printf("Invalid option: %s\n", argv[i]);
            exit(1);
        }

        switch (argv[i][1]) {
            case 'v':
                if (i + 1 >= argc) {
                    printf("Missing value for -v\n");
                    exit(1);
                }
                v = atoi(argv[++i]);
                config.verbose_level = static_cast<Verbosity>(v);
                break;

            case 'n':
                if (i + 1 >= argc) {
                    printf("Missing value for -n\n");
                    exit(1);
                }
                config.num_codewords = atoi(argv[++i]);
                break;

            case 's':
                if (i + 1 >= argc) {
                    printf("Missing value for -s\n");
                    exit(1);
                }
                config.random_seed = atoi(argv[++i]);
                break;

            case 'r':
                if (i + 1 >= argc) {
                    printf("Missing value for -r\n");
                    exit(1);
                }
                config.error_rate = atof(argv[++i]);
                break;

            default:
                printf("Unknown option: -%c\n", argv[i][1]);
                exit(1);
        }
    }

    return config;
}

int main(int argc, char *argv[])
{
    RsTestConfig config = parse_args(argc, argv);
    kk = nn - 2*tt;

    RS_Init();

    /* Generate the Galois field GF(2^8) */
    // printf("Generating the Galois field GF(%d)\n\n", n);
    prs = new RS_ENCODER(tt);

    GF error_byte;
    GF received[1000], cw[1000];
    int i, j, no_cws, nn_short, kk_short, no_ch_errs, iter, error_flag;
    int decode_flag, no_decoder_errors;
    long seed;
    const char* cw_file = "cws.txt";
    const char* file_in = "cws.txt";;
    FILE *fp_out, *fp;
    double x;

#ifdef DEBUG
    printf("Look-up tables for GF(2^8)\n");
    printf("i\t\t Pow2Poly[i]\t Poly2Pow[i]\n");
    for (i=0; i < n; i++)
        printf("%3d \t\t %#x \t\t %d\n",i,Pow2Poly[i],Poly2Pow[i]);
    printf("\n The polynomial-form representation of @^i, 0 <= i < %d \n",n);
    printf("is obtained as follows: It is given by the binary representation of \n");
    printf("the integer Pow2Poly[i] with LS-Bit corresponding to @^0\n");
    printf("and the next bit to @^1 and so on.\n");
    printf("\n The power of a Galois field element x  whose polynomial-form \n");
    printf("representation is the binary representation of integer \n");
    printf("j, 0 <= j < %d is given by the integer Poly2Pow[j]\n",n);
    printf("Therefore, @^j = x \n");
    printf("\n EXAMPLES IN GF(256) \n");
    printf("The polynomial-form of @^8 is the binary representation of the integer Pow2Poly[8] viz. 0x1d. i.e., @^8 = @^4 + @^3 + @^2 + @^0 \n");
    printf("The power of x whose polynomial-form is 0x72 (integer 155) is Poly2Pow[155] = 217, so @^217 = x \n");
#endif

#ifdef DEBUG
    printf("\n g(x) of the %d-error correcting RS (%d,%d) code: \n",tt,nn,kk);
    printf("Coefficient in power form with @ being the primitive element\n");
    for (i = 0; i <= nn-kk; i++)
    {
        printf("%3d x^%02d ", gg[i], i);
        if (i < nn-kk)
            printf("+ ");
        if (i && ((i+1) % 6 == 0))
            printf("\n"); /* 7 coefficients per line */
    }
    printf("\n");

    printf("\nCoefficient in polynomial form {@^7,...,@,1}\n");
    for (i = 0; i <= nn-kk; i++)
    {
        printf("0x%02x x^%02d ", Pow2Poly[gg[i]], i);
        if (i < nn-kk)
            printf("+ ");
        if (i && ((i+1) % 6 == 0))
            printf("\n"); /* 7 coefficients per line */
    }
    printf("\n\n");
#endif

    if (ENCODING == mode)
    {
        nn_short = 255;
        if ((nn_short < 2*tt) || (nn_short > nn))
        {
            printf("Invalid entry %d for shortened length\n",nn_short);
            exit(0);
        }

        /* compute dimension of shortened code */
        kk_short = kk - (nn-nn_short);
        no_cws = 10000;

#ifndef PROFILE
    	if (!rs_fopen(&fp_out, cw_file, "wb") || nullptr == fp_out)
        {
            printf("Could not open %s\n", cw_file);
            exit(0);
        }
#endif

        /**** BEGIN: Encoding random information vectors ****/
        /* Pad with zeros rightmost (kk-kk_short) positions */
        memset(data + kk_short, 0, kk - kk_short);
        for (i = 0; i < no_cws; i++)
        {
            for ( j = 0; j < (int)kk_short; j++)
                data[j] = (int) (nrand48() % 256);

            prs->RSEncode(data, bb);

#ifndef PROFILE
            fwrite(bb, sizeof(GF), nn_short - kk_short, fp_out);
            fwrite(data, sizeof(GF), kk_short, fp_out);
#if 0
            for (j=0;j < (int)nn_short-(int)kk_short;j++) /* parity bytes first */
                fprintf(fp_out,"%02x ",bb[j]);
            for (j=0;j < (int)kk_short;j++) /* info bytes last */
                fprintf(fp_out,"%02x ",data[j]);
            fprintf(fp_out,"\n");
#endif
#endif
        }
#ifndef PROFILE
        fclose(fp_out);
#endif
    }
    else if (DECODING == mode)
    {
        DecodingStats stats = {
            .total_codewords = 0,
            .total_errors_injected = 0,
            .errors_corrected = 0,
            .decode_failures = 0,
            .uncorrectable_errors = 0
        };
        memset(stats.errors_by_count, 0, sizeof(stats.errors_by_count));

        config.error_rate = config.error_rate > 0 ? config.error_rate : 0.01;
        print_debug(config, "Error rate set to: %f\n", config.error_rate);

        nn_short = 255;
        if ((nn_short < 2*tt) || (nn_short > nn))
        {
            printf("Invalid entry %d for shortened length\n",nn_short);
            exit(0);
        }
        kk_short = kk - (nn-nn_short); /* compute dimension of shortened code */
        //printf("The (%d,%d) %d-error correcting RS code\n",nn_short,kk_short,tt);
        //printf("Enter a random positive integer\n");
        seed = 1093; //scanf("%ld",&seed);
        srand(seed);

        no_cws = 10000; // printf("Enter no of codewords \n"); scanf("%d",&no_cws);

        if (!rs_fopen(&fp, file_in, "rb") || nullptr == fp)
        {
            printf("Could not open %s\n", file_in);
            exit(0);
        }

        no_decoder_errors = 0;
        iter = -1;
        while (++iter < no_cws) {
            stats.total_codewords++;
            fread(cw, sizeof(GF), nn_short, fp);

            print_debug(config, "\n\n\n\n Transmitting codeword %d \n", iter);
            print_debug(config, "Channel caused following errors Location (Error): \n");
            no_ch_errs = 0;
            for (i=0;i < nn_short;i++)
            {
                x = erand48();
                print_debug(config, "Random value: %f vs threshold: %f\n", x, config.error_rate);
                if (x < config.error_rate)
                {
                    error_byte = (int) (nrand48() % 255) + 1;
                    print_debug(config, "Injecting error: %d\n", error_byte);
                    received[i] = cw[i] ^ error_byte;
                    ++no_ch_errs;
                }
                else
                    received[i] = cw[i];
            }
            stats.total_errors_injected += no_ch_errs;
            stats.errors_by_count[no_ch_errs]++;
            print_debug(config, "Channel caused %d errors\n", no_ch_errs);

            // Pad with zeros and decode as if was a (255,kk) tt-error correcting
            // code.
            //
            for (i=0;i < nn-kk;i++)
                recd[i] = received[i]; /* parity bytes */
            for (i=nn-kk;i < nn-kk+kk_short;i++)
                recd[i] = received[i]; /* info bytes */
            for (i=nn-kk+kk_short;i < nn;i++)
                recd[i] = 0; /* zero padding */

#if 1
            decode_flag = prs->RSDecode(recd);
#else
            decode_flag = rs.RSDecodeErasures(recd, 0, 0);
#endif
            print_debug(config, "decode_rs() returned %d\n", decode_flag);
            error_flag = 0;
            for (i=0; i < kk_short; i++)
            {
                if (recd[i+nn-kk] != cw[i+nn-kk])
                {
                    if (no_ch_errs <= tt)
                    {
                        printf("Position %d miscorrected %x,%x=%x.\n", i+nn-kk,
                        recd[i+nn-kk], cw[i+nn-kk], recd[i+nn-kk]^cw[i+nn-kk]);
                    }
                    error_flag = 1;
                    //break;
                }
            }

            if (decode_flag == 0) {
                stats.errors_corrected++;
            }
            if (error_flag == 1 && no_ch_errs <= tt) {
                stats.decode_failures++;
            }
            if (no_ch_errs > tt) {
                stats.uncorrectable_errors++;
            }        

            if (decode_flag == RS_ERROR_LAMDA_ERROR && no_ch_errs <= tt)
            {
                printf("%d ch errs  <=  max # correctable errs but \n", no_ch_errs);
                printf("\n DECODER ERROR: deg(lambda) unequal to #roots \n");
                exit(2); /* DECODER ERROR condition */
            }
            else if (decode_flag == RS_ERROR_CHIEN_SEARCH && no_ch_errs <= tt)
            {
                printf(" %d ch errs  <= max # correctable errs but \n", no_ch_errs);
                printf("\n deg(lambda) > 2*tt \n");
                exit(3);/* DECODER ERROR condition */
            }
            if (error_flag == 1 && no_ch_errs <= tt)
            {
                printf("DECODER FAILED TO CORRECT %d ERRORS\n", no_ch_errs);
                ++no_decoder_errors;
#ifndef NO_PRINT
                printf("The Transmitted Codeword\n");
                for (i = 0; i < nn_short; i++)
                    printf("%02x ", cw[i]);
                printf("\n");
#endif
            }

        }
        /* closing iteration loop */
        print_progress(config, "\nDecoding Summary:\n");
        print_progress(config, "Total codewords processed: %d\n", stats.total_codewords);
        print_progress(config, "Total errors injected: %d (avg %.2f per codeword)\n", 
            stats.total_errors_injected, 
            (double)stats.total_errors_injected / stats.total_codewords);
        print_progress(config, "Successful corrections: %d\n", stats.errors_corrected);
        if (stats.decode_failures > 0) {
            print_progress(config, "Decoder failures: %d\n", stats.decode_failures);
        }
        print_progress(config, "Uncorrectable error patterns: %d\n", stats.uncorrectable_errors);
    }
    return 0;
}
