// Copyright (C) 1996-1997 Solid Vertical Domains, Ltd. All rights reserved.
//
// Solid Vertical Domains, Ltd.
// 11410 NE 124th ST #162, Kirkland WA 98034-4399 USA
// (206) 488-0378
// Stephen V. Dennis <SDennis@accessone.com>
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rs.h"

#if 1
long int nrand48(unsigned short seed[3])
{
    return rand();
}

double erand48(unsigned short seed[3])
{
    return ((double)rand())/0xFFFFUL;
}
#endif

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

int mode = 0;
#define ENCODING 1
#define DECODING 2

RS_ENCODER *prs;

int main(int argc, char *argv[])
{
    if (argc == 3 && argv[1][0] == '-' && (argv[1][1] == 'e' || argv[1][1] == 'E'))
    {
        mode = ENCODING;
        tt = atoi(argv[2]);
        kk = nn - 2*tt;
    }
    else if (argc == 3 && argv[1][0] == '-' && (argv[1][1] == 'd' || argv[1][1] == 'D'))
    {
        mode = DECODING;
        tt = atoi(argv[2]);
        kk = nn - 2*tt;
    }
    else
    {
        printf("Usage: rstest [-e|-d] tt\n");
        exit(0);
    }

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
    /**** Storage for random seed ****/
    unsigned short int store[3];
    double x, prob_symb_err;

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
        no_cws = 1000000;
        for (i=0;i < 3;i++)
        {
            store[i] = 109-i*5; // store[i] = tmp;
        }
      
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
                data[j] = (int) (nrand48(store) % 256);
    
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
        // printf("Enter length of the shortened code\n");
        nn_short = 255; // scanf("%d",&nn_short);
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

        no_cws = 1000000; // printf("Enter no of codewords \n"); scanf("%d",&no_cws);

        if (rs_fopen(&fp, file_in, "rb") || nullptr == fp)
        {
            printf("Could not open %s\n", file_in);
            exit(0);
        }
    
        prob_symb_err = (double)tt/ (8.0 * (double)nn_short);
        /* CHANGE above probablity of symbol error if desired */
        /* printf("Probability of symbol error = %6e\n",prob_symb_err);*/
    
        no_decoder_errors = 0;
        iter = -1;
        while (++iter < no_cws)
        {
            fread(cw, sizeof(GF), nn_short, fp);
#if 0
            /* read in codewords from file */
            for (i=0;i < nn_short;i++)
            {
                fscanf(fp,"%02x ",&byte);
                cw[i] = byte;
            }
            /* printf("The Transmitted Codeword\n");
            for (i=0;i < nn_short;i++) printf("%02x ",cw[i]); printf("\n");*/
#endif
        
#ifndef NO_PRINT
            printf("\n\n\n\n Transmitting codeword %d \n",iter);
            printf("Channel caused following errors Location (Error): \n");
#endif
            no_ch_errs = 0;
            for (i=0;i < nn_short;i++)
            {
#if 1
                x = erand48(store);
                if (x < prob_symb_err)
                {
                    error_byte = (int) (nrand48(store) % 255) + 1;
                    received[i] = cw[i] ^ error_byte;
                    ++no_ch_errs;
#ifndef NO_PRINT
                    printf("%d (%#x) ", i, error_byte);
#endif
                }
                else
#endif
                    received[i] = cw[i];
            }
#ifndef NO_PRINT
            printf("\n");
            printf("Channel caused %d errors\n", no_ch_errs);
            /* for (i=0;i < nn_short;i++) printf("%02x ",received[i]); printf("\n");*/
#endif

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
#ifndef NO_PRINT
            printf("decode_rs() returned %d\n", decode_flag);
#endif
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
    
        } /* closing iteration loop */
        if (no_decoder_errors > 0)
            printf("The number of decoder errors were %d\n",no_decoder_errors);
        else
            printf("Decoder corrected all occurrences of %d or less errors\n", tt);
    }
    return 0;
}
