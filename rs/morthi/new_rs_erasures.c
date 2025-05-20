/*
 * File:     general_RS.c
 * Title:    Encoder/decoder for RS codes in C
 * Authors:  Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu)
 *           and Hari Thirumoorthy (harit@spectra.eng.hawaii.edu)
 * Date:     Aug 1995
 *
 *
 * ===============  Encoder/Decoder for RS codes in C =================
 *
 *
 * The encoding and decoding methods used in this program are based on the
 * book "Error Control Coding: Fundamentals and Applications", by Lin and
 * Costello, Prentice Hall, 1983. Most recently on "Theory and Practice of 
 * Error Control Codes", by R.E. Blahut.
 *
 *
 * NOTE:    
 *          The authors are not responsible for any malfunctioning of
 *          this program, nor for any damage caused by it. Please include the
 *          original program along with these comments in any redistribution.
 *
 * Portions of this program are from a Reed-Solomon encoder/decoder
 * in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89.
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may 
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * Copyright (c) 1995.  Robert Morelos-Zaragoza and Hari Thirumoorthy.
 *                      All rights reserved.
 *
 */

/*
  Program computes the generator polynomial of a RS code. Also performs 
encoding and decoding of the RS code or a shortened RS code. Compile
using one of the following options:

  In the commands below, using  the [-DNO_PRINT] option will ensure that
 the program runs without printing a lot of detail that is generally
 unnecessary. If the program is compiled without using the option a lot
 of detail will be printed. Some other options, for example changing the
 probability of symbol error are available by directly editing the program 
 itself. 

 1. cc general_RS.c -DNO_PRINT -lm -o rsgp
        "rsgp" computes the generator polynomial of a RS code whose
        length is "nn", dimension is "kk" and error correcting capability
        is "tt". The generator polynomial is displayed in two formats.
        In one the coefficients of the generator polynomial are given
        as exponents of the primitive element @ of the field. In the other
        the coefficients are given by their representations in the
        basis {@^7,...,@,1}.

 2. cc general_RS.c -DDEBUG -DNO_PRINT -lm -o rsfield
        "rsfield" does all that rsgp does and in addition it also
        displays the field GF(n) and the two representations of the
        field elements namely the "power of @" representation and
        the "Basis" representation.

 3. cc general_RS.c -DENCODE -DNO_PRINT -lm -o rsenc
        "rsenc" does all that rsgp does and in addition it prompts the
        user for the following data. (1) The length of the shortened
        RS code. Of course, if shortening is not desired then, the integer
        nn is entered as the length. (2) The number of codewords that are
        to be generated. (3) The name of the file into which they are to
        be stored. (4) Random number seeds that are needed to ensure true
        randomness in generation of the codewords. (The DowJones closing
        index is a good number to enter!).  These codewords can then be
        transmitted over a noisy channel and decoded using the rsdec program.

        The error correcting capability of the code cannot be set by running
        the program. It must be set in this file "general_RS.c" by changing
        the "tt" parameter. Each time the "tt" is
        changed, the dimension "kk" of the unshortened RS code changes and
        this is given by kk = nn - 2 * tt. This value of kk must be set in
        the program below. After both "tt" and "kk" are set, the program must
        be recompiled.

 
 4.  cc general_RS.c -DDECODE -DNO_PRINT -lm -o rsdec
        "rsdec" does all that rsgp does and in addition it prompts the user
        for the following data. (1) The length of the RS code. (2) The
        file from which the codewords are to be read. (3) The number of code
        words to read. (4) The probability of symbol error. The program then
        reads in codewords and transmits each codeword over a discrete
        memoryless channel that has a symbol error probability specified
        above. The received corrupted word is then decoded and the results
        of decoding displayed. This basically verifies that the decoding
        procedure is capable of producing the transmitted codeword provided
        the channel caused at most tt errors. Almost always when the channel
        causes more than tt errors, the decoder fails to produce a codeword.
        In such a case the decoder outputs the uncorrected information bytes
        as they were received from the channel.

5.  cc general_RS.c -DDECODE -DDECODER_DEBUG -lm -o rsddebug
        "rsddebug" is similar to rsdec and can be used to debug the
        the decoder and verify its operation. It reads in one codeword
        from a file specified by the user and allows the user to specify
        the number of errors to be caused and their locations. Then the
        decoder transmits the codeword and causes the specified number
        of errors at the specified locations. The actual value of the
        error byte is randomnly generated. The received word is the sum
        of the transmitted codeword and the error word. The decoder then
        decodes this received word and displays all the intermediate
        results and the final result as explained below.
 
        (1) The syndrome of the received word is shown. This consists of
        2*tt bytes {S(1),S(2), ... ,S(2*tt)}. Each S(i) is an element in
        the field GF(n). NOTE: Iff all the S(i)'s are zero, the received
        word is a codeword in the code. In such a case, the channel
        presumably caused no errors.
 
  INPUT:
	Change the Galois Field and error correcting capability by 
	setting the appropriate global variables below including
	1. primitive polynomial 2. error correcting capability 
	3. dimension 4. length.
  
  FUNCTIONS:
		generate_gf() generates the field.
		gen_poly() generates the generator polynomial.
		encode_rs() encodes in systematic form.
		decode_rs() errors-only decodes a vector assumed to be encoded in systematic form.
	        eras_dec_rs() errors+erasures decoding of a vector assumed to be encoded in 
		 	      systematic form.
*/

#include <math.h>
#include <stdio.h>
#define maxrand 0x7fffffff
#define mm  8           /* RS code over GF(2**mm) - change to suit */
#define n   256   	/* n = size of the field */
#define nn  255         /* nn=2**mm -1   length of codeword */
#define tt  16          /* number of errors that can be corrected */
#define kk  223         /* kk = nn-2*tt  */ /* Degree of g(x) = 2*tt */

/**** Primitive polynomial ****/
int pp [mm+1] = { 1, 0, 1, 1, 1, 0, 0, 0, 1}; /* 1+x^2+x^3+x^4+x^8 */

/* int pp [mm+1] = { 1, 1, 0, 0, 1}; */ /* 1 + x + x^4 */

/* int pp[mm+1] = {1,1,0,1}; */ /* 1 + x + x^3 */

/* generator polynomial, tables for Galois field */
int alpha_to [nn+1], index_of [nn+1], gg [nn-kk+1];
int b0;  /* g(x) has roots @**b0, @**(b0+1), ... ,@^(b0+2*tt-1) */

/* data[] is the info vector, bb[] is the parity vector, recd[] is the 
  noise corrupted received vector  */
int recd[nn], data[kk], bb[nn-kk] ;

/* generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
   lookup tables:  index->polynomial form   alpha_to[] contains j=alpha**i;
                   polynomial form -> index form  index_of[j=alpha**i] = i
   alpha=2 is the primitive element of GF(2**m)
   HARI's COMMENT: (4/13/94) alpha_to[] can be used as follows:
        Let @ represent the primitive element commonly called "alpha" that
   is the root of the primitive polynomial p(x). Then in GF(2^m), for any
   0 <= i <= 2^m-2,
        @^i = a(0) + a(1) @ + a(2) @^2 + ... + a(m-1) @^(m-1)
   where the binary vector (a(0),a(1),a(2),...,a(m-1)) is the representation
   of the integer "alpha_to[i]" with a(0) being the LSB and a(m-1) the MSB. Thus for
   example the polynomial representation of @^5 would be given by the binary
   representation of the integer "alpha_to[5]".
                   Similarily, index_of[] can be used as follows:
        As above, let @ represent the primitive element of GF(2^m) that is
   the root of the primitive polynomial p(x). In order to find the power
   of @ (alpha) that has the polynomial representation
        a(0) + a(1) @ + a(2) @^2 + ... + a(m-1) @^(m-1)
   we consider the integer "i" whose binary representation with a(0) being LSB
   and a(m-1) MSB is (a(0),a(1),...,a(m-1)) and locate the entry
   "index_of[i]". Now, @^index_of[i] is that element whose polynomial 
    representation is (a(0),a(1),a(2),...,a(m-1)).
   NOTE:
        The element alpha_to[2^m-1] = 0 always signifying that the
   representation of "@^infinity" = 0 is (0,0,0,...,0).
        Similarily, the element index_of[0] = -1 always signifying
   that the power of alpha which has the polynomial representation
   (0,0,...,0) is "infinity".
 
*/

void generate_gf()
 {
   register int i, mask ;

  mask = 1 ;
  alpha_to[mm] = 0 ;
  for (i=0; i<mm; i++)
   { alpha_to[i] = mask ;
     index_of[alpha_to[i]] = i ;
     if (pp[i]!=0) /* If pp[i] == 1 then, term @^i occurs in poly-repr of @^mm */
       alpha_to[mm] ^= mask ;  /* Bit-wise EXOR operation */
     mask <<= 1 ; /* single left-shift */
   }
  index_of[alpha_to[mm]] = mm ;
  /* Have obtained poly-repr of @^mm. Poly-repr of @^(i+1) is given by 
     poly-repr of @^i shifted left one-bit and accounting for any @^mm 
     term that may occur when poly-repr of @^i is shifted. */
  mask >>= 1 ;
  for (i=mm+1; i<nn; i++)
   { if (alpha_to[i-1] >= mask)
        alpha_to[i] = alpha_to[mm] ^ ((alpha_to[i-1]^mask)<<1) ;
     else alpha_to[i] = alpha_to[i-1]<<1 ;
     index_of[alpha_to[i]] = i ;
   }
  index_of[0] = -1 ;
 }


void gen_poly()
/* Obtain the generator polynomial of the tt-error correcting, length
  nn=(2**mm -1) Reed Solomon code from the product of (X+@**(b0+i)), i = 0, ... ,(2*tt-1)
  Examples: 	If b0 = 1, tt = 1. deg(g(x)) = 2*tt = 2.
 	g(x) = (x+@) (x+@**2)
		If b0 = 0, tt = 2. deg(g(x)) = 2*tt = 4.
	g(x) = (x+1) (x+@) (x+@**2) (x+@**3)  	
*/
 {
   register int i,j ;

   gg[0] = alpha_to[b0];
   gg[1] = 1 ;    /* g(x) = (X+@**b0) initially */
   for (i=2; i <= nn-kk; i++)
    { gg[i] = 1 ;
      /* Below multiply (gg[0]+gg[1]*x + ... +gg[i]x^i) by (@**(b0+i-1) + x) */
      for (j=i-1; j>0; j--)
        if (gg[j] != 0)  gg[j] = gg[j-1]^ alpha_to[((index_of[gg[j]])+b0+i-1)%nn] ;
        else gg[j] = gg[j-1] ;
      gg[0] = alpha_to[((index_of[gg[0]])+b0+i-1)%nn] ;     /* gg[0] can never be zero */
    }
   /* convert gg[] to index form for quicker encoding */
   for (i=0; i <= nn-kk; i++)  gg[i] = index_of[gg[i]] ;
 }


void encode_rs()
/* take the string of symbols in data[i], i=0..(k-1) and encode systematically
   to produce 2*tt parity symbols in bb[0]..bb[2*tt-1]
   data[] is input and bb[] is output in polynomial form.
   Encoding is done by using a feedback shift register with appropriate
   connections specified by the elements of gg[], which was generated above.
   Codeword is   c(X) = data(X)*X**(nn-kk)+ b(X)          */
 {
   register int i,j ;
   int feedback ;

   for (i=0; i<nn-kk; i++)   bb[i] = 0 ;
   for (i=kk-1; i>=0; i--)
    {  feedback = index_of[data[i]^bb[nn-kk-1]] ;
       if (feedback != -1) /* feedback term is non-zero */
        { for (j=nn-kk-1; j>0; j--)
            if (gg[j] != -1)
              bb[j] = bb[j-1]^alpha_to[(gg[j]+feedback)%nn] ;
            else
              bb[j] = bb[j-1] ;
          bb[0] = alpha_to[(gg[0]+feedback)%nn] ;
        }
       else /* feedback term is zero. encoder becomes a single-byte shifter */
        { for (j=nn-kk-1; j>0; j--)
            bb[j] = bb[j-1] ;
          bb[0] = 0 ;
        } ;
    } ;
 }


/* Performs ERRORS-ONLY decoding of RS codes. If decoding is successful,
 writes the codeword into recd[] itself. Otherwise recd[] is unaltered. 
 If channel caused no more than "tt" errors, the tranmitted codeword will
 be recovered.
*/ 

int decode_rs()
 {
   register int i,j,u,q ;
   int elp[nn-kk+2][nn-kk], d[nn-kk+2], l[nn-kk+2], u_lu[nn-kk+2], s[nn-kk+1] ;
   int count=0, syn_error=0, root[tt], loc[tt], z[tt+1], err[nn], reg[tt+1] ;

/* recd[] is in polynomial form, convert to index form */
   for (i=0;i < nn;i++) recd[i] = index_of[recd[i]];

/* first form the syndromes; i.e., evaluate recd(x) at roots of g(x) namely
 @**(b0+i), i = 0, ... ,(2*tt-1) */
   for (i=1; i <= nn-kk; i++)
    { s[i] = 0 ;
      for (j=0; j < nn; j++)
        if (recd[j] != -1)
          s[i] ^= alpha_to[(recd[j]+(b0+i-1)*j)%nn] ;      /* recd[j] in index form */
/* convert syndrome from polynomial form to index form  */
      if (s[i] != 0)  syn_error = 1 ;   /* set flag if non-zero syndrome => error */
      s[i] = index_of[s[i]] ;
    };

#ifdef DECODER_DEBUG
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("The Syndromes of the received vector in polynomial-form are:\n");
    for (i=1;i <= nn-kk;i++){
        if (s[i] == -1)
           printf("s(%d) = %#x ",i,0);
	else
           printf("s(%d) = %#x ",i,s[i]);
    }
    printf("\n NOTE: The basis is {@^%d,...,@,1}\n",mm-1);
#endif   


#ifdef DECODER_DEBUG
    printf("The Syndromes of the received vector in index-form are:\n");
    for (i=1;i <= nn-kk;i++){
        printf("s(%d) = %d ",i,s[i]);
    } 
    printf("\n NOTE: Since 0 is not expressible as a power of the primitive element\n");
    printf("alpha (denoted @),  -1 is printed if the corresponding syndrome is zero.\n");
    printf("---------------------------------------------------------------\n");
#endif

   if (syn_error)       /* if errors, try and correct */
    {
#ifdef DECODER_DEBUG
    printf("Received vector is not a codeword. Channel has caused at least one error\n");
#endif

/* compute the error location polynomial via the Berlekamp iterative algorithm,
   following the terminology of Lin and Costello :   d[u] is the 'mu'th
   discrepancy, where u = 'mu'+ 1 and 'mu' (the Greek letter!) is the step number
   ranging from -1 to 2*tt (see L&C),  l[u] is the
   degree of the elp at that step, and u_l[u] is the difference between the
   step number and the degree of the elp.
  
   The notation is the same as that in Lin and Costello's book; pages 155-156 and 175. 

*/
/* initialise table entries */
      d[0] = 0 ;           /* index form */
      d[1] = s[1] ;        /* index form */
      elp[0][0] = 0 ;      /* index form */
      elp[1][0] = 1 ;      /* polynomial form */
      for (i=1; i<nn-kk; i++)
        { elp[0][i] = -1 ;   /* index form */
          elp[1][i] = 0 ;   /* polynomial form */
        }
      l[0] = 0 ;
      l[1] = 0 ;
      u_lu[0] = -1 ;
      u_lu[1] = 0 ;
      u = 0 ;

      do
      {
        u++ ;
        if (d[u]==-1)
          { l[u+1] = l[u] ;
            for (i=0; i<=l[u]; i++)
             {  elp[u+1][i] = elp[u][i] ;
                elp[u][i] = index_of[elp[u][i]] ;
             }
          }
        else
/* search for words with greatest u_lu[q] for which d[q]!=0 */
          { q = u-1 ;
            while ((d[q]==-1) && (q>0)) q-- ;
/* have found first non-zero d[q]  */
            if (q>0)
             { j=q ;
               do
               { j-- ;
                 if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                   q = j ;
               }while (j>0) ;
             } ;

/* have now found q such that d[u]!=0 and u_lu[q] is maximum */
/* store degree of new elp polynomial */
            if (l[u]>l[q]+u-q)  l[u+1] = l[u] ;
            else  l[u+1] = l[q]+u-q ;

/* form new elp(x) */
            for (i=0; i<nn-kk; i++)    elp[u+1][i] = 0 ;
            for (i=0; i<=l[q]; i++)
              if (elp[q][i]!=-1)
                elp[u+1][i+u-q] = alpha_to[(d[u]+nn-d[q]+elp[q][i])%nn] ;
            for (i=0; i<=l[u]; i++)
              { elp[u+1][i] ^= elp[u][i] ;
                elp[u][i] = index_of[elp[u][i]] ;  /*convert old elp value to index*/
              }
          }
        u_lu[u+1] = u-l[u+1] ;

/* form (u+1)th discrepancy */
        if (u<nn-kk)    /* no discrepancy computed on last iteration */
          {
            if (s[u+1]!=-1)
                   d[u+1] = alpha_to[s[u+1]] ;
            else
              d[u+1] = 0 ;
            for (i=1; i<=l[u+1]; i++)
              if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                d[u+1] ^= alpha_to[(s[u+1-i]+index_of[elp[u+1][i]])%nn] ;
            d[u+1] = index_of[d[u+1]] ;    /* put d[u+1] into index form */
          }
      } while ((u < nn-kk) && (l[u+1] <= tt)) ;


      u++ ;
      if (l[u] <= tt)         /* can correct error */
       {
#ifdef DECODER_DEBUG
    printf("---------------------------------------------------------------\n");
      printf("The Final Error Locator Polynomial Lambda(x) in polynomial-form is: \n");
      for (i=0;i < l[u];i++){
	  printf("%#x x^%d + ",elp[u][i],i);
          if (i && (i % 7 == 0)) printf("\n"); /* 8 coefficients per line */
      }
      printf("%#x x^%d",elp[u][i],i);
      printf("\n");
#endif

/* put elp into index form */
         for (i=0; i<=l[u]; i++)   elp[u][i] = index_of[elp[u][i]] ;

#ifdef DECODER_DEBUG
      printf("The Final Error Locator Polynomial Lambda(x) in index-form is: \n");
      for (i=0;i < l[u];i++){
          printf("%d x^%d + ",elp[u][i],i);
          if (i && (i % 7 == 0)) printf("\n"); /* 8 coefficients per line */
      }  
      printf("%d x^%d",elp[u][i],i);
      printf("\n");
    printf("---------------------------------------------------------------\n");
#endif


/* find roots of the error location polynomial */
         for (i=1; i <= l[u]; i++)
           reg[i] = elp[u][i] ;
         count = 0 ;
         for (i=1; i <= nn; i++)
          {  q = 1 ;
             for (j=1; j<=l[u]; j++)
              if (reg[j]!=-1)
                { reg[j] = (reg[j]+j)%nn ;
                  q ^= alpha_to[reg[j]] ;
                } ;
             if (!q)        /* store root and error location number indices */
              { root[count] = i;
                loc[count] = nn-i ;
                count++ ;
              };
          } ;
#ifndef NO_PRINT
	  printf("Computed Error Locations\n");
	  for (i=0;i < l[u];i++) printf("%d ",loc[i]);
	  printf("\n");
#endif
         if (count == l[u])    /* no. roots = degree of elp hence <= tt errors */
          {
#ifdef DECODER_DEBUG
 	printf("\n No of roots of Lambda(x) found by Chien search equals with \n");
 	printf("the degree of Lambda(x). Hence channel presumably caused <= %d errors\n",tt);
	printf("The roots found are (in index-form) : \n");
	for (i=0;i < count;i++)
	    printf("%d ",root[i]);
        printf("\n");
#endif 
/* form polynomial z(x) */
           for (i=1; i <= l[u]; i++)        /* Z[0] = 1 always - do not need */
            { if ((s[i]!=-1) && (elp[u][i]!=-1))
                 z[i] = alpha_to[s[i]] ^ alpha_to[elp[u][i]] ;
              else if ((s[i]!=-1) && (elp[u][i]==-1))
                      z[i] = alpha_to[s[i]] ;
                   else if ((s[i]==-1) && (elp[u][i]!=-1))
                          z[i] = alpha_to[elp[u][i]] ;
                        else
                          z[i] = 0 ;
              for (j=1; j<i; j++)
                if ((s[j]!=-1) && (elp[u][i-j]!=-1))
                   z[i] ^= alpha_to[(elp[u][i-j] + s[j])%nn] ;
              z[i] = index_of[z[i]] ;         /* put into index form */
            } ;

#ifdef DECODER_DEBUG
    printf("---------------------------------------------------------------\n");
      printf(" \n The Final Error Evaluator Polynomial Omega(x) in index-form is: \n");
      printf("0 x^0 + ");
      for (i=1;i <= l[u];i++){
          printf("%d x^%d + ",z[i],i);
          if (i && (i % 7 == 0)) printf("\n"); /* 8 coefficients per line */
      }
      printf("%d x^%d",z[i],i);
      printf("\n");
#endif

#ifdef DECODER_DEBUG
      printf(" \n The Final Error Evaluator Polynomial Omega(x) in polynomial-form is: \n");
      printf("%#x x^%d + ",1,0);
      for (i=1;i <= l[u];i++){
          printf("%#x x^%d + ",z[i],i);
          if (i && (i % 7 == 0)) printf("\n"); /* 8 coefficients per line */
      }  
      printf("%#x x^%d",z[i],i);
      printf("\n");
    printf("---------------------------------------------------------------\n");
#endif

  /* evaluate errors at locations given by error location numbers loc[i] */
           for (i=0; i<nn; i++)
             { err[i] = 0 ;
               if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                 recd[i] = alpha_to[recd[i]] ;
               else  recd[i] = 0 ;
             }
           for (i=0; i < l[u]; i++)    /* compute numerator of error term first */
            { err[loc[i]] = 1;       /* accounts for z[0] */
              for (j=1; j<=l[u]; j++){
                if (z[j]!=-1)
                  err[loc[i]] ^= alpha_to[(z[j]+j*root[i])%nn] ;
	      } /* z(x) evaluated at X(l)**(-1) */
              if (err[loc[i]] != 0) /* term X(l)**(1-b0) */
		 err[loc[i]] = alpha_to[(index_of[err[loc[i]]]+root[i]*(b0+nn-1))%nn];
              if (err[loc[i]]!=0)
               {
		 err[loc[i]] = index_of[err[loc[i]]] ;
                 q = 0 ;     /* form denominator of error term */
                 for (j=0; j<l[u]; j++)
                   if (j!=i)
                     q += index_of[1^alpha_to[(loc[j]+root[i])%nn]] ;
                 q = q % nn ;
                 err[loc[i]] = alpha_to[(err[loc[i]]-q+nn)%nn] ;
                 recd[loc[i]] ^= err[loc[i]] ;  /*recd[i] must be in polynomial form */
               }
            }
#ifdef DECODER_DEBUG
    		printf("---------------------------------------------------------------\n");
		printf("The computed Error Value Bytes are (in polynomial-form) :\n");
		printf("Location (Error) \n");
		for (i=0;i < l[u];i++){
		    printf("%d (%#x) ",loc[i],err[loc[i]]);
		}
		printf("\n");
    		printf("---------------------------------------------------------------\n");
#endif

	    return(1);
          }
          else{    /* no. roots != degree of elp => >tt errors and cannot solve */
#ifdef DECODER_DEBUG
		printf("No of roots of Lambda(x) found by Chien search\n");
		printf("does not equal the degree of the Lambda(x).\n");
		printf("Hence channel has definitely caused >= %d errors\n",tt);
#endif
           	for (i=0; i<nn; i++){        /* could return error flag if desired */
               	   if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                 	recd[i] = alpha_to[recd[i]] ;
                   else  
			recd[i] = 0 ;     /* just output received word as is */
		}
		return(2);
	  }
       }
      else{         /* elp has degree has degree >tt hence cannot solve */
#ifdef DECODER_DEBUG
		printf("Degree of the Lambda(x) > %d. \n",tt);
		printf("Hence channel has definitely caused >= %d errors\n",tt);
#endif
       	for (i=0; i<nn; i++){       /* could return error flag if desired */
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                recd[i] = alpha_to[recd[i]] ;
            else  
		recd[i] = 0 ;     /* just output received word as is */
	}
	return(3);
      }
#ifdef DECODER_DEBUG
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
   }
   else{       /* no non-zero syndromes => no errors: output received codeword */
#ifdef DECODER_DEBUG
    printf("Received vector is a codeword. Channel has presumably caused no error\n");
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
       	for (i=0; i < nn; i++){
       	   if (recd[i] != -1)        /* convert recd[] to polynomial form */
         	recd[i] = alpha_to[recd[i]] ;
           else  
	 	recd[i] = 0 ;
 	}
	return(0);
   }
 }

/* Performs ERRORS+ERASURES decoding of RS codes. If decoding is successful, writes
the codeword into recd[] itself. Otherwise recd[] is unaltered except in the
erased positions where it is set to zero. 
  First "no_eras" erasures are declared by the calling program. Then, the 
maximum # of errors correctable is t_after_eras = floor((2*tt-no_eras)/2). If the number
of channel errors is not greater than "t_after_eras" the transmitted codeword
 will be recovered. Details of algorithm can be found
in R. Blahut's "Theory ... of Error-Correcting Codes". */

int eras_dec_rs(eras_pos,no_eras)
int eras_pos[2*tt],no_eras;
{
  register int i,j,r,u,q,tmp,tmp1,tmp2,num1,num2,den,pres_root,pres_loc;
  int phi[2*tt+1],tmp_pol[2*tt+1]; /* The erasure locator in polynomial form */
  int U,syn_err,discr_r,deg_phi,deg_lambda,L,deg_omega,t_after_eras;
  int lambda[2*tt+1],s[2*tt+1],lambda_pr[2*tt+1];/* Err+Eras Locator poly and syndrome poly */ 
  int b[2*tt+1],T[2*tt+1],omega[2*tt+1];
  int syn_error=0, root[2*tt], z[tt+1], err[nn], reg[2*tt+1] ;
  int loc[2*tt],count = 0;

/* Maximum # ch errs correctable after "no_eras" erasures */
   t_after_eras = (int)floor((2.0*tt-no_eras)/2.0);

/* Compute erasure locator polynomial phi[x] */
  for (i=0;i < nn-kk+1;i++) tmp_pol[i] = 0;
  for (i=0;i < nn-kk+1;i++) phi[i] = 0; 
  if (no_eras > 0){
      phi[0] = 1; /* index form */
      phi[1] = alpha_to[eras_pos[0]];
      for (i=1;i < no_eras;i++){
          U = eras_pos[i];
	  for (j=1;j < i+2;j++){
	     tmp1 = index_of[phi[j-1]];
	     tmp_pol[j] = (tmp1 == -1)?  0 : alpha_to[(U+tmp1)%nn];
	  }
	  for (j=1;j < i+2;j++)
	      phi[j] = phi[j]^tmp_pol[j];
      }
  /* put phi[x] in index form */
     for (i=0;i < nn-kk+1;i++) phi[i] = index_of[phi[i]];
#ifdef ERASURE_DEBUG
/* find roots of the erasure location polynomial */
         for (i=1; i <= no_eras; i++) reg[i] = phi[i] ;
         count = 0 ;
         for (i=1; i <= nn; i++)
          {  q = 1 ;
             for (j=1; j <= no_eras; j++)
              if (reg[j]!=-1)
                { reg[j] = (reg[j]+j)%nn ;
                  q ^= alpha_to[(reg[j])%nn] ;
                } ;
             if (!q)        /* store root and error location number indices */
              { root[count] = i;
                loc[count] = nn-i ;
                count++ ;
              };
          } ;
	  if (count != no_eras){ printf("\n phi(x) is WRONG\n"); exit(1);}
#ifndef NO_PRINT 
          printf("\n Erasure positions as determined by roots of Eras Loc Poly:\n");
          for (i=0;i < count;i++) printf("%d ",loc[i]);
	  printf("\n");
#endif
#endif
  }

/* recd[] is in polynomial form, convert to index form */
   for (i=0;i < nn;i++) recd[i] = index_of[recd[i]];

/* first form the syndromes; i.e., evaluate recd(x) at roots of g(x) namely
 @**(b0+i), i = 0, ... ,(2*tt-1) */
   for (i=1; i <= nn-kk; i++)
    { s[i] = 0 ;
      for (j=0; j < nn; j++)
        if (recd[j] != -1)
          s[i] ^= alpha_to[(recd[j]+(b0+i-1)*j)%nn] ;      /* recd[j] in index form */
/* convert syndrome from polynomial form to index form  */
      if (s[i] != 0)  syn_error = 1 ;   /* set flag if non-zero syndrome => error */
      s[i] = index_of[s[i]] ;
    };


   if (syn_error){ /* if syndrome is zero, modified recd[] is a codeword */
/* Begin Berlekamp-Massey algorithm to determine error+erasure locator polynomial */
   	r = no_eras;
   	deg_phi = no_eras;
   	L = no_eras;
   	if (no_eras > 0){
        	/* Initialize lambda(x) and b(x) (in poly-form) to phi(x) */
        	for (i=0;i < deg_phi+1;i++) lambda[i] = (phi[i] == -1)? 0 : alpha_to[phi[i]];
        	for (i=deg_phi+1;i < 2*tt+1;i++) lambda[i] = 0;
        	deg_lambda = deg_phi;
        	for (i=0;i < 2*tt+1;i++) b[i] = lambda[i];
   	}
   	else{
		lambda[0] = 1;
		for (i=1;i < 2*tt+1;i++) lambda[i] = 0;
		for (i=0;i < 2*tt+1;i++) b[i] = lambda[i];
   	}
	while (++r <= 2*tt){ /* r is the step number */
	        /* Compute discrepancy at the r-th step in poly-form */
	        discr_r = 0;
	        for (i=0;i < 2*tt+1;i++){
		    if ((lambda[i] != 0) && (s[r-i] != -1)){
		       tmp = alpha_to[(index_of[lambda[i]]+s[r-i])%nn];
		       discr_r ^= tmp;
		    }
	        }
	    if (discr_r == 0){
	       /* 3 lines below: B(x) <-- x*B(x) */
	       tmp_pol[0] = 0;
	       for (i=1;i < 2*tt+1;i++) tmp_pol[i] = b[i-1];
	       for (i=0;i < 2*tt+1;i++) b[i] = tmp_pol[i];
	    }
	    else{
	       /* 5 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
	       T[0] = lambda[0];
	       for (i=1;i < 2*tt+1;i++){
		    tmp =  (b[i-1] == 0)? 0 : alpha_to[(index_of[discr_r]+index_of[b[i-1]])%nn];
		    T[i] = lambda[i]^tmp;
	       }

	       if (2*L <= r+no_eras-1){
		    L = r+no_eras-L;
		    /* 2 lines below: B(x) <-- inv(discr_r) * lambda(x) */
		    for (i=0;i < 2*tt+1;i++)
		        b[i] = (lambda[i] == 0)? 0 : alpha_to[(index_of[lambda[i]]-index_of[discr_r]+nn)%nn];
		    for (i=0;i < 2*tt+1;i++) lambda[i] = T[i]; 
	       }
	       else{
		    for (i=0;i < 2*tt+1;i++) lambda[i] = T[i]; 
	            /* 3 lines below: B(x) <-- x*B(x) */
	            tmp_pol[0] = 0;
	            for (i=1;i < 2*tt+1;i++) tmp_pol[i] = b[i-1];
	            for (i=0;i < 2*tt+1;i++) b[i] = tmp_pol[i];
	       }
	    }
          }
/* Put lambda(x) into index form */
    for (i=0; i < 2*tt+1; i++)
	lambda[i] = index_of[lambda[i]];
    
/* Compute deg(lambda(x)) */
    deg_lambda = 2*tt;
    while ((lambda[deg_lambda] == -1) && (deg_lambda > 0)) 
	--deg_lambda;
    if (deg_lambda <= 2*tt){
/* Find roots of the error+erasure locator polynomial. By Chien Search */
         for (i=1; i < 2*tt+1; i++) reg[i] = lambda[i] ;
         count = 0 ; /* Number of roots of lambda(x) */
         for (i=1; i <= nn; i++)
          {  q = 1 ;
             for (j=1; j <= deg_lambda; j++)
              if (reg[j] != -1)
                { reg[j] = (reg[j]+j)%nn ;
                  q ^= alpha_to[(reg[j])%nn] ;
                } ;
             if (!q)        /* store root (index-form) and error location number */
              { root[count] = i;
                loc[count] = nn-i;
                count++;
              };
          } ;
#ifndef NO_PRINT
	  printf("\n Final error positions:\t");
	  for (i=0;i < count;i++) printf("%d ",loc[i]);
	  printf("\n");
#endif 
	if (deg_lambda == count){ /* correctable error */
/* Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo x**(nn-kk)). in poly-form */
	     for (i=0;i < 2*tt;i++){
		omega[i] = 0;
		for (j=0;(j < deg_lambda+1) && (j < i+1);j++){
		    if ((s[i+1-j] != -1) && (lambda[j] != -1))
		       	tmp = alpha_to[(s[i+1-j]+lambda[j])%nn];
		    else
			tmp = 0;
		    omega[i] ^= tmp;	
		}
	     }
	     omega[2*tt] = 0;
/* Compute lambda_pr(x) = formal derivative of lambda(x) in poly-form */
	for (i=0;i < tt;i++){
	   lambda_pr[2*i+1] = 0;
	   lambda_pr[2*i] = (lambda[2*i+1] == -1)? 0 : alpha_to[lambda[2*i+1]];
	}
	lambda_pr[2*tt] = 0;
/* Compute deg(omega(x)) */
    deg_omega = 2*tt;
    while ((omega[deg_omega] == 0) && (deg_omega > 0)) 
	--deg_omega;
/* Compute error values in poly-form. num1 = omega(inv(X(l))), 
  num2 = inv(X(l))**(b0-1) and den = lambda_pr(inv(X(l))) all in poly-form */
        for (j=0;j < count;j++){
	   pres_root = root[j];
	   pres_loc = loc[j];
	   num1 = 0;
	   for (i=0;i < deg_omega+1;i++){
	      if (omega[i] != 0) 
		tmp = alpha_to[(index_of[omega[i]]+i*pres_root)%nn];
	      else
		tmp = 0;
	      num1 ^= tmp;
	   }
	   num2 = alpha_to[(pres_root*(b0-1))%nn];
	   den = 0;
	   for (i=0;i < deg_lambda+1;i++){
	      if (lambda_pr[i] != 0)
		tmp = alpha_to[(index_of[lambda_pr[i]]+i*pres_root)%nn];
	      else
		tmp = 0;
	      den ^= tmp;
	   }
	   if (den == 0){
		printf("\n ERROR: denominator = 0\n");
	   }
	   err[pres_loc] = 0;
	   if (num1 != 0){
		err[pres_loc] = alpha_to[(index_of[num1]+index_of[num2]+(nn-index_of[den]))%nn];
	   }
	 }
/* Correct word by subtracting out error bytes. First convert recd[] into poly-form */
		for (i=0;i < nn;i++) recd[i] = (recd[i] == -1)? 0 : alpha_to[recd[i]];
		for (j=0;j < count;j++)
		    recd[loc[j]] ^= err[loc[j]];
		return(1);
	}
	else /* deg(lambda) unequal to number of roots => uncorrectable error detected */
	     return(2);
     }
     else /* deg(lambda) > 2*tt => uncorrectable error detected */
	 return(3);
   }
   else{  /* no non-zero syndromes => no errors: output received codeword */
       for (i=0; i < nn; i++){
          if (recd[i] != -1)        /* convert recd[] to polynomial form */
           recd[i] = alpha_to[(recd[i])%nn];
          else
           recd[i] = 0;
       }
       return(0);
   }
}


/* 
    Generates the Galois field and then the generator polynomial.
  Must be recompiled for each different generator polynomial after 
  setting parameters at the top of this file. 
*/

main()
{
  int i,j,no_cws,nn_short,kk_short,byte,error_byte,no_ch_errs,iter,error_flag;
  int received[1000],cw[1000],decode_flag,no_decoder_errors;
  long seed;
  unsigned short int tmp;
  char cw_file[30],file_in[30];
  FILE *fp_out,*fp;
  /**** Storage for random seed ****/
  unsigned short int store[3];
  long int nrand48();
  double x,erand48(),prob_symb_err;
  void srand48();

 /* Generate the Galois field GF(2**mm) */
 printf("Generating the Galois field GF(%d)\n\n",n); 
 generate_gf();

#ifdef DEBUG
  printf("Look-up tables for GF(2**%2d)\n",mm) ;
  printf("i\t\t alpha_to[i]\t index_of[i]\n") ;
  for (i=0; i < n; i++)
   printf("%3d \t\t %#x \t\t %d\n",i,alpha_to[i],index_of[i]) ;
  printf("\n The polynomial-form representation of @**i, 0 <= i < %d \n",n);
  printf("is obtained as follows: It is given by the binary representation of \n");
  printf("the integer alpha_to[i] with LS-Bit corresponding to @**0\n");
  printf("and the next bit to @**1 and so on.\n");
  printf("\n The index of a Galois field element x  whose polynomial-form \n");
  printf("representation is the binary representation of integer \n");
  printf("j, 0 <= j < %d is given by the integer index_of[j]\n",n);
  printf("Therefore, @**j = x \n");
  printf("\n EXAMPLES IN GF(256) \n");
  printf("The polynomial-form of @^8 is the binary representation of the integer alpha_to[8] viz. 0x1d. i.e., @**8 = @**4 + @**3 + @**2 + @**0 \n");
  printf("The index of x whose polynomial-form is 0x72 (integer 155) is index_of[155] = 217, so @**217 = x \n");
	 
#endif
 
 
  printf("Enter b0 = index of the first root \n");
  scanf("%d",&b0);

  gen_poly();

  printf("\n g(x) of the %d-error correcting RS (%d,%d) code: \n",tt,nn,kk);
  printf("Coefficient is the exponent of @ = primitive element\n");
  for (i=0;i <= nn-kk ;i++){
      printf("%d x^%d ",gg[i],i);
      if (i < nn-kk) printf(" + ");
      if (i && (i % 7 == 0)) printf("\n"); /* 8 coefficients per line */
  }
  printf("\n");

  printf("\n Coefficient is the representation in basis {@^%d,...,@,1}\n",mm-1);
  for (i=0;i <= nn-kk ;i++){
      printf("%#x x^%d ",alpha_to[gg[i]],i);
      if (i < nn-kk) printf(" + ");
      if (i && (i % 7 == 0)) printf("\n"); /* 8 coefficients per line */
  }
  printf("\n\n");

#ifdef ENCODE
  printf("Enter length of the shortened code\n");
  scanf("%d",&nn_short);
  if ((nn_short < 2*tt) || (nn_short > nn)){
     printf("Invalid entry %d for shortened length\n",nn_short);
     exit();
  }
  kk_short = kk - (nn-nn_short); /* compute dimension of shortened code */
  printf("The (%d,%d) %d-error correcting RS code\n",nn_short,kk_short,tt);
  printf("Enter the number of codewords desired \n");
  scanf("%d",&no_cws);
  printf("Enter the filename for storing cws\n");
  scanf("%s",cw_file);
  printf("Enter 3 positive integers as seed \n");
  for (i=0;i < 3;i++){ scanf("%hu",&tmp);store[i] = tmp;}
  
  if ((fp_out = fopen(cw_file,"w")) == NULL){
 	printf("Could not open %s\n",cw_file);
	exit();
  }
  
  /**** BEGIN: Encoding random information vectors ****/
     /* Pad with zeros rightmost (kk-kk_short) positions */
  for (j=kk_short;j < kk;j++) data[j] = 0;
  for (i=0;i < no_cws;i++){
      for (j=0;j < (int)kk_short;j++) data[j] = (int) (nrand48(store) % 256);

      encode_rs();

      for (j=0;j < (int)nn_short-(int)kk_short;j++) /* parity bytes first */
	 fprintf(fp_out,"%02x ",bb[j]);
      for (j=0;j < (int)kk_short;j++) /* info bytes last */
	 fprintf(fp_out,"%02x ",data[j]);
      fprintf(fp_out,"\n");
  }
  fclose(fp_out);
#endif

#ifdef DECODE
  printf("Enter length of the shortened code\n");
  scanf("%d",&nn_short);
  if ((nn_short < 2*tt) || (nn_short > nn)){
     printf("Invalid entry %d for shortened length\n",nn_short);
     exit();
  }
  kk_short = kk - (nn-nn_short); /* compute dimension of shortened code */
  printf("The (%d,%d) %d-error correcting RS code\n",nn_short,kk_short,tt);
  printf("Enter a random positive integer\n");
  scanf("%ld",&seed);
  srand48(seed);

  printf("Enter filename from which to read codewords\n"); scanf("%s",file_in);
  printf("Enter no of codewords \n"); scanf("%d",&no_cws);

  if ((fp = fopen(file_in,"r")) == NULL){             
        printf("Could not open %s\n",file_in);
        exit();
  }

  prob_symb_err = (double)tt/ (8.0 * (double)nn_short);
  /* CHANGE above probablity of symbol error if desired */
  /* printf("Probability of symbol error = %6e\n",prob_symb_err);*/

  no_decoder_errors = 0;
  iter = -1;
  while (++iter < no_cws){
     	/* read in codewords from file */
  	for (i=0;i < nn_short;i++){
      	    fscanf(fp,"%02x ",&byte);
      	    cw[i] = byte;
  	}
  	/* printf("The Transmitted Codeword\n");
  	for (i=0;i < nn_short;i++) printf("%02x ",cw[i]); printf("\n");*/

	
#ifndef NO_PRINT
	 printf("\n\n\n\n Transmitting codeword %d \n",iter);
	 printf("Channel caused following errors Location (Error): \n");
#endif
	 no_ch_errs = 0;
	 for (i=0;i < nn_short;i++){
	      x = erand48(store);
	      if (x < prob_symb_err){
		 error_byte = (int) (lrand48() % 256);
		 received[i] = cw[i] ^ error_byte;
		 ++no_ch_errs;
#ifndef NO_PRINT
		 printf("%d (%#x) ",i,error_byte);
#endif
	      }
	      else
		 received[i] = cw[i];
	  }
#ifndef NO_PRINT
	  printf("\n");
	  printf("Channel caused %d errors\n",no_ch_errs);
	  /* for (i=0;i < nn_short;i++) printf("%02x ",received[i]); printf("\n");*/
#endif

	  /* Pad with zeros and decode as if in (255,kk) tt-error correcting code */
          for (i=0;i < nn-kk;i++)  recd[i] = received[i]; /* parity bytes */
	  for (i=nn-kk;i < nn-kk+kk_short;i++) recd[i] = received[i]; /* info bytes */ 
	  for (i=nn-kk+kk_short;i < nn;i++) recd[i] = 0; /* zero padding */
	
	  decode_flag = decode_rs();
#ifndef NO_PRINT
	  printf("decode_rs() returned %d\n",decode_flag);
          error_flag = 0;
          for (i=0; i < kk_short; i++){ 
	      if (recd[i+nn-kk] != cw[i+nn-kk]){
		 error_flag = 1;
		 break;
	      }
	  }
#endif

	  if (decode_flag == 2){
             if (no_ch_errs <= max_errs){
                 printf("%d ch errs  <=  max # correctable errs but \n",no_ch_errs,max_errs);
                 printf("\n DECODER ERROR: deg(lambda) unequal to #roots \n");
                 exit(2); /* DECODER ERROR condition */
             }
	  }
	  else if (decode_flag == 3){
	     if (no_ch_errs <= max_errs){
                printf(" %d ch errs  <= max # correctable errs but \n",no_ch_errs,max_errs);
       		printf("\n deg(lambda) > 2*tt \n");
             	exit(3);/* DECODER ERROR condition */
             }
	  }
	  else{
	     if ((no_ch_errs <= max_errs) && (error_flag == 1)){
		printf("DECODER FAILED TO CORRECT %d ERRORS\n",no_ch_errs);
		++no_decoder_errors;
#ifndef NO_PRINT
		printf("The Transmitted Codeword\n");
		for (i=0;i < nn_short;i++) printf("%02x ",cw[i]); printf("\n");
#endif 
	     }
	  }

   } /* closing iteration loop */
   if (no_decoder_errors > 0)
      printf("The number of decoder errors were %d\n",no_decoder_errors);
   else
      printf("Decoder corrected all occurrences of %d or less errors\n",tt);
#endif

}

