#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rs.h"


void
fill_eras(int eras_pos[],int n)
{
	int i,j,t,work[NN];

	for(i=0;i<NN;i++)
		work[i] = i;
	for(j=NN-1;j>0;j--){
		i = random() % j;	/* not really uniform, I know */
		t = work[i];
		work[i] = work[j];
		work[j] = t;
	}
#ifdef notdef
	for(i=0;i<NN;i++)
	  printf(" %d",work[i]);
	printf("\n");
#endif
	for(i=0;i<n;i++)
		eras_pos[i] = work[i];
}

/* Return non-zero random number in range 0 - NN (NN power of 2 minus 1) */
int
randomnz(void)
{
	int i;

	while((i = random() & NN) == 0)
		;
	return i;
}

int
main(int argc,char *argv[])
{

	dtype data[NN];
	dtype tdata[NN];
	dtype ddata[NN];
	int eras_pos[NN];
	int i,trials,k;
	long t;
	int nerrors,nerase,ntrials,verbose,timetest;
	int detfails,fails;
	extern char *optarg;

	nerrors = nerase = 0;
	timetest = verbose = 0;
	ntrials = 100;
	while((i = getopt(argc,argv,"e:E:n:vt")) != EOF){
		switch(i){
		case 'e':	/* Number of errors per block */
			nerrors = atoi(optarg);
			break;
		case 'E':	/* Number of erasures per block */
			nerase = atoi(optarg);
			break;
		case 'n':	/* Number of trials */
			ntrials = atoi(optarg);
			break;
		case 'v':	/* Be verbose */
			verbose = 1;
			break;
		case 't':	/* Repeatedly decode the same block */
			timetest = 1;
			break;
		default:
			printf("usage: %s [-v] [-t] [-e errors] [-E erasures] [-n ntrials]\n",argv[0]);
			exit(1);
		}
	}
	printf("Reed-Solomon code is (%d,%d) over GF(%d)\n",
	   NN,KK,NN+1);
	printf("test erasures: %d errors %d\n",nerase,nerrors);
	if(2*nerrors + nerase > NN-KK)
		printf("Warning: %d errors and %d erasures exceeds the correction ability of the code\n",nerrors,nerase);
	time(&t);
	srandom(t);
	init_rs();

	if(timetest){
		printf("Speed timing test (repeated decoding of same block)\n");
		for(i=0;i<KK;i++)
			data[i] = random() & NN;
		encode_rs(data,&data[KK]);
		fill_eras(eras_pos,nerase+nerrors);
		if(verbose && nerase){
			printf("erasing:");
			for(i=0;i<nerase;i++)
				printf(" %d",eras_pos[i]);
			printf("\n");
		}
		if(verbose && nerrors){
			printf("erroring:");
			for(i=nerase;i<nerase+nerrors;i++)
				printf(" %d",eras_pos[i]);
			printf("\n");
		}
		memcpy(tdata,data,sizeof(data));
		for(i=0;i<nerase+nerrors;i++)
			tdata[eras_pos[i]] ^= randomnz();

		for(k=0;k<ntrials;k++){
			memcpy(ddata,tdata,sizeof(tdata));
			eras_dec_rs(ddata,eras_pos,nerase);
		}
		exit(0);
	}
	fails = detfails = 0;
	for(trials=0;trials < ntrials;trials++){
		if(verbose)
			printf("Trial %d:",trials);
		for(i=0;i<KK;i++)
			data[i] = random() & NN;
		encode_rs(data,&data[KK]);
		fill_eras(eras_pos,nerase+nerrors);
		if(verbose && nerase){
			printf(" erasing:");
			for(i=0;i<nerase;i++)
				printf(" %d",eras_pos[i]);
		}
		if(verbose && nerrors){
			printf(" erroring:");
			for(i=nerase;i<nerase+nerrors;i++)
				printf(" %d",eras_pos[i]);
		}
		if(verbose)
			printf("\n");
		memcpy(ddata,data,sizeof(data));
		for(i=0;i<nerase+nerrors;i++)
			ddata[eras_pos[i]] ^= randomnz();

		i = eras_dec_rs(ddata,eras_pos,nerase);
		if(verbose){
			printf("errs + erasures corrected: %d\n",i);
		}
		if(i == -1){
			detfails++;
			printf("RS decoder detected failure\n");
		} else if(memcmp(ddata,data,NN) != 0){
			fails++;
			printf("Undetected decoding failure!\n");
		}
	}
	printf("Trials: %d decoding failures: %d; not detected by decoder: %d\n",
	  ntrials,detfails,fails);
	return 0;
}
