/* User include file for libfec
 */

#pragma once
#ifndef _FEC_H_
#define _FEC_H_

/* r=1/2 k=7 convolutional encoder polynomials
 * The NASA-DSN convention is to use V27POLYA inverted, then V27POLYB
 * The CCSDS/NASA-GSFC convention is to use V27POLYB, then V27POLYA inverted
 */
#define	V27POLYA	0x6d
#define	V27POLYB	0x4f

void *create_viterbi27(int len);
void set_viterbi27_polynomial(int polys[2]);
int init_viterbi27(void *vp,int starting_state);
int update_viterbi27_blk(void *vp,unsigned char sym[],int npairs);
int chainback_viterbi27(void *vp, unsigned char *data,unsigned int nbits,unsigned int endstate);
void delete_viterbi27(void *vp);

void *create_viterbi27_sse2(int len);
void set_viterbi27_polynomial_sse2(int polys[2]);
int init_viterbi27_sse2(void *p,int starting_state);
int chainback_viterbi27_sse2(void *p,unsigned char *data,unsigned int nbits,unsigned int endstate);
void delete_viterbi27_sse2(void *p);

void *create_viterbi27_port(int len);
void set_viterbi27_polynomial_port(int polys[2]);
int init_viterbi27_port(void *p,int starting_state);
int chainback_viterbi27_port(void *p,unsigned char *data,unsigned int nbits,unsigned int endstate);
void delete_viterbi27_port(void *p);
int update_viterbi27_blk_port(void *p,unsigned char *syms,int nbits);

/* Tables to map from conventional->dual (Taltab) and
 * dual->conventional (Tal1tab) bases
 */
extern unsigned char Taltab[],Tal1tab[];


/* CPU SIMD instruction set available */
extern enum cpu_mode {UNKNOWN=0,PORT,MMX,SSE,SSE2,ALTIVEC} Cpu_mode;
void find_cpu_mode(void); /* Call this once at startup to set Cpu_mode */

/* Determine parity of argument: 1 = odd, 0 = even */
#ifdef __i386__
static inline int parityb(unsigned char x){
  __asm__ __volatile__ ("test %1,%1;setpo %0" : "=g"(x) : "r" (x));
  return x;
}
#else
void partab_init();

static inline int parityb(unsigned char x){
  extern unsigned char Partab[256];
  extern int P_init;
  if(!P_init){
    partab_init();
  }
  return Partab[x];
}
#endif

static inline int parity(int x){
  /* Fold down to one byte */
  x ^= (x >> 16);
  x ^= (x >> 8);
  return parityb(x);
}

/* Useful utilities for simulation */
double normal_rand(double mean, double std_dev);
unsigned char addnoise(int sym,double amp,double gain,double offset,int clip);

extern int Bitcnt[];

/* Low-level data structures and routines */

extern "C" {
    int cpu_features(void);
    //int update_viterbi27_blk_sse2(void *p,unsigned char *syms,int nbits);
}

#endif // _FEC_H_
