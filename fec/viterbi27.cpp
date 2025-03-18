// K=7 r=1/2 Viterbi decoder
//
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "fec.h"

void *create_viterbi27(int len) {
    find_cpu_mode();

    if (Cpu_mode == SSE2) {
        return create_viterbi27_sse2(len);
    } else {
        return create_viterbi27_port(len);
    }
}

void set_viterbi27_polynomial(int polys[2]) {
    if (Cpu_mode == SSE2) {
        set_viterbi27_polynomial_sse2(polys);
    } else {
        set_viterbi27_polynomial_port(polys);
    }
}

int init_viterbi27(void *p, int starting_state) {
    if (Cpu_mode == SSE2) {
        return init_viterbi27_sse2(p, starting_state);
    } else {
        return init_viterbi27_port(p, starting_state);
    }
}

int chainback_viterbi27(void *p, unsigned char *data, unsigned int nbits, unsigned int endstate) {
    if (Cpu_mode == SSE2) {
        return chainback_viterbi27_sse2(p, data, nbits, endstate);
    } else {
        return chainback_viterbi27_port(p, data, nbits, endstate);
    }
}

void delete_viterbi27(void *p) {
    if (Cpu_mode == SSE2) {
        delete_viterbi27_sse2(p);
    } else {
        delete_viterbi27_port(p);
    }
}

// nbits is the number of decoded data bits.
//
int update_viterbi27_blk(void *p, unsigned char syms[], int nbits) {
    if (p == nullptr)
        return -1;

    if (Cpu_mode == SSE2) {
        update_viterbi27_blk_sse2(p, syms, nbits);
    } else {
        update_viterbi27_blk_port(p, syms, nbits);
    }
    return 0;
}
