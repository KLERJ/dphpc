#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>


#pragma once


#define BM_ERROR 1
#define BM_SUCCESS 0



typedef struct bm_handle {
    uint32_t num_iters;
    uint32_t iter;
    

    double * event_lengths;

    double event_start;
} bm_handle;

#ifdef NO_BENCHMARK

inline int bm_init(bm_handle *bm, uint32_t num_iters){}

inline void bm_start(bm_handle *bm){}

inline void bm_stop(bm_handle *bm){}

inline void bm_print_events(bm_handle *bm, FILE *ptr){}
inline void bm_destroy(bm_handle *bm){}


#else 
int bm_init(bm_handle *bm, uint32_t num_iters);

void bm_start(bm_handle *bm);

void bm_stop(bm_handle *bm);

void bm_print_events(bm_handle *bm, FILE *ptr);
void bm_destroy(bm_handle *bm);


#endif