#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


#pragma once

#define BM_ERROR 1
#define BM_SUCCESS 0

typedef struct bm_handle {
  uint32_t num_iters;
  uint32_t iter;

  double *event_lengths;

  double event_start;
} bm_handle;

#ifdef NO_BENCHMARK

inline int bm_init(bm_handle *bm, uint32_t num_iters) {}

inline void bm_start(bm_handle *bm) {}
inline void bm_resume(bm_handle *bm) {}
inline void bm_pause(bm_handle *bm) {}
inline void bm_stop(bm_handle *bm) {}

inline void bm_print_events(bm_handle *bm, FILE *ptr) {}
inline void bm_destroy(bm_handle *bm) {}



inline void bm_resume(bm_handle *bm){};
inline void bm_pause(bm_handle *bm){};
inline void bm_print_events(bm_handle *bm, FILE *ptr){}
inline void bm_destroy(bm_handle *bm){}


#else 

static inline double rtclock() {
  struct timeval Tp;
  gettimeofday(&Tp, NULL);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}


int bm_init(bm_handle *bm, uint32_t num_iters);


static inline void bm_start(bm_handle *bm) { bm->event_start = rtclock(); }

static inline void bm_stop(bm_handle *bm) {
  double diff = rtclock() - bm->event_start;

  bm->event_lengths[bm->iter] = diff;
  bm->iter++;
}

static inline void bm_resume(bm_handle *bm) { bm->event_start = rtclock(); }

static inline void bm_pause(bm_handle *bm) {
  double diff = rtclock() - bm->event_start;
  bm->event_lengths[bm->iter] += diff;
}

void bm_print_events(bm_handle *bm, FILE *ptr);
void bm_destroy(bm_handle *bm);

#endif
