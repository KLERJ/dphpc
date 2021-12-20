
#include "bm.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

static double rtclock() {
  struct timeval Tp;
  int stat;
  stat = gettimeofday(&Tp, NULL);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

#ifdef NO_BENCHMARK

#else

int bm_init(bm_handle *bm, uint32_t num_iters) {
  bm->iter = 0;
  bm->num_iters = num_iters;

  bm->event_lengths = calloc(sizeof(double), num_iters);

  if (bm->event_lengths == NULL) {
    return BM_ERROR;
  }

  return BM_SUCCESS;
}

void bm_start(bm_handle *bm) { bm->event_start = rtclock(); }

void bm_resume(bm_handle *bm) { bm->event_start = rtclock(); }

void bm_pause(bm_handle *bm) {
  double diff = rtclock() - bm->event_start;
  bm->event_lengths[bm->iter] += diff;
}

void bm_stop(bm_handle *bm) {
  double diff = rtclock() - bm->event_start;

  bm->event_lengths[bm->iter] = diff;
  bm->iter++;
}

void bm_print_events(bm_handle *bm, FILE *ptr) {
  uint32_t i;
  for (i = 0; i < bm->iter - 1; i++) {
    fprintf(ptr, "%f,", bm->event_lengths[i]);
  }
  fprintf(ptr, "%f\n",
          bm->event_lengths[i]); // Print last element with new line
}

void bm_destroy(bm_handle *bm) { free(bm->event_lengths); }

#endif
