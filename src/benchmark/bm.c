
#include "bm.h"
#include <stdio.h>
#include <stdlib.h>

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
