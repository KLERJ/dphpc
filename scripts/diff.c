#include <math.h>
#include <stdio.h>

#define FEPS 1e-10

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: ./diff bin1, bin2\n");
    return 1;
  }

  FILE *f1 = fopen(argv[1], "rb");
  FILE *f2 = fopen(argv[2], "rb");

  if (f1 == NULL) {
    fprintf(stderr, "Couldn't read %s\n", argv[1]);
    return 1;
  }

  if (f2 == NULL) {
    fprintf(stderr, "Couldn't read %s\n", argv[2]);
    return 1;
  }

  int error = 0;

  double num_f1, num_f2;
  for (int i = 0;; i++) {
    int len1 = fread(&num_f1, sizeof(double), 1, f1);
    int len2 = fread(&num_f2, sizeof(double), 1, f2);

    if (len1 != len2) {
      fprintf(stderr, "Arrays don't have same length at %d\n", i);
      return 1;
    }

    if (FEPS < fabs(num_f1 - num_f2)) {
      printf("%i: %f != %f\n", i, num_f1, num_f2);
      error++;
    }

    if (len1 == 0) {
      break;
    }
  }

  if (error > 0) {
    printf("Found %d missmatches\n");
    return 1;
  }

  return 0;
}
