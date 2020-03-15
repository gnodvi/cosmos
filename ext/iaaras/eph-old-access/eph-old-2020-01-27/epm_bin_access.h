#ifndef __epm_bin_access_h__
#define __epm_bin_access_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

typedef struct
{
  int    jd;            // initial Julian date (integer part)
  double jd_frac;       // initial Julian date (fractional part)
  double subinterval;
  int num_coefficients;
  int dimension;
  int num_subintervals;
  FILE *file;
} EpmFile;

EpmFile * epm_file_open (const char *filename);
void epm_file_close (EpmFile *epmf);

int epm_calc_raw_values (EpmFile *epmf, int jd, double jd_frac, double *pos, double *vel);

#ifdef __cplusplus
}
#endif

#endif

