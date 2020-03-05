#include "epm_bin_access.h"

#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#define MAX_NCOEF 20

#ifdef _WIN32
#define snprintf _snprintf
#define strcasecmp _stricmp
#endif

EpmFile * epm_file_open (const char *filename)
{
  EpmFile *epmf = (EpmFile *)malloc(sizeof(EpmFile));

  if (epmf == NULL)
  {
    fprintf(stderr, "no memory?\n");
    return NULL;
  }
  
  epmf->file = fopen(filename, "rb");

  if (epmf->file == NULL)
  {
    fprintf(stderr, "error opening %s: %s\n", filename, strerror(errno));
    free(epmf);
    return NULL;
  }

#define READ(field) \
  do { \
    if (fread(&epmf->field, sizeof(epmf->field), 1, epmf->file) < 1)            \
    {                                                                           \
      fprintf(stderr, "error reading field %s from file %s", #field, filename); \
      fclose(epmf->file);                                                       \
      free(epmf);                                                               \
      return NULL;                                                              \
    } } while (0)

  READ(jd);
  READ(jd_frac);
  READ(subinterval);
  READ(num_coefficients);
  READ(dimension);
  READ(num_subintervals);
#undef READ

  if (epmf->dimension != 3 && epmf->dimension != 1)
  {
    fprintf(stderr, "EPM dimension is %d, can accept only 1 or 3\n", epmf->dimension);
    fclose(epmf->file);
    free(epmf);
    return NULL;
  }

  if (epmf->num_coefficients > MAX_NCOEF)
  {
    fprintf(stderr, "number of coefficients %d is too large\n", epmf->num_coefficients);
    fclose(epmf->file);
    free(epmf);
    return NULL;
  }

  return epmf;
}

void epm_file_close (EpmFile *epmf)
{
  fclose(epmf->file);
  free(epmf);
}

double calc_pos (int ncoef, double *antiderivative,
                 double *coef, double delta)
{
  int i;
  double x = 0.0;

  for (i = ncoef - 1; i >= 0; i--)
    x = x + coef[i] * antiderivative[i]; 

  x = 0.5 * delta * x + coef[ncoef];
  return x;
}

double calc_vel (int ncoef, double *polynomials,  double *coef)
{
  double v = 0.0;
  int i;

  for (i = ncoef - 1; i >= 0; i--)
    v += coef[i] * polynomials[i];
  return v;
}


// t is in the interval [-1 .. +1]
void calc_cheb_pol (int ncoef, double t, double *polynomials, double *antiderivatives)
{       
  double d;
  int i, j;
  int flag;
  
  polynomials[0] = 1.0;	
  polynomials[1] = t; 

  for (i = 2; i <= ncoef; i++) 
    polynomials[i] = 2.0 * polynomials[i - 1] * t - polynomials[i - 2];
  
  antiderivatives[0] = t;
  antiderivatives[1] = (polynomials[2] + polynomials[0]) * 0.25;
  for (i = 2; i < ncoef; i++) 
    antiderivatives[i] = 0.5 * (polynomials[i + 1] / (i + 1) -
                                polynomials[i - 1] / (i - 1));
  
  // Eliminate constant terms
  flag = 0; 
  d = 2.0 * t; 

  for (i = 3, j = 1; i < ncoef; i += 2, j++)
  { 
    d = 0.25 / j + 0.25 / (j + 1); 
    flag = 1 - flag; 
    if (flag)
      d = -d; 
    antiderivatives[i] += d;
  }
}

int epm_calc_raw_values (EpmFile *epmf, int jd, double jd_frac, double *pos, double *vel)
{
  double diff_int;
  double diff;
  int i_subinterval;
  double t;
  double coefs[MAX_NCOEF];
  double polynomials[MAX_NCOEF + 1];
  double antiderivatives[MAX_NCOEF];
  int i;

  if (epmf == NULL)
  {
    fprintf(stderr, "null EpmFile given\n");
    return 0;
  }

  diff_int = jd - epmf->jd;
  diff = (diff_int + (jd_frac - epmf->jd_frac)) / epmf->subinterval;

  // -----------------------------------
  i_subinterval = (int)floor(diff);
  // -----------------------------------
  
  t = 2.0 * (diff_int - i_subinterval * epmf->subinterval + (jd_frac - epmf->jd_frac))
          / epmf->subinterval - 1.0;

  if (i_subinterval < 0 || i_subinterval >= epmf->num_subintervals)
  {
    fprintf(stderr, "date %d+%lf not within the range of the theory\n", jd, jd_frac);
    return 0;
  }

  calc_cheb_pol(epmf->num_coefficients - 1, t, polynomials, antiderivatives);

  for (i = 0; i < epmf->dimension; i++)
  {
    if (fseek(epmf->file, 32 + i_subinterval * epmf->dimension * 8 *
              epmf->num_coefficients + i * 8 * epmf->num_coefficients, SEEK_SET) == -1)
    {
      fprintf(stderr, "error seeking: %s\n", strerror(errno));
      return 0;
    }
    
    if (fread(coefs, 8, epmf->num_coefficients, epmf->file) < epmf->num_coefficients)
    {
      fprintf(stderr, "error reading coefficients: %s\n", strerror(errno));
      return 0;
    }
    
    pos[i] = calc_pos(epmf->num_coefficients - 1, antiderivatives, coefs, epmf->subinterval);
    vel[i] = calc_vel(epmf->num_coefficients - 1, polynomials, coefs);
  }
  return 1;
}
