#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "../fewbody.h"

int calc_units(fb_hier_t hier, fb_units_t *units)
{
  /* mass unit is total system mass */
  units->m = hier.hier[hier.hi[1]+0].m + hier.hier[hier.hi[1]+1].m +	\
    hier.hier[hier.hi[1]+2].m + hier.hier[hier.hi[1]+3].m +		\
    hier.hier[hier.hi[1]+4].m;
  /* length unit is one AU */
  units->l = FB_CONST_AU;
  /* everything else is derived */
  units->E = FB_CONST_G * fb_sqr(units->m) / units->l;
  units->v = sqrt(units->E/units->m);
  units->t = units->l / units->v;
  return(0);
}

int main(int argc, char *argv[])
{
  int j;
  double t;
  fb_hier_t hier;
  fb_input_t input;
  fb_ret_t retval;
  fb_units_t units;
  char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
  
  /* set input parameters */
  input.ks = 0; /* turn K-S regularization off */
  input.tstop = 1.0e4; /* stopping time in units of units.t */
  input.Dflag = 0; /* don't output dynamical info to stdout */
  input.dt = 0.0; /* irrelevant when Dflag=0 */
  input.tcpustop = 120.0; /* stopping CPU time in seconds */
  input.absacc = 1.0e-9; /* integrator absolute accuracy */
  input.relacc = 1.0e-9; /* integrator relative accuracy */
  input.ncount = 500; /* number of integration steps between calls to fb_classify() */
  input.tidaltol = 1.0e-6; /* tidal perturbation required to force numerical */
                           /* integration of a binary node */
  input.fexp = 1.0; /* radius expansion factor of merger products */
  fb_debug = 0;
  
  /* initialize a few things */
  t = 0.0;
  hier.nstarinit = 5;
  hier.nstar = 5;
  fb_malloc_hier(&hier);
  fb_init_hier(&hier);
  
  /* set stellar properties */
  for (j=0; j<hier.nstar; j++) {
    hier.hier[hier.hi[1]+j].ncoll = 1;
    hier.hier[hier.hi[1]+j].id[0] = j;
    snprintf(hier.hier[hier.hi[1]+j].idstring, FB_MAX_STRING_LENGTH, "%d", j);
    hier.hier[hier.hi[1]+j].n = 1;
    hier.hier[hier.hi[1]+j].obj[0] = NULL;
    hier.hier[hier.hi[1]+j].obj[1] = NULL;
    hier.hier[hier.hi[1]+j].Eint = 0.0;
    hier.hier[hier.hi[1]+j].Lint[0] = 0.0;
    hier.hier[hier.hi[1]+j].Lint[1] = 0.0;
    hier.hier[hier.hi[1]+j].Lint[2] = 0.0;
    hier.hier[hier.hi[1]+j].R = 1.0e-3 * FB_CONST_RSUN;
    hier.hier[hier.hi[1]+j].m = FB_CONST_MSUN;
    hier.hier[hier.hi[1]+j].x[2] = 0.0; /* everything in x-y plane */
    hier.hier[hier.hi[1]+j].v[2] = 0.0; 
  }
  
  /* set some positions and velocities */
  hier.hier[hier.hi[1]+0].x[0] = 10.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+0].x[1] = 0.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+0].v[0] = 10.0e5;
  hier.hier[hier.hi[1]+0].v[1] = 0.0;
  
  hier.hier[hier.hi[1]+1].x[0] = 10.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+1].x[1] = 10.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+0].v[0] = -10.0e5;
  hier.hier[hier.hi[1]+0].v[1] = 0.0;
  
  hier.hier[hier.hi[1]+2].x[0] = -30.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+2].x[1] = -5.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+0].v[0] = 0.0;
  hier.hier[hier.hi[1]+0].v[1] = 5.0e5;
  
  hier.hier[hier.hi[1]+3].x[0] = 0.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+3].x[1] = 27.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+0].v[0] = -5.0e5;
  hier.hier[hier.hi[1]+0].v[1] = -4.0e5;
  
  hier.hier[hier.hi[1]+4].x[0] = 10.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+4].x[1] = -32.0 * FB_CONST_AU;
  hier.hier[hier.hi[1]+0].v[0] = 5.0e5;
  hier.hier[hier.hi[1]+0].v[1] = -1.0e5;
  
  /* set units with our own routine, then normalize */
  calc_units(hier, &units);
  fb_normalize(&hier, units);
  
  /* call Fewbody to evolve system */
  retval = fewbody(input, &hier, &t);
  
  /* all the rest is parsing the output */
  if (retval.retval == 1) {
    fprintf(stderr, "encounter complete.\n");
  } else {
    fprintf(stderr, "encounter NOT complete.\n");
  }
  
  fprintf(stderr, "final configuration:  %s  (%s)\n",
	  fb_sprint_hier(hier, string1),
	  fb_sprint_hier_hr(hier, string2));
  
  fprintf(stderr, "t_final=%.6g (%.6g yr)  t_cpu=%.6g s\n",	\
	  t, t*units.t/FB_CONST_YR, retval.tcpu);
  
  fprintf(stderr, "DeltaL/L0=%.6g  DeltaL=%.6g\n", retval.DeltaLfrac, retval.DeltaL);
  fprintf(stderr, "DeltaE/E0=%.6g  DeltaE=%.6g\n", retval.DeltaEfrac, retval.DeltaE);
  fprintf(stderr, "Rmin=%.6g (%.6g RSUN)  Rmin_i=%d  Rmin_j=%d\n",	\
	  retval.Rmin, retval.Rmin*units.l/FB_CONST_RSUN, retval.Rmin_i, retval.Rmin_j);
  fprintf(stderr, "Nosc=%d (%s)\n", retval.Nosc, 
	  (retval.Nosc>=1?"resonance":"non-resonance"));

  fprintf(stderr, "orbital parameters of outermost binaries:\n");
  for (j=0; j<hier.nobj; j++) {
    if (hier.obj[j]->n >= 2) {
      fprintf(stderr, "j=%d  a=%.6g AU  e=%.6g\n", j, 
	      hier.obj[j]->a * units.l / FB_CONST_AU, hier.obj[j]->e);
    }
  }

  fb_free_hier(hier);
  return(0);
}
