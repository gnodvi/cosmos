/* -*- linux-c -*- */
/* binsingle.c

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "fewbody.h"
#include "binsingle_simple.h"

/* calculate the units used */
int calc_units(fb_obj_t *single, fb_obj_t *binary, fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(single->m + binary->m)/(single->m * binary->m) * \
			(binary->obj[0]->m * binary->obj[1]->m / binary->a));
	units->l = binary->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
	
	return(0);
}

/* the main attraction */
int main(int argc, char *argv[])
{
	int i, n;
	unsigned long int seed;
	double m0, m10, m11, r0, r10, r11, a1, e1;
	double rtid, vinf, b, m1, M, mu, Ei, t;
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;

	/* set parameters to default values */
	m0 = FB_M0;
	m10 = FB_M10;
	m11 = FB_M11;
	r0 = FB_R0;
	r10 = FB_R10;
	r11 = FB_R11;
	a1 = FB_A1;
	e1 = FB_E1;
	vinf = FB_VINF;
	b = FB_B;
	input.ks = FB_KS;
	input.tstop = FB_TSTOP;
	input.Dflag = 0;
	input.dt = FB_DT;
	input.tcpustop = FB_TCPUSTOP;
	input.absacc = FB_ABSACC;
	input.relacc = FB_RELACC;
	input.ncount = FB_NCOUNT;
	input.tidaltol = FB_TIDALTOL;
	input.fexp = FB_FEXP;
	seed = FB_SEED;
	fb_debug = FB_DEBUG;

	/* initialize a few things for integrator */
	t = 0.0;
	fbui_new_hier(&hier, 3);

	/* initialize GSL rng */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, seed);

	/* create binary */
	fbui_make_pair(fbui_hierarchy_binary(&hier, 0), 
		       fbui_hierarchy_single(&hier, 1), fbui_hierarchy_single(&hier, 2));
	fbui_obj_t_set(fbui_hierarchy_binary(&hier, 0), t);

	/* give the objects some properties */
	fbui_initialize_single(fbui_hierarchy_single(&hier, 0), 0, "0");
	fbui_initialize_single(fbui_hierarchy_single(&hier, 1), 1, "1");
	fbui_initialize_single(fbui_hierarchy_single(&hier, 2), 2, "2");
	
	fbui_obj_radius_set(fbui_hierarchy_single(&hier, 0), r0);
	fbui_obj_radius_set(fbui_hierarchy_single(&hier, 1), r10);
	fbui_obj_radius_set(fbui_hierarchy_single(&hier, 2), r11);

	fbui_obj_mass_set(fbui_hierarchy_single(&hier, 0), m0);
	fbui_obj_mass_set(fbui_hierarchy_single(&hier, 1), m10);
	fbui_obj_mass_set(fbui_hierarchy_single(&hier, 2), m11);

	fbui_obj_mass_set(fbui_hierarchy_binary(&hier, 0), m10+m11);
	
	fbui_obj_a_set(fbui_hierarchy_binary(&hier, 0), a1);
	fbui_obj_e_set(fbui_hierarchy_binary(&hier, 0), e1);

	/* get the units and normalize */
	calc_units(fbui_hierarchy_single(&hier, 0), fbui_hierarchy_binary(&hier, 0), &units);
	fb_normalize(&hier, units);
	
	/* move hierarchies analytically in from infinity along hyperbolic orbit */
	m0 = fbui_obj_mass_get(fbui_hierarchy_single(&hier, 0));
	m1 = fbui_obj_mass_get(fbui_hierarchy_binary(&hier, 0));
	M = m0 + m1;
	mu = m0 * m1 / M;

	Ei = 0.5 * mu * fb_sqr(vinf);

	a1 = fbui_obj_a_get(fbui_hierarchy_binary(&hier, 0));
	e1 = fbui_obj_e_get(fbui_hierarchy_binary(&hier, 0));
	m10 = fbui_obj_mass_get(fbui_obj_left_child_get(fbui_hierarchy_binary(&hier, 0)));
	m11 = fbui_obj_mass_get(fbui_obj_right_child_get(fbui_hierarchy_binary(&hier, 0)));

	rtid = pow(2.0*(m0+m1)/(m1*input.tidaltol), 1.0/3.0) * a1 * (1.0+e1);

	fb_init_scattering(fbui_hierarchy_single(&hier, 0), fbui_hierarchy_binary(&hier, 0), vinf, b, rtid);

	/* trickle down the binary properties, then back up */
	fb_randorient(fbui_hierarchy_binary(&hier, 0), rng);
	fb_downsync(fbui_hierarchy_binary(&hier, 0), t);
	
	/* call fewbody! */
	retval = fewbody(input, &hier, &t);

	/* print information to screen */
	fprintf(stderr, "OUTCOME:\n");
	if (retval.retval == 1) {
		fprintf(stderr, "  encounter complete:  t=%.6g (%.6g yr)  %s  (%s)\n\n",
			t, t * units.t/FB_CONST_YR,
			fb_sprint_hier(hier, string1),
			fb_sprint_hier_hr(hier, string2));
	} else {
		fprintf(stderr, "  encounter NOT complete:  t=%.6g (%.6g yr)  %s  (%s)\n\n",
			t, t * units.t/FB_CONST_YR,
			fb_sprint_hier(hier, string1),
			fb_sprint_hier_hr(hier, string2));
	}

	fb_dprintf("there were %ld integration steps\n", retval.count);
	fb_dprintf("fb_classify() was called %ld times\n", retval.iclassify);
	
	fprintf(stderr, "FINAL:\n");
	fprintf(stderr, "  t_final=%.6g (%.6g yr)  t_cpu=%.6g s\n", \
		t, t*units.t/FB_CONST_YR, retval.tcpu);

	fprintf(stderr, "  DeltaL/L0=%.6g  DeltaL=%.6g\n", retval.DeltaLfrac, retval.DeltaL);
	fprintf(stderr, "  DeltaE/E0=%.6g  DeltaE=%.6g\n", retval.DeltaEfrac, retval.DeltaE);
	fprintf(stderr, "  Rmin=%.6g (%.6g RSUN)  Rmin_i=%d  Rmin_j=%d\n", \
		retval.Rmin, retval.Rmin*units.l/FB_CONST_RSUN, retval.Rmin_i, retval.Rmin_j);
	fprintf(stderr, "  Nosc=%d (%s)\n", retval.Nosc, (retval.Nosc>=1?"resonance":"non-resonance"));
	
	fprintf(stderr, "Resulting objects:\n");
	for (i=0; i<hier.nobj; i++) {
		n = fbui_obj_n_get(fbui_tree(&hier, i));
		if (n == 1) {
			fprintf(stderr, "  single: m=%g idstring=%s\n", 
				fbui_obj_mass_get(fbui_tree(&hier, i)) * units.m / FB_CONST_MSUN,
				fbui_obj_idstring_get(fbui_tree(&hier, i)));
		} else if (n == 2) {
			fprintf(stderr, "  binary: m_i=%g,%g a=%g e=%g idstring=%s\n",
				fbui_obj_mass_get(fbui_obj_left_child_get(fbui_tree(&hier, i))) * units.m / FB_CONST_MSUN,
				fbui_obj_mass_get(fbui_obj_right_child_get(fbui_tree(&hier, i))) * units.m / FB_CONST_MSUN,
				fbui_obj_a_get(fbui_tree(&hier, i)) * units.l / FB_CONST_AU, 
				fbui_obj_e_get(fbui_tree(&hier, i)), 
				fbui_obj_idstring_get(fbui_tree(&hier, i)));
		}
	}

	/* free GSL stuff */
	gsl_rng_free(rng);

	/* free our own stuff */
	fbui_delete_hier(&hier);

	/* done! */
	return(0);
}
