/* -*- linux-c -*- */
/* fewbody_ui.c

   Copyright (C) 2010 John M. Fregeau
   
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

   A rudimentary user interface for Fewbody.
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "fewbody.h"

/* memory management, initialization */
int fbui_new_hier(fb_hier_t *hier, int n)
{
	if (n <= 0) {
		fprintf(stderr, "n must be > 0!\n");
		exit(1);
	}
	hier->nstarinit = n;
	hier->nstar = n;
	fb_malloc_hier(hier);
	fb_init_hier(hier);
	return(0);
}

int fbui_delete_hier(fb_hier_t *hier)
{
	fb_free_hier(*hier);
	return(0);
}

/* element access */
fb_obj_t *fbui_hierarchy_element(fb_hier_t *hier, int n, int m)
{
	if (n > hier->nstarinit) {
		fprintf(stderr, "n > initial number of stars: no hierarchies that large possible!\n");
		exit(1);
	}
	
	if (m < 0) {
		fprintf(stderr, "m must be >= 0!  Elements are indexed relative to 0.\n");
		exit(1);
	}

	if (hier->hi[n]+m >= hier->hi[n+1]) {
		fprintf(stderr, "m > number of possible hierarchies of size %d!\n", n);
		exit(1);
	}

	return(&(hier->hier[hier->hi[n]+m]));
}

/* shortcuts for accessing certain hierarchy sizes */
fb_obj_t *fbui_hierarchy_single(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 1, m));
}

fb_obj_t *fbui_hierarchy_binary(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 2, m));
}

fb_obj_t *fbui_hierarchy_triple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 3, m));
}

fb_obj_t *fbui_hierarchy_quadruple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 4, m));
}

fb_obj_t *fbui_hierarchy_quintuple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 5, m));
}

fb_obj_t *fbui_hierarchy_sextuple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 6, m));
}

fb_obj_t *fbui_hierarchy_septuple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 7, m));
}

fb_obj_t *fbui_hierarchy_octuple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 8, m));
}

fb_obj_t *fbui_hierarchy_nonuple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 9, m));
}

fb_obj_t *fbui_hierarchy_decuple(fb_hier_t *hier, int m)
{
	return(fbui_hierarchy_element(hier, 10, m));
}

int fbui_make_pair(fb_obj_t *parentobj, fb_obj_t *obj1, fb_obj_t *obj2)
{
	parentobj->obj[0] = obj1;
	parentobj->obj[1] = obj2;
	return(0);
}

int fbui_initialize_single(fb_obj_t *single, long id, char idstring[FB_MAX_STRING_LENGTH])
{
	single->ncoll = 1;
	single->id[0] = id;
	snprintf(single->idstring, FB_MAX_STRING_LENGTH, "%s", idstring);
	single->n = 1;
	single->obj[0] = NULL;
	single->obj[1] = NULL;
	single->Eint = 0.0;
	single->Lint[0] = 0.0;
	single->Lint[1] = 0.0;
	single->Lint[2] = 0.0;
	return(0);
}

/* let's call hier.obj the tree list to make it simpler */
fb_obj_t *fbui_tree(fb_hier_t *hier, int n)
{
	if (n >= hier->nstar) {
		fprintf(stderr, "n=%d > possible number of trees\n", n);
		exit(1);
	}

	return(hier->obj[n]);
}

/* ncoll */
int fbui_obj_ncoll_set(fb_obj_t *obj, int ncoll)
{
	obj->ncoll = ncoll;
	return(0);
}

int fbui_obj_ncoll_get(fb_obj_t *obj)
{
	return(obj->ncoll);
}

/* id array */
int fbui_obj_id_set(fb_obj_t *obj, int index, long id)
{
	obj->id[index] = id;
	return(0);
}

long fbui_obj_id_get(fb_obj_t *obj, int index)
{
	return(obj->id[index]);
}

/* idstring */
int fbui_obj_idstring_set(fb_obj_t *obj, char idstring[FB_MAX_STRING_LENGTH])
{
	snprintf(obj->idstring, FB_MAX_STRING_LENGTH, "%s", idstring);
	return(0);
}

char *fbui_obj_idstring_get(fb_obj_t *obj)
{
	return(obj->idstring);
}

/* mass */
int fbui_obj_mass_set(fb_obj_t *obj, double mass)
{
	obj->m = mass;
	return(0);
}

double fbui_obj_mass_get(fb_obj_t *obj)
{
	return(obj->m);
}

/* radius */
int fbui_obj_radius_set(fb_obj_t *obj, double radius)
{
	obj->R = radius;
	return(0);
}

double fbui_obj_radius_get(fb_obj_t *obj)
{
	return(obj->R);
}

/* Eint */
int fbui_obj_Eint_set(fb_obj_t *obj, double Eint)
{
	obj->Eint = Eint;
	return(0);
}

double fbui_obj_Eint_get(fb_obj_t *obj)
{
	return(obj->Eint);
}

/* Lint */
int fbui_obj_Lint_set(fb_obj_t *obj, double Lint[3])
{
	obj->Lint[0] = Lint[0];
	obj->Lint[1] = Lint[1];
	obj->Lint[2] = Lint[2];
	return(0);
}

int fbui_obj_Linti_set(fb_obj_t *obj, double Lint0, double Lint1, double Lint2)
{
	obj->Lint[0] = Lint0;
	obj->Lint[1] = Lint1;
	obj->Lint[2] = Lint2;
	return(0);
}

double *fbui_obj_Lint_get(fb_obj_t *obj)
{
	return(obj->Lint);
}

/* position */
int fbui_obj_x_set(fb_obj_t *obj, double x[3])
{
	obj->x[0] = x[0];
	obj->x[1] = x[1];
	obj->x[2] = x[2];
	return(0);
}

int fbui_obj_xi_set(fb_obj_t *obj, double x0, double x1, double x2)
{
	obj->x[0] = x0;
	obj->x[1] = x1;
	obj->x[2] = x2;
	return(0);
}

double *fbui_obj_x_get(fb_obj_t *obj)
{
	return(obj->x);
}

/* velocity */
int fbui_obj_v_set(fb_obj_t *obj, double v[3])
{
	obj->v[0] = v[0];
	obj->v[1] = v[1];
	obj->v[2] = v[2];
	return(0);
}

int fbui_obj_vi_set(fb_obj_t *obj, double v0, double v1, double v2)
{
	obj->v[0] = v0;
	obj->v[1] = v1;
	obj->v[2] = v2;
	return(0);
}

double *fbui_obj_v_get(fb_obj_t *obj)
{
	return(obj->v);
}

/* n */
int fbui_obj_n_set(fb_obj_t *obj, int n)
{
	obj->n = n;
	return(0);
}

int fbui_obj_n_get(fb_obj_t *obj)
{
	return(obj->n);
}

/* obj array */
int fbui_obj_obj_set(fb_obj_t *obj, int i, fb_obj_t *objtopointto)
{
	obj->obj[i] = objtopointto;
	return(0);
}

int fbui_obj_left_child_set(fb_obj_t *obj, fb_obj_t *objtopointto)
{
	obj->obj[0] = objtopointto;
	return(0);
}

int fbui_obj_right_child_set(fb_obj_t *obj, fb_obj_t *objtopointto)
{
	obj->obj[1] = objtopointto;
	return(0);
}

fb_obj_t *fbui_obj_obj_get(fb_obj_t *obj, int i)
{
	return(obj->obj[i]);
}

fb_obj_t *fbui_obj_left_child_get(fb_obj_t *obj)
{
	return(obj->obj[0]);
}

fb_obj_t *fbui_obj_right_child_get(fb_obj_t *obj)
{
	return(obj->obj[1]);
}

/* a */
int fbui_obj_a_set(fb_obj_t *obj, double a)
{
	obj->a = a;
	return(0);
}

double fbui_obj_a_get(fb_obj_t *obj)
{
	return(obj->a);
}

/* e */
int fbui_obj_e_set(fb_obj_t *obj, double e)
{
	obj->e = e;
	return(0);
}

double fbui_obj_e_get(fb_obj_t *obj)
{
	return(obj->e);
}

/* Lhat */
int fbui_obj_Lhat_set(fb_obj_t *obj, double Lhat[3])
{
	obj->Lhat[0] = Lhat[0];
	obj->Lhat[1] = Lhat[1];
	obj->Lhat[2] = Lhat[2];
	return(0);
}

int fbui_obj_Lhati_set(fb_obj_t *obj, double Lhat0, double Lhat1, double Lhat2)
{
	obj->Lhat[0] = Lhat0;
	obj->Lhat[1] = Lhat1;
	obj->Lhat[2] = Lhat2;
	return(0);
}

double *fbui_obj_Lhat_get(fb_obj_t *obj)
{
	return(obj->Lhat);
}

/* Ahat */
int fbui_obj_Ahat_set(fb_obj_t *obj, double Ahat[3])
{
	obj->Ahat[0] = Ahat[0];
	obj->Ahat[1] = Ahat[1];
	obj->Ahat[2] = Ahat[2];
	return(0);
}

int fbui_obj_Ahati_set(fb_obj_t *obj, double Ahat0, double Ahat1, double Ahat2)
{
	obj->Ahat[0] = Ahat0;
	obj->Ahat[1] = Ahat1;
	obj->Ahat[2] = Ahat2;
	return(0);
}

double *fbui_obj_Ahat_get(fb_obj_t *obj)
{
	return(obj->Ahat);
}

/* t */
int fbui_obj_t_set(fb_obj_t *obj, double t)
{
	obj->t = t;
	return(0);
}

double fbui_obj_t_get(fb_obj_t *obj)
{
	return(obj->t);
}

/* mean_anom */
int fbui_obj_mean_anom_set(fb_obj_t *obj, double mean_anom)
{
	obj->mean_anom = mean_anom;
	return(0);
}

double fbui_obj_mean_anom_get(fb_obj_t *obj)
{
	return(obj->mean_anom);
}
