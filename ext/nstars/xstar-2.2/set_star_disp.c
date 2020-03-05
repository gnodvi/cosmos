#include "xstar.h"
#include "xstar_ext.h"


/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com)
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/



/*
 * Set up the star display variables
 */

void set_star_disp( star_disp_type *sd, sys_param_type *sp )
{
    /* set up default erase_hstep adjust points */
    sd->num_disp_pt_32 = num_disp_pt / 32;
    sd->num_disp_pt_16 = num_disp_pt / 16;
    sd->num_disp_pt_8 = num_disp_pt / 8;
    sd->num_disp_pt_4 = num_disp_pt / 4;

    /* handle special cases */
    if( sp->do_bounce && sp->star_line )
	sd->hsteps = sd->num_steps = MAX_HSTEP/2;
	   
    else if( sp->do_bounce )
    {
	if( sp->num_collapsar == 0 )
	    if( sp->no_speed )
		sd->hsteps = sd->num_steps = 400 / sp->live_stars + 20;
	    else
		sd->hsteps = sd->num_steps = 80 / sp->live_stars + 10;
	else
	    sd->hsteps = sd->num_steps = 20;

	sd->num_disp_pt_32 = num_disp_pt * .875 + num_disp_pt / (32*4);
	sd->num_disp_pt_16 = num_disp_pt * .875 + num_disp_pt / (16*4);
	sd->num_disp_pt_8 = num_disp_pt * .875 + num_disp_pt / (8*4);
	sd->num_disp_pt_4 = num_disp_pt * .875 + num_disp_pt / (4*4);
    }
    else if( sp->few_stars )
    {
	sd->hsteps = sd->num_steps = 80 / sp->live_stars + 10;

	sd->num_disp_pt_32 = num_disp_pt * .75 + num_disp_pt / (32*4);
	sd->num_disp_pt_16 = num_disp_pt * .75 + num_disp_pt / (16*4);
	sd->num_disp_pt_8 = num_disp_pt * .75 + num_disp_pt / (8*4);
	sd->num_disp_pt_4 = num_disp_pt * .75 + num_disp_pt / (4*4);
    }
    else if( sp->star_circle )
	sd->hsteps = sd->num_steps = 400;

    else
	sd->hsteps = sd->num_steps = 200;


    /* if we don't have many display points, then let the tails grow slowly */
    if( sd->num_disp_pt_4 - sd->num_disp_pt_8 < 1024 )
	sd->num_disp_pt_4 += 1024 - (sd->num_disp_pt_4 - sd->num_disp_pt_8);
    
    
    /* set up other misc star_disp variables */
    sd->erase_disps = 0;
    sd->live_erase = sp->live_stars;

    sd->next_free = 1;
    sd->next_disp = 1;
    sd->next_erase = 0;
    sd->erase_hstep = 0;
    sd->updt_hstep = sd->hsteps;

    sd->hstep_scale = fv_inv * (1./4);

    sd->raw_num_steps = sd->num_steps / sd->hstep_scale;
    sd->raw_hsteps = sd->hsteps / sd->hstep_scale;
}



void set_buffer_factor( star_disp_type *sd, sys_param_type *sp )
{
    int		total_stars = sp->live_stars + sp->num_collapsar;
    
    sd->buffer_factor = max_disp_skip - total_stars*total_stars;
}

