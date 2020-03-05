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



/* update the screen with anything that hasn't been plotted */
void update_screen( sys_param_type *sp, star_disp_type *sd, disp *display )
{

    if( num_disp_pt )
    {
	int		i;

	int		num_pts, pp, next_color;
	int		pts_to_plot;
	int		free_pts;
	int		stop_erase, stop_err;
	int		new_steps;
	

	/*
	 * erase the old points
	 */

	/* need to have at least a minimal amount of slots free */
	stop_erase = sd->next_erase;
	
	free_pts = DIFF( sd->next_erase, sd->next_free );
	if( free_pts < sp->live_stars*2+2 )
	{
	    stop_erase = SUM( stop_erase, (sp->live_stars*2+2) - free_pts );

	    if( sd->erase_hstep < sd->disp_pts[ stop_erase ].hstep )
		sd->erase_hstep = sd->disp_pts[ stop_erase ].hstep;
	}

	new_steps = sd->num_steps - sd->updt_hstep;
	sd->erase_hstep += new_steps;

	/* modify erase_hstep based on how full/empty the buffer is */
	if( free_pts < sp->live_stars*8+2 )
	    sd->erase_hstep += new_steps << 2;

	else if( free_pts < sd->num_disp_pt_32 )
	    sd->erase_hstep += new_steps << 1;

	else if( free_pts < sd->num_disp_pt_16 )
	    sd->erase_hstep += new_steps;

	else if( free_pts < sd->num_disp_pt_8 )
	    sd->erase_hstep += new_steps >> 2;

	else if( free_pts < sd->num_disp_pt_4 )
	    sd->erase_hstep -= (new_steps >> 5) + 1;

	else
	    sd->erase_hstep -= new_steps >> 1;


	    
	/* find where to stop erasing stars at */
	stop_err = PREV( sd->next_erase );
	for( ;; stop_erase = NEXT( stop_erase ) )
	{
	    if( stop_erase == stop_err )
	    {
		if( verbose > 0 )
		    fprintf( stderr, "%s:%d step %7d  Error:  buffer overflow.\n", __FILE__, __LINE__, sd->raw_num_steps );
		stop_erase = PREV( stop_erase );
		sd->erase_hstep = sd->disp_pts[ stop_erase ].hstep;
		break;
	    }
	    
	    if( sd->disp_pts[ stop_erase ].star == DISP_PT_UNUSED )
		continue;

	    if( sd->disp_pts[ stop_erase ].hstep <= sd->erase_hstep )
		continue;

	    break;
	}
	

	/* erase the stars */
	if( stop_erase != sd->next_erase )
	{
	    num_pts = 0;
	    for (pp = sd->next_erase; pp != stop_erase; pp = NEXT(pp) )
	    {
		if( sd->disp_pts[pp].star == DISP_PT_UNUSED )
		    continue;
	    
		sd->points[num_pts++] = sd->disp_pts[pp].pt;

		/* the following is not needed and screws up the redraw code */
/*		sd->disp_pts[pp].star = DISP_PT_UNUSED; */
		
		/*
		 * delete this point out of the hash table
		 */
	    {
		int		hash_loc;
		int		hashed_idx;
		int		found = FALSE;
		int		xpt = sd->disp_pts[pp].pt.x;
		int		ypt = sd->disp_pts[pp].pt.y;
		
		
		for( hash_loc = PT_HASH( xpt, ypt ),
		    hashed_idx = sd->hash_index[ hash_loc ];
		    
		    hashed_idx != HASH_UNUSED;
		    
		    hash_loc = NEXT_H( hash_loc ),
		    hashed_idx = sd->hash_index[ hash_loc ] )
		{
		    if( hashed_idx == HASH_SEARCH )
			continue;
		    
		    if( sd->disp_pts[ hashed_idx ].pt.x != xpt
		       || sd->disp_pts[ hashed_idx ].pt.y != ypt
		       )
			continue;
		    
		    /* we have found the entry, now delete it */
		    found = TRUE;
		    if( sd->hash_index[ NEXT_H( hash_loc ) ] == HASH_UNUSED )
		    {
			do
			{
			    sd->hash_index[ hash_loc ] = HASH_UNUSED;
			    hash_loc = PREV_H( hash_loc );
			}
			while( sd->hash_index[ hash_loc ] == HASH_SEARCH );
		    }
		    else
			sd->hash_index[ hash_loc ] = HASH_SEARCH;
			
			
		    break;
		}

		if( !found && verbose > 0 )
		    fprintf( stderr, "%s:%d step %7d  Error: (%d,%d) was not found in the hash table!  hash_loc=%ld\n", __FILE__, __LINE__, sd->raw_num_steps, xpt, ypt, PT_HASH( xpt, ypt ) );
		
	    }
		
		/* see if the buffer is too full */
		if( num_pts >= sd->max_points )
		{
		    XDrawPoints( display->dpy, display->win, display->erase_gc,
				sd->points, num_pts, CoordModeOrigin );

		    num_pts = 0;
		}
	    }
	    
	    if( num_pts )
	    {
		XDrawPoints( display->dpy, display->win, display->erase_gc,
			    sd->points, num_pts, CoordModeOrigin );
		num_pts = 0;
	    }

	    if( sd->next_erase == sd->next_disp )
		sd->next_disp = pp;
	    sd->next_erase = pp;

	    free_pts = DIFF( sd->next_erase, sd->next_free );
	    if( verbose > 0 &&  free_pts < sp->live_stars*2+2 )
		fprintf( stderr, "%s:%d step%7d  free_pts=%d  min=%d  live_stars=%d\n", __FILE__, __LINE__, sd->num_steps, free_pts, sp->live_stars*2+2, sp->live_stars );
	}


	
	/*
	 * display the new points
	 */
	pts_to_plot = DIFF( sd->next_free, sd->next_disp );
	if( pts_to_plot > 0 )
	{
	    
	    for( i = 0; i < NUM_COLORS; i = next_color )
	    {
		next_color = NUM_COLORS;
		
		num_pts = 0;
		for (pp = sd->next_disp; pp != sd->next_free; pp = NEXT(pp) )
		{
		    if( sd->disp_pts[pp].star == DISP_PT_UNUSED )
			continue;
		    
		    if( rotate_colors )
		    {
			int	color = sd->disp_pts[pp].color;

			if( i == color )
			    sd->points[num_pts++] = sd->disp_pts[pp].pt;
			else if( i < color && color < next_color )
			    next_color = color;
		    }
		    else if( multi_colors )
		    {
			int	color = sd->star_color[ sd->disp_pts[pp].star ];
			if( i == color )
			    sd->points[num_pts++] = sd->disp_pts[pp].pt;
			else if( i < color && color < next_color )
			    next_color = color;
		    }
		    else
		    {
			sd->points[num_pts++] = sd->disp_pts[pp].pt;
		    }
		    
		    
		    if( num_pts >= sd->max_points )
		    {
			if( rotate_colors || multi_colors )
			    XDrawPoints( display->dpy, display->win,
					color_gcs[i], sd->points, num_pts,
					CoordModeOrigin );
			else 
			    XDrawPoints( display->dpy, display->win,
					display->star_gc, sd->points, num_pts,
					CoordModeOrigin );
			
			num_pts = 0;
		    }
		}
		
		if( num_pts )
		{
		    if( rotate_colors || multi_colors )
			XDrawPoints( display->dpy, display->win,
				    color_gcs[i], sd->points, num_pts,
				    CoordModeOrigin );
		    else 
			XDrawPoints( display->dpy, display->win,
				    display->star_gc, sd->points, num_pts,
				    CoordModeOrigin );
		    num_pts = 0;
		}
	    }
	    sd->next_disp = sd->next_free;
	    
	    
	    /* if we're not generating much data, then force it out quickly */
	    if( sd->num_disp_skipped > 5 + sd->buffer_factor
	       || (pts_to_plot >= 4*sd->num_visible
		   && sd->num_disp_skipped*4 > sd->buffer_factor
		   )
	       )
		XFlush( display->dpy );
	    
	    sd->points_plotted += pts_to_plot;
	    
	}

	sd->num_disp_skipped = 0;
	sd->updt_hstep = sd->num_steps;
    }	


    else
    {
	sd->erase_disps++;
	
	if( rotate_colors )
	    XDrawPoints(display->dpy, display->win,
			color_gcs[ sd->color_number ],
			sd->points, sd->points_used, CoordModeOrigin );
	else if( multi_colors )
	{
	    int		num_pts, pp, next_color;
	    int		i;
	    
	    
	    for( i = 0; i < NUM_COLORS; i = next_color )
	    {
		next_color = NUM_COLORS;
		
		num_pts = 0;
		for (pp = 0; pp < sd->points_used; pp++)
		{
		    if( i == sd->pixels[pp] )
			sd->tmp_pts[num_pts++] = sd->points[pp];
		    else if( i < sd->pixels[pp]
			    && sd->pixels[pp] < next_color )
			next_color = sd->pixels[pp];
		}
		
		if( num_pts )
		{
		    XDrawPoints(display->dpy, display->win,
				color_gcs[ i ],
				sd->tmp_pts, num_pts, CoordModeOrigin );
		    num_pts = 0;
		}
	    }
	}
	else
	    XDrawPoints(display->dpy, display->win,
			display->star_gc,
			sd->points, sd->points_used, CoordModeOrigin );
	    
	
	/* if we're not generating much data, then force it out quickly */
	if( sd->num_disp_skipped > 5 + sd->buffer_factor
	   || (sd->points_used >= 5*sd->num_visible
	       && sd->num_disp_skipped*4 > sd->buffer_factor
	       )
	   )
	    XFlush( display->dpy );
	
	sd->points_plotted += sd->points_used;
	sd->points_used = 0;
	sd->num_disp_skipped = 0;
    }
}



/* clear out the old data */
void clear_disp_pt( star_disp_type *sd )
{
    int		i;

    if( num_disp_pt )
    {
	/* mark the entries as unused */
	for( i = 0; i < num_disp_pt; i++ )
	    sd->disp_pts[i].star = DISP_PT_UNUSED;

	for( i = 0; i < HASH_TABLE_SIZE; i++ )
	    sd->hash_index[i] = HASH_UNUSED;
    }
    else
	sd->points_used = 0;

    sd->num_disp_skipped = 0;
    sd->num_poll_skipped = 0;
    sd->num_seen = 1;

    sd->points_plotted = 0;
    sd->total_points = 0;
}




/* redraw the screen from scratch */
void redraw_screen( point_2d *cur, double *m, sys_param_type *sp, star_disp_type *sd, disp *display )
{
    int		i;
    
    int		num_pts, pp, next_color;
    

    /* redraw the collapsars */
    plot_collapsars( cur, m, sp, sd, display );
    
	    
    /* redraw the stars */
    if( !num_disp_pt )
	return;
    
    for( i = 0; i < NUM_COLORS; i = next_color )
    {
	next_color = NUM_COLORS;
	
	num_pts = 0;
	for (pp = sd->next_erase; pp != sd->next_free; pp = NEXT(pp) )
	{
	    XPoint	*pt = &sd->disp_pts[pp].pt;

	    if( sd->disp_pts[pp].star == DISP_PT_UNUSED )
		continue;
		    
	    if( pt->x < sd->redraw_min_x || pt->x > sd->redraw_max_x
	       || pt->y < sd->redraw_min_y || pt->y > sd->redraw_max_y
	       )
		continue;

	    
	    if( rotate_colors )
	    {
		int	color = sd->disp_pts[pp].color;
		
		if( i == color )
		    sd->points[num_pts++] = sd->disp_pts[pp].pt;
		else if( i < color && color < next_color )
		    next_color = color;
	    }
	    else if( multi_colors )
	    {
		int	color = sd->star_color[ sd->disp_pts[pp].star ];

		if( i == color )
		    sd->points[num_pts++] = sd->disp_pts[pp].pt;
		else if( i < color && color < next_color )
		    next_color = color;
	    }
	    else
		sd->points[num_pts++] = sd->disp_pts[pp].pt;
	    
	    
	    if( num_pts >= sd->max_points )
	    {
		if( rotate_colors || multi_colors )
		    XDrawPoints( display->dpy, display->win,
				color_gcs[i], sd->points, num_pts,
				CoordModeOrigin );
		else 
		    XDrawPoints( display->dpy, display->win,
				display->star_gc, sd->points, num_pts,
				CoordModeOrigin );
		
		num_pts = 0;
	    }
	}
	
	if( num_pts )
	{
	    if( rotate_colors || multi_colors )
		XDrawPoints( display->dpy, display->win,
			    color_gcs[i], sd->points, num_pts,
			    CoordModeOrigin );
	    else 
		XDrawPoints( display->dpy, display->win,
			    display->star_gc, sd->points, num_pts,
			    CoordModeOrigin );
	    num_pts = 0;
	}
    }
    
    XFlush( display->dpy );
}


void set_redraw_full( star_disp_type *sd )
{
    sd->redraw_min_x = min_x + center_x;
    sd->redraw_max_x = max_x + center_x;
    sd->redraw_min_y = min_y + center_y;
    sd->redraw_max_y = max_y + center_y;
}

