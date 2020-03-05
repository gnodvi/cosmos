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



void init_sys_8( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp )
{
    int		i, k;

    sp->star_distrib  = 0;
    sp->size          = 3.0300000000000000e+02;
    sp->num_collapsar = 0;
    sp->mono_stars    = 8;

    sp->star_circle   = 0;
    sp->cir_dist      = 0;

    sp->star_line     = 0;

    sp->do_bounce     = 0;
    sp->no_speed      = 0;
    sp->few_stars     = 0;
    sp->drift         = 0;
    sp->min_angular_mom  = 0;

    sp->min_stars     = 1;
    sp->num_add       = 0;
    sp->live_stars    = 8;

/* star 0 */

    m[0]        = 1.1468208092485550e+00;
    cur[0].x    = -2.7038883602489072e+00;
    cur[0].y    = 6.4089177018335889e+01;
    prev[0].x   = -2.7038883602489072e+00;
    prev[0].y   = 6.4089177018335889e+01;
    vk[0][0].x  = 1.3226751051071546e-02;
    vk[0][0].y  = -2.2643948440514610e-03;
    ak[0][0].x  = -2.1461637019132762e-06;
    ak[0][0].y  = -3.8329229061867005e-06;

/* star 1 */

    m[1]        = 2.2566473988439308e+00;
    cur[1].x    = 6.0566200308691172e+01;
    cur[1].y    = 5.0696810406423374e+01;
    prev[1].x   = 6.0566200308691172e+01;
    prev[1].y   = 5.0696810406423374e+01;
    vk[0][1].x  = 1.1326416168617813e-02;
    vk[0][1].y  = -1.2983888267850539e-02;
    ak[0][1].x  = -3.7978097624060328e-06;
    ak[0][1].y  = -4.9601977809216730e-06;

/* star 2 */

    m[2]        = 1.0358381502890175e+00;
    cur[2].x    = -1.3325368021270076e+01;
    cur[2].y    = 4.3795680059746633e+00;
    prev[2].x   = -1.3325368021270076e+01;
    prev[2].y   = 4.3795680059746633e+00;
    vk[0][2].x  = -1.4657477777562395e-03;
    vk[0][2].y  = 8.3098025452551811e-03;
    ak[0][2].x  = 1.4837373960419059e-06;
    ak[0][2].y  = 3.2846997147130252e-06;

/* star 3 */

    m[3]        = 2.7560693641618497e+00;
    cur[3].x    = 5.0296086642545227e+01;
    cur[3].y    = 3.3756236820588086e+00;
    prev[3].x   = 5.0296086642545227e+01;
    prev[3].y   = 3.3756236820588086e+00;
    vk[0][3].x  = -3.7597216177549879e-03;
    vk[0][3].y  = -8.3720019727285024e-03;
    ak[0][3].x  = -1.9157531920650712e-06;
    ak[0][3].y  = 3.8658856402276140e-06;

/* star 4 */

    m[4]        = 2.1271676300578037e+00;
    cur[4].x    = 6.8270243736874789e+01;
    cur[4].y    = -1.1595121023145049e+02;
    prev[4].x   = 6.8270243736874789e+01;
    prev[4].y   = -1.1595121023145049e+02;
    vk[0][4].x  = -1.4273624489641538e-02;
    vk[0][4].y  = -8.7326800300004132e-03;
    ak[0][4].x  = -9.7523640637911302e-07;
    ak[0][4].y  = 1.7221688018067123e-06;

/* star 5 */

    m[5]        = 2.0531791907514454e+00;
    cur[5].x    = -6.1922907126650045e+01;
    cur[5].y    = 8.2905813493962910e+01;
    prev[5].x   = -6.1922907126650045e+01;
    prev[5].y   = 8.2905813493962910e+01;
    vk[0][5].x  = 1.4057051300353198e-02;
    vk[0][5].y  = 9.7559665347273592e-03;
    ak[0][5].x  = 4.3248908536138410e-06;
    ak[0][5].y  = -6.5951618951828798e-06;

/* star 6 */

    m[6]        = 1.9236994219653181e+00;
    cur[6].x    = -4.8033093344416187e+01;
    cur[6].y    = 4.7578502331255180e+01;
    prev[6].x   = -4.8033093344416187e+01;
    prev[6].y   = 4.7578502331255180e+01;
    vk[0][6].x  = 3.8568338976656485e-03;
    vk[0][6].y  = 9.4050873297924203e-03;
    ak[0][6].x  = 2.2103324433326556e-06;
    ak[0][6].y  = 3.6631219712243281e-06;

/* star 7 */

    m[7]        = 2.7005780346820809e+00;
    cur[7].x    = -6.8161132128329925e+01;
    cur[7].y    = -8.0295095516064634e+01;
    prev[7].x   = -6.8161132128329925e+01;
    prev[7].y   = -8.0295095516064634e+01;
    vk[0][7].x  = -1.2873818798120700e-02;
    vk[0][7].y  = 9.9296098370418636e-03;
    ak[0][7].x  = 1.3764904099135169e-06;
    ak[0][7].y  = 1.6155758780025582e-06;


    /* scale things to the current accuracy */

    for( i = 0; i < max_stars; i++ )
    {
	if( m[i] <= 0 )
	    continue;
	
	m[i] *= fm;

	for( k = 0; k < K; k++ )
	{
	    vk[k][i].x *= fv_inv;
	    vk[k][i].y *= fv_inv;

	    ak[k][i].x *= fv_inv*fv_inv;
	    ak[k][i].y *= fv_inv*fv_inv;
	}
    }
}
