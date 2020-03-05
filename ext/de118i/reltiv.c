#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#include "const.h"
#include "ssystem.h"
#endif

/* Relativistic corrections.  Reference:
 *
 * Newhall, XX, E. M. Standish, and J. G. Williams, "DE102: a
 * numerically integrated ephemeris of the Moon and planets
 * spanning forty-four centuries," Astronomy and Astrophysics
 * 125, 150-167 (1983).
 */
	
extern DOUBLE C;
extern DOUBLE GMs[];
extern DOUBLE Rij[NTOTAL][NTOTAL];
extern DOUBLE Rij3[NTOTAL][NTOTAL];
extern DOUBLE vnewt[6*NTOTAL];
DOUBLE v2j[NTOTAL];

reltiv( y, v )
DOUBLE y[], v[];
{
DOUBLE xd, yd, zd, xs, ys, zs, xv, yv, zv, csqi;
DOUBLE rc, r, s, temp, rci;

int i, j, k, ii, jj;

csqi = One/(C*C);

/* velocity terms */
jj = 6*FMASS;
for(j=FMASS; j<NTOTAL; j++)
	{
	xv = y[jj];
	yv = y[jj+2];
	zv = y[jj+4];
	v2j[j] = (xv * xv + yv * yv + zv * zv) * csqi;
	jj += 6;
	}

ii = 6*FMASS;
for(i=FMASS; i<IAROIDS; i++)
	{
	v[ii] = Zero;
	v[ii+2] = Zero;
	v[ii+4] = Zero;
	rci = Zero;
	for(k=FMASS; k<NTOTAL; k++)
		{
		if( k == i )
			continue;
		rci += GMs[k]*Rij[i][k];
		}
	xs = Zero;
	ys = Zero;
	zs = Zero;
	jj = 6*FMASS;
	for(j=FMASS; j<NTOTAL; j++)
		{
		if( j == i )
			{
			jj += 6;
			continue; /* skip to next j if i = j */
			}

		xd = y[jj+1] - y[ii+1];
		yd = y[jj+3] - y[ii+3];
		zd = y[jj+5] - y[ii+5];

		rc = -Four * csqi * rci;

		s = Zero;
		for(k=FMASS; k<NTOTAL; k++)
			{
			if( k == j )
				continue;
			s += GMs[k]*Rij[j][k];
			}
		rc -= csqi * s;

		rc += v2j[i];
		rc += Two * v2j[j];

		r = y[ii] * y[jj] + y[ii+2] * y[jj+2] + y[ii+4] * y[jj+4];
		rc -= Four * csqi * r;

		s = -xd * y[jj] - yd * y[jj+2] - zd * y[jj+4];
		s *= Rij[i][j];
		rc -= OneandaHalf * csqi * s * s;

		s = xd * vnewt[jj] + yd * vnewt[jj+2] + zd * vnewt[jj+4];
		rc += Half * csqi * s;

		temp = GMs[j] * Rij3[i][j];

/*
 *		rc = temp * (One + rc );
 */
		rc = temp * rc;
		xs += xd * rc;
		ys += yd * rc;
		zs += zd * rc;

		s = -xd * (Four*y[ii] - Three*y[jj])
			- yd * (Four*y[ii+2] - Three*y[jj+2])
			- zd * (Four*y[ii+4] - Three*y[jj+4]);
		s = s * temp * csqi;
		xs += s * (y[ii] - y[jj]);
		ys += s * (y[ii+2] - y[jj+2]);
		zs += s * (y[ii+4] - y[jj+4]);

		temp = ThreeandaHalf * csqi * GMs[j] * Rij[i][j];
		xs += temp * vnewt[jj];
		ys += temp * vnewt[jj+2];
		zs += temp * vnewt[jj+4];
		jj += 6;
		}
	v[ii] = xs;
	v[ii+2] = ys;
	v[ii+4] = zs;
	v[ii+1] = y[ii];
	v[ii+3] = y[ii+2];
	v[ii+5] = y[ii+4];
	ii += 6;
	}
}

