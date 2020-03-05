#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#include "ssystem.h"
#include "const.h"
#endif

#define FIDEBUG 0

extern DOUBLE C;
extern DOUBLE Rij[NTOTAL][NTOTAL];
extern DOUBLE GMs[];
extern DOUBLE EMRAT;
extern DOUBLE yw[];
extern DOUBLE ymoon[];



/* This routine adds a constant vector to move the
 * barycenter to the origin.
 */
extern DOUBLE mustar[NTOTAL];

findcenter( t, yin )
DOUBLE t;
DOUBLE yin[];
{
int i, j, k, ii, jj;
DOUBLE xx, yy, zz, vx, vy, vz, csqi, s, mu;

#if MOON
fromemb( yin, yw );
#endif
#if AROIDS
aroids( t, yw );
#endif
csqi = Half / (C*C);
for( k=0; k<3; k++ )
	{ /* Iterate to find solution of implicit expressions. */

	distances(yw);
#if DOREL
/* Relativistic GM of each body */
	ii = 6*FMASS;
	for( i=FMASS; i<NTOTAL; i++ )
		{
		vx = yw[ii]; /* velocity */
		vy = yw[ii+2];
		vz = yw[ii+4];
		s = vx * vx + vy * vy + vz * vz; /* square of velocity */
		for( j=FMASS; j<NTOTAL; j++ )
			{
			if( j == i )
				continue;
			s -= GMs[j]*Rij[i][j];
			}
		mustar[i] = GMs[i] * (One + csqi * s);
		ii += 6;
		}
#endif

/* Compute relativistic (or Newtonian) center of mass. */
	mu = Zero; /* total mass of the system */
	xx = Zero;
	yy = Zero;
	zz = Zero;
	vx = Zero;
	vy = Zero;
	vz = Zero;
	ii = 6*FMASS;
	for( i=FMASS; i<NTOTAL; i++ )
		{
#if DOREL
		s = mustar[i];
#else
		s = GMs[i];
#endif
		mu += s;
		xx += yw[ii+1] * s; /* position */
		yy += yw[ii+3] * s;
		zz += yw[ii+5] * s;
		vx += yw[ii] * s; /* velocity */
		vy += yw[ii+2] * s;
		vz += yw[ii+4] * s;
		ii += 6;
		}
/* Offset the coordinates of all the objects equally. */
	xx /= mu;
	yy /= mu;
	zz /= mu;
	vx /= mu;
	vy /= mu;
	vz /= mu;
	jj = 6*FMASS;
	for( i=FMASS; i<NTOTAL; i++ )
		{
		yw[jj+1] -= xx;
		yw[jj+3] -= yy;
		yw[jj+5] -= zz;
		yw[jj] -= vx;
		yw[jj+2] -= vy;
		yw[jj+4] -= vz;
		jj += 6;
		}
#if FIDEBUG
	if( k == 0 )
		{
		xx = xx*xx + yy*yy + zz*zz;
		vx = vx*vx + vy*vy + vz*vz;
		printf( "Barycenter offset ^2 = %.5Le au, %.5Le au/d\n",
			xx, vx );
		}
#endif
	}

#if MOON
jj = 6*FMASS;
for( i=FMASS; i<NTOTAL; i++ )
	{
	if( i == IMOON )
		jj += 6;
	else
		{
		for( j=0; j<6; j++ )
			{
			yin[jj] = yw[jj];
			jj += 1;
			}
		}
	}
/* Convert barycentric Earth and Moon to output EMB variables. */
ii = 6*IEARTH;
jj = 6*IMOON;
for( i=0; i<6; i++ )
	{
	xx = yw[ii+i]; /* SS barycentric Earth */
	yy = yw[jj+i]; /* SS barycentric Moon */
	yin[ii+i] = (EMRAT * xx + yy)/(EMRAT+One); /* EMB */
/*	yin[jj+i] = yy - xx; */ /* M = Moon - Earth */
	yin[jj+i] = ymoon[i];
	}
#endif

#if FIDEBUG
printf( "Sun's state:\n" );
for( i=0; i<6; i++ )
	{
	xx = yw[6*ISUN+i];
	printf( "%27.18Le,\n", xx );
	}
#endif
}


#if FIDEBUG
/* Initial heliocentric accelerations, from DE118
 * JED 2440400.5
 */
DOUBLE iniaccel[] = {
/* Libs */
-2.33840106202781591e-6L,
 3.76106870289525089e-6L,
-8.36366071885551190e-7L,
/* Mercury */
-1.93862993974248726e-3L,
 5.20550192322063647e-4L,
 4.77930057349665243e-4L,

-4.62592548867567062e-4L,
 2.72876752608682010e-4L,
 1.52182594409640823e-4L,

-2.92050891691249361e-5L,
 2.61318585392418325e-4L,
 1.13316577018344841e-4L,

 1.24770701161065107e-5L,
 1.24988885615999961e-4L,
 5.70375473684709056e-5L,

 9.85358573091297212e-6L,
 1.40990437369746227e-6L,
 3.64122840684491628e-7L,

-2.95364938480744481e-6L,
-1.67622531774450414e-6L,
-5.65095232380545449e-7L,

 8.90424390777303382e-7L,
 4.87032002931296781e-8L,
 8.78848553654263136e-9L,

 1.80994923472092847e-7L,
 2.55021868335872720e-7L,
 1.00082692642110886e-7L,

 2.88404461789641877e-7L,
 7.37133065237483562e-9L,
-8.27405630423931700e-8L,
/* Moon */
 5.40548681996079620e-5L,
 1.26326595513576372e-4L,
 6.91101026282511274e-5L,
/* Sun */
0.0e0L,
0.0e0L,
0.0e0L,
};

chkacc(v)
DOUBLE v[];
{
DOUBLE asunx, asuny, asunz, dx, dy, dz;
double px, py, pz;
int i, ii, jj;

/* estimate the initial acceleration of the Sun
 * from the given heliocentric acceleration of Pluto
 */
ii = 6*(IMOON-1);
jj = 3*(IMOON-1);
asunx = v[ii] - iniaccel[jj];
asuny = v[ii+2] - iniaccel[jj+1];
asunz = v[ii+4] - iniaccel[jj+2];

printf( "\n" );
ii = 0;
jj = 0;
for( i=0; i<11; i++ )
	{
	dx = v[ii]   - iniaccel[jj];
	dy = v[ii+2] - iniaccel[jj+1];
	dz = v[ii+4] - iniaccel[jj+2];
	if( (i != 0) && (i != IMOON) )
		{
		dx -= asunx;
		dy -= asuny;
		dz -= asunz;
		}
	px = dx;
	py = dy;
	pz = dz;
	printf( "%11.3e %11.3e %11.3e\n", px, py, pz );
	ii += 6;
	jj += 3;
	}
}
#endif
