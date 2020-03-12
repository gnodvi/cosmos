/*							tanl.c
 *
 *	Circular tangent, long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * long double x, y, tanl();
 *
 * y = tanl( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular tangent of the radian argument x.
 *
 * Range reduction is modulo pi/4.  A rational function
 *       x + x**3 P(x**2)/Q(x**2)
 * is employed in the basic interval [0, pi/4].
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     +-1.07e9       30000     1.9e-19     4.8e-20
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * tan total loss   x > 2^39                0.0
 *
 */
/*							cotl.c
 *
 *	Circular cotangent, long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * long double x, y, cotl();
 *
 * y = cotl( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular cotangent of the radian argument x.
 *
 * Range reduction is modulo pi/4.  A rational function
 *       x + x**3 P(x**2)/Q(x**2)
 * is employed in the basic interval [0, pi/4].
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     +-1.07e9      30000      1.9e-19     5.1e-20
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * cot total loss   x > 2^39                0.0
 * cot singularity  x = 0                  MAXNUM
 *
 */

/*
Cephes Math Library Release 2.2:  December, 1990
Copyright 1984, 1990 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#ifndef NOINCS
#include "mconf.h"
#endif

/* long double polevll(), p1evll(), floorl(), ldexpl(); */
extern long double p1evll ( long double x, long double *P, int n );
extern long double polevll ( long double x, long double *P, int n );
extern long double floorl ( long double x );
extern long double ldexpl ( long double x, int pw2 );

#ifdef UNK
long double tanP[] = {
-1.3093693918138377764608E4L,
 1.1535166483858741613983E6L,
-1.7956525197648487798769E7L,
};
long double tanQ[] = {
/* 1.0000000000000000000000E0L,*/
 1.3681296347069295467845E4L,
-1.3208923444021096744731E6L,
 2.5008380182335791583922E7L,
-5.3869575592945462988123E7L,
};
long double tanDP1 = 7.853981554508209228515625E-1L;
long double tanDP2 = 7.946627356147928367136046290398E-9L;
long double tanDP3 = 3.061616997868382943065164830688E-17L;
long double tanlossth = 5.49755813888e11L; /* 2^39 */
#endif


#ifdef IBMPC
short tanP[] = {
0xbc1c,0x79f9,0xc692,0xcc96,0xc00c, XPD
0xe5b1,0xe4ee,0x652f,0x8ccf,0x4013, XPD
0xaf9a,0x4c8b,0x5699,0x88ff,0xc017, XPD
};
short tanQ[] = {
/*0x0000,0x0000,0x0000,0x8000,0x3fff,*/
0x8ed4,0x9b2b,0x2f75,0xd5c5,0x400c, XPD
0xadcd,0x55e4,0xe2c1,0xa13d,0xc013, XPD
0x7adf,0x56c7,0x7e17,0xbecc,0x4017, XPD
0x86f6,0xf2d1,0x01e5,0xcd7f,0xc018, XPD
};
short tanP1[] = {0x0000,0x0000,0xda80,0xc90f,0x3ffe, XPD};
short tanP2[] = {0x0000,0x0000,0xa300,0x8885,0x3fe4, XPD};
short tanP3[] = {0x3707,0xa2e0,0x3198,0x8d31,0x3fc8, XPD};
#define tanDP1 *(long double *)tanP1
#define tanDP2 *(long double *)tanP2
#define tanDP3 *(long double *)tanP3
short tanlosi[] = {0x0000,0x0000,0x0000,0x8000,0x4026, XPD};
#define tanlossth *(long double *)tanlosi
#endif

#ifdef MIEEE
long tanP[] = {
0xc00c0000,0xcc96c692,0x79f9bc1c,
0x40130000,0x8ccf652f,0xe4eee5b1,
0xc0170000,0x88ff5699,0x4c8baf9a,
};
long tanQ[] = {
/*0x3fff0000,0x80000000,0x00000000,*/
0x400c0000,0xd5c52f75,0x9b2b8ed4,
0xc0130000,0xa13de2c1,0x55e4adcd,
0x40170000,0xbecc7e17,0x56c77adf,
0xc0180000,0xcd7f01e5,0xf2d186f6,
};
long tanP1[] = {0x3ffe0000,0xc90fda80,0x00000000};
long tanP2[] = {0x3fe40000,0x8885a300,0x00000000};
long tanP3[] = {0x3fc80000,0x8d313198,0xa2e03707};
#define tanDP1 *(long double *)tanP1
#define tanDP2 *(long double *)tanP2
#define tanDP3 *(long double *)tanP3
long tanlosi[] = {0x40260000,0x80000000,0x00000000};
#define tanlossth *(long double *)tanlosi
#endif

extern long double PIO4L, MACHEPL, MAXNUML;
extern long double Zero, One;

long double tancotl();

long double tanl(x)
long double x;
{

return( tancotl(x,0) );
}


long double cotl(x)
long double x;
{

if( x == 0 )
	{
	mtherr( "cotl", SING );
	return( MAXNUML );
	}
return( tancotl(x,1) );
}


long double tancotl( xx, cotflg )
long double xx;
int cotflg;
{
long double x, y, z, zz;
int j, sign;

/* make argument positive but save the sign */
if( xx < 0 )
	{
	x = -xx;
	sign = -1;
	}
else
	{
	x = xx;
	sign = 1;
	}

if( x > tanlossth )
	{
	if( cotflg )
		mtherr( "cotl", TLOSS );
	else
		mtherr( "tanl", TLOSS );
	return(Zero);
	}

/* compute x mod PIO4 */
y = floorl( x/PIO4L );

/* strip high bits of integer part */
z = ldexpl( y, -4 );
z = floorl(z);		/* integer part of y/16 */
z = y - ldexpl( z, 4 );  /* y - 16 * (y/16) */

/* integer and fractional part modulo one octant */
j = (int )z;

/* map zeros and singularities to origin */
if( j & 1 )
	{
	j += 1;
	y += One;
	}

z = ((x - y * tanDP1) - y * tanDP2) - y * tanDP3;

zz = z * z;

if( zz > MACHEPL )
	y = z  +  z * (zz * polevll( zz, tanP, 2 )/p1evll(zz, tanQ, 4));
else
	y = z;
	
if( j & 2 )
	{
	if( cotflg )
		y = -y;
	else
		y = -One/y;
	}
else
	{
	if( cotflg )
		y = One/y;
	}

if( sign < 0 )
	y = -y;

return( y );
}
