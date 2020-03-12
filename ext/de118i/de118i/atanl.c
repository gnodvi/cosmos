/*							atanl.c
 *
 *	Inverse circular tangent, long double precision
 *      (arctangent)
 *
 *
 *
 * SYNOPSIS:
 *
 * long double x, y, atanl();
 *
 * y = atanl( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between -pi/2 and +pi/2 whose tangent
 * is x.
 *
 * Range reduction is from four intervals into the interval
 * from zero to  tan( pi/8 ).  The approximant uses a rational
 * function of degree 3/4 of the form x + x**3 P(x)/Q(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10, 10    150000       1.3e-19     3.0e-20
 *
 */
/*							atan2l()
 *
 *	Quadrant correct inverse circular tangent,
 *	long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * long double x, y, z, atan2l();
 *
 * z = atan2l( y, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle whose tangent is y/x.
 * Define compile time symbol ANSIC = 1 for ANSI standard,
 * range -PI < z <= +PI, args (y,x); else ANSIC = 0 for range
 * 0 to 2PI, args (x,y).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10, 10     60000       1.7e-19     3.2e-20
 * See atan.c.
 *
 */

/*							atan.c */


/*
Cephes Math Library Release 2.2:  December, 1990
Copyright 1984, 1990 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#ifndef NOINCS
#include "mconf.h"
#endif

/* long double polevll(), p1evll(); */
extern long double p1evll ( long double x, long double *P, int n );
extern long double polevll ( long double x, long double *P, int n );
long double atanl ( long double x );

#ifdef UNK
long double atanP[] = {
-8.6863818178092187535440E-1L,
-1.4683508633175792446076E1L,
-6.3976888655834347413154E1L,
-9.9988763777265819915721E1L,
-5.0894116899623603312185E1L,
};
long double atanQ[] = {
/* 1.00000000000000000000E0L,*/
 2.2981886733594175366172E1L,
 1.4399096122250781605352E2L,
 3.6144079386152023162701E2L,
 3.9157570175111990631099E2L,
 1.5268235069887081006606E2L,
};

/* tan( 3*pi/8 ) */
long double T3P8 =  2.41421356237309504880169L;

/* tan( pi/8 ) */
long double TP8 =  4.1421356237309504880169L;
#endif


#ifdef IBMPC
short atanP[] = {
0x8ece,0xce53,0x1266,0xde5f,0xbffe, XPD
0x07e6,0xa061,0xa6bf,0xeaef,0xc002, XPD
0x53ee,0xf291,0x557f,0xffe8,0xc004, XPD
0xf9d6,0xeda6,0x3f3e,0xc7fa,0xc005, XPD
0xb6c3,0x6abc,0x9361,0xcb93,0xc004, XPD
};
short atanQ[] = {
/*0x0000,0x0000,0x0000,0x8000,0x3fff,*/
0x54d4,0x894e,0xe76e,0xb7da,0x4003, XPD
0x76b9,0x7a46,0xafa2,0x8ffd,0x4006, XPD
0xe3a9,0xe9c0,0x6bee,0xb4b8,0x4007, XPD
0xabc1,0x50a7,0xb098,0xc3c9,0x4007, XPD
0x891c,0x100d,0xae89,0x98ae,0x4006, XPD
};

/* tan( 3*pi/8 ) = 2.41421356237309504880 */
short T3P8A[] = {0x3242,0xfcef,0x7999,0x9a82,0x4000, XPD};
#define T3P8 *(long double *)T3P8A

/* tan( pi/8 ) = 0.41421356237309504880 */
short TP8A[] = {0x9211,0xe779,0xcccf,0xd413,0x3ffd, XPD};
#define TP8 *(long double *)TP8A
#endif

#ifdef MIEEE
long atanP[] = {
0xbffe0000,0xde5f1266,0xce538ece,
0xc0020000,0xeaefa6bf,0xa06107e6,
0xc0040000,0xffe8557f,0xf29153ee,
0xc0050000,0xc7fa3f3e,0xeda6f9d6,
0xc0040000,0xcb939361,0x6abcb6c3,
};
long atanQ[] = {
/*0x3fff0000,0x80000000,0x00000000,*/
0x40030000,0xb7dae76e,0x894e54d4,
0x40060000,0x8ffdafa2,0x7a4676b9,
0x40070000,0xb4b86bee,0xe9c0e3a9,
0x40070000,0xc3c9b098,0x50a7abc1,
0x40060000,0x98aeae89,0x100d891c,
};

/* tan( 3*pi/8 ) = 2.41421356237309504880 */
long T3P8A[] = {0x40000000,0x9a827999,0xfcef3242};
#define T3P8 *(long double *)T3P8A

/* tan( pi/8 ) = 0.41421356237309504880 */
long TP8A[] = {0x3ffd0000,0xd413cccf,0xe7799211};
#define TP8 *(long double *)TP8A
#endif
extern long double Zero, One, Two, Three;

long double atanl(x)
long double x;
{
extern long double PIO2L, PIO4L;
long double y, z;
short sign;

/* make argument positive and save the sign */
sign = 1;
if( x < 0 )
	{
	sign = -1;
	x = -x;
	}

/* range reduction */
if( x > T3P8 )
	{
	y = PIO2L;
	x = -( One/x );
	}

else if( x > TP8 )
	{
	y = PIO4L;
	x = (x-One)/(x+One);
	}
else
	y = Zero;

/* rational form in x**2 */
z = x * x;
y = y + ( polevll( z, atanP, 4 ) / p1evll( z, atanQ, 5 ) ) * z * x + x;

if( sign < 0 )
	y = -y;

return(y);
}

/*							atan2	*/


extern long double PIL, PIO2L;

#if ANSIC
long double atan2l( y, x )
#else
long double atan2l( x, y )
#endif
long double x, y;
{
long double z, w;
short code;


code = 0;

if( x < 0 )
	code = 2;
if( y < 0 )
	code |= 1;

if( x == 0 )
	{
	if( code & 1 )
		{
#if ANSIC
		return( -PIO2L );
#else
		return( Three*PIO2L );
#endif
		}
	if( y == 0 )
		return( Zero );
	return( PIO2L );
	}

if( y == 0 )
	{
	if( code & 2 )
		return( PIL );
	return( Zero );
	}

#if ANSIC
if( code < 2 )
	w = Zero;
else if( code == 2 )
	w = PIL;
else
	w = -PIL;
#else
if( code == 0 )
	w = Zero;
else if( code == 1 )
	w = Two * PIL;
else
	w = PIL;
#endif

z = atanl( y/x );

return( w + z );
}
