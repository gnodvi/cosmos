/*							sinl.c
 *
 *	Circular sine, long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * long double x, y, sinl();
 *
 * y = sinl( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of pi/4.  The reduction
 * error is nearly eliminated by contriving an extended precision
 * modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by the Cody
 * and Waite polynomial form
 *      x + x**3 P(x**2) .
 * Between pi/4 and pi/2 the cosine is represented as
 *      1 - .5 x**2 + x**4 Q(x**2) .
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE     +-5.5e11      200,000    1.2e-19     2.9e-20
 * 
 * ERROR MESSAGES:
 *
 *   message           condition        value returned
 * sin total loss   x > 2**39               0.0
 *
 * Loss of precision occurs for x > 2**39 = 5.49755813888e11.
 * The routine as implemented flags a TLOSS error for
 * x > 2**39 and returns 0.0.
 */
/*							cosl.c
 *
 *	Circular cosine, long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * long double x, y, cosl();
 *
 * y = cosl( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of pi/4.  The reduction
 * error is nearly eliminated by contriving an extended precision
 * modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the cosine is approximated by
 *      1 - .5 x**2 + x**4 Q(x**2) .
 * Between pi/4 and pi/2 the sine is represented by the Cody
 * and Waite polynomial form
 *      x  +  x**3 P(x**2) .
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE     +-5.5e11       50000      1.2e-19     2.9e-20
 */

/*							sin.c	*/

/*
Cephes Math Library Release 2.2:  December, 1990
Copyright 1985, 1990 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#ifndef NOINCS
#include "mconf.h"
#endif

/* long double polevll(), floorl(), ldexpl(); */

extern long double polevll ( long double x, long double *P, int n );
extern long double floorl ( long double x );
extern long double ldexpl ( long double x, int pw2 );

#ifdef UNK
long double sincof[7] = {
-7.5785404094842805756289E-13L,
 1.6058363167320443249231E-10L,
-2.5052104881870868784055E-8L,
 2.7557319214064922217861E-6L,
-1.9841269841254799668344E-4L,
 8.3333333333333225058715E-3L,
-1.6666666666666666640255E-1L,
};
long double coscof[7] = {
 4.7377507964246204691685E-14L,
-1.1470284843425359765671E-11L,
 2.0876754287081521758361E-9L,
-2.7557319214999787979814E-7L,
 2.4801587301570552304991E-5L,
-1.3888888888888872993737E-3L,
 4.1666666666666666609054E-2L,
};
long double sinDP1 = 7.853981554508209228515625E-1L;
long double sinDP2 = 7.946627356147928367136046290398E-9L;
long double sinDP3 = 3.061616997868382943065164830688E-17L;
long double sinlossth = 5.49755813888e11L; /* 2^39 */
#endif

#ifdef IBMPC
short sincof[] = {
0x4e27,0xe1d6,0x2389,0xd551,0xbfd6, XPD
0x64d7,0xe706,0x4623,0xb090,0x3fde, XPD
0x01b1,0xbf34,0x2946,0xd732,0xbfe5, XPD
0xc8f7,0x9845,0x1d29,0xb8ef,0x3fec, XPD
0x6514,0x0c53,0x00d0,0xd00d,0xbff2, XPD
0x569a,0x8888,0x8888,0x8888,0x3ff8, XPD
0xaa97,0xaaaa,0xaaaa,0xaaaa,0xbffc, XPD
};
short coscof[] = {
0x7436,0x6f99,0x8c3a,0xd55e,0x3fd2, XPD
0x2f37,0x58f4,0x920f,0xc9c9,0xbfda, XPD
0x5350,0x659e,0xc648,0x8f76,0x3fe2, XPD
0x4d2b,0xf5c6,0x7dba,0x93f2,0xbfe9, XPD
0x53ed,0x0c66,0x00d0,0xd00d,0x3fef, XPD
0x7b67,0x0b60,0x60b6,0xb60b,0xbff5, XPD
0xaa9a,0xaaaa,0xaaaa,0xaaaa,0x3ffa, XPD
};
short sinP1[] = {0x0000,0x0000,0xda80,0xc90f,0x3ffe, XPD};
short sinP2[] = {0x0000,0x0000,0xa300,0x8885,0x3fe4, XPD};
short sinP3[] = {0x3707,0xa2e0,0x3198,0x8d31,0x3fc8, XPD};
#define sinDP1 *(long double *)sinP1
#define sinDP2 *(long double *)sinP2
#define sinDP3 *(long double *)sinP3
short sinlosi[] = {0x0000,0x0000,0x0000,0x8000,0x4026, XPD};
#define sinlossth *(long double *)sinlosi
#endif

#ifdef MIEEE
long sincof[] = {
0xbfd60000,0xd5512389,0xe1d64e27,
0x3fde0000,0xb0904623,0xe70664d7,
0xbfe50000,0xd7322946,0xbf3401b1,
0x3fec0000,0xb8ef1d29,0x9845c8f7,
0xbff20000,0xd00d00d0,0x0c536514,
0x3ff80000,0x88888888,0x8888569a,
0xbffc0000,0xaaaaaaaa,0xaaaaaa97,
};
long coscof[] = {
0x3fd20000,0xd55e8c3a,0x6f997436,
0xbfda0000,0xc9c9920f,0x58f42f37,
0x3fe20000,0x8f76c648,0x659e5350,
0xbfe90000,0x93f27dba,0xf5c64d2b,
0x3fef0000,0xd00d00d0,0x0c6653ed,
0xbff50000,0xb60b60b6,0x0b607b67,
0x3ffa0000,0xaaaaaaaa,0xaaaaaa9a,
};
long sinP1[] = {0x3ffe0000,0xc90fda80,0x00000000};
long sinP2[] = {0x3fe40000,0x8885a300,0x00000000};
long sinP3[] = {0x3fc80000,0x8d313198,0xa2e03707};
#define sinDP1 *(long double *)sinP1
#define sinDP2 *(long double *)sinP2
#define sinDP3 *(long double *)sinP3
long sinlosi[] = {0x40260000,0x80000000,0x00000000};
#define sinlossth *(long double *)sinlosi
#endif

extern long double PIO4L, Zero, One;


long double sinl(x)
long double x;
{
long double y, z, zz;
int j, sign;

/* make argument positive but save the sign */
sign = 1;
if( x < 0 )
	{
	x = -x;
	sign = -1;
	}

if( x > sinlossth )
	{
	mtherr( "sinl", TLOSS );
	return(Zero);
	}

y = floorl( x/PIO4L ); /* integer part of x/PIO4 */

/* strip high bits of integer part to prevent integer overflow */
z = ldexpl( y, -4 );
z = floorl(z);           /* integer part of y/8 */
z = y - ldexpl( z, 4 );  /* y - 16 * (y/16) */

j = (int )z; /* convert to integer for tests on the phase angle */
/* map zeros to origin */
if( j & 1 )
	{
	j += 1;
	y += One;
	}
j = j & 07; /* octant modulo 360 degrees */
/* reflect in x axis */
if( j > 3)
	{
	sign = -sign;
	j -= 4;
	}

/* Extended precision modular arithmetic */
z = ((x - y * sinDP1) - y * sinDP2) - y * sinDP3;

zz = z * z;
if( (j==1) || (j==2) )
	{
	y = One - ldexpl(zz,-1) + zz * zz * polevll( zz, coscof, 6 );
	}
else
	{
	y = z  +  z * (zz * polevll( zz, sincof, 6 ));
	}

if(sign < 0)
	y = -y;

return(y);
}





long double cosl(x)
long double x;
{
long double y, z, zz;
long i;
int j, sign;


/* make argument positive */
sign = 1;
if( x < 0 )
	x = -x;

if( x > sinlossth )
	{
	mtherr( "cosl", TLOSS );
	return(Zero);
	}

y = floorl( x/PIO4L );
z = ldexpl( y, -4 );
z = floorl(z);		/* integer part of y/8 */
z = y - ldexpl( z, 4 );  /* y - 16 * (y/16) */

/* integer and fractional part modulo one octant */
i = (int )z;
if( i & 1 )	/* map zeros to origin */
	{
	i += 1;
	y += One;
	}
j = i & 07;
if( j > 3)
	{
	j -=4;
	sign = -sign;
	}

if( j > 1 )
	sign = -sign;

/* Extended precision modular arithmetic */
z = ((x - y * sinDP1) - y * sinDP2) - y * sinDP3;

zz = z * z;
if( (j==1) || (j==2) )
	{
	y = z  +  z * (zz * polevll( zz, sincof, 6 ));
	}
else
	{
	y = One - ldexpl(zz,-1) + zz * zz * polevll( zz, coscof, 6 );
	}

if(sign < 0)
	y = -y;

return(y);
}
