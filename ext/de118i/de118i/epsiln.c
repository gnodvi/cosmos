/* Obliquity of the ecliptic at Julian date J
 *
 * IAU Coefficients are from:
 * J. H. Lieske, T. Lederle, W. Fricke, and B. Morando,
 * "Expressions for the Precession Quantities Based upon the IAU
 * (1976) System of Astronomical Constants,"  Astronomy and Astrophysics
 * 58, 1-16 (1977).
 *
 * Before or after 200 years from J2000, the formula used is from:
 * J. Laskar, "Secular terms of classical planetary theories
 * using the results of general theory," Astronomy and Astrophysics
 * 157, 59070 (1986).
 *
 *  See precess.c and page B18 of the Astronomical Almanac.
 * 
 */

#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#include "const.h"
#endif

#if LDOUBLE
#if IBMPC
short epsAi[] = {
0xb229,0x50d6,0x2f6a,0xeda2,0x3ff5, XPD
0xedd0,0x8d25,0x3ad1,0x9aaa,0xbff4, XPD
0xc28f,0x28f5,0x8f5c,0xbb42,0xc004, XPD
0x4dd3,0x1062,0xb958,0xa4ce,0x400f, XPD
};
short epsNi[] = {
0xfbf4,0xcdfe,0x138b,0xed5f,0x3ff5, XPD
0xedd0,0x8d25,0x3ad1,0x9aaa,0xbff4, XPD
0xe148,0x147a,0x47ae,0xbb61,0xc004, XPD
0xa4a9,0x404e,0x2113,0xa4e6,0x400f, XPD
};
short epsLi[] = {
0x12eb,0x0788,0xaf45,0x86b0,0x3fdf, XPD
0xf75f,0xd639,0x60eb,0xc6f1,0x3fe3, XPD
0x43e0,0x4bdb,0x3c80,0x95a0,0x3fe9, XPD
0x3646,0xb013,0x4476,0xbf20,0x3fea, XPD
0x7f66,0x2345,0x9e44,0xa3c9,0xbff0, XPD
0xad85,0x117e,0xacd9,0xa39f,0xbff6, XPD
0xde1a,0xc1ac,0xaafb,0xa85c,0xbff7, XPD
0x8106,0x4395,0x6c8b,0xffe7,0x3fff, XPD
0xc083,0xa1ca,0xb645,0xfdf3,0xbff8, XPD
0x9581,0x8b43,0xe76c,0xea0b,0xc007, XPD
0x4dd3,0x1062,0xb958,0xa4ce,0x400f, XPD
};
#define epsA ((long double *)epsAi)
#define epsN ((long double *)epsNi)
#define epsL ((long double *)epsLi)
#endif
#if MIEEE
long epsAi[] = {
0x3ff50000,0xeda22f6a,0x50d6b229,
0xbff40000,0x9aaa3ad1,0x8d25edd0,
0xc0040000,0xbb428f5c,0x28f5c28f,
0x400f0000,0xa4ceb958,0x10624dd3,
};
long epsNi[] = {
0x3ff50000,0xed5f138b,0xcdfefbf4,
0xbff40000,0x9aaa3ad1,0x8d25edd0,
0xc0040000,0xbb6147ae,0x147ae148,
0x400f0000,0xa4e62113,0x404ea4a9,
};
long epsLi[] = {
0x3fdf0000,0x86b0af45,0x078812eb,
0x3fe30000,0xc6f160eb,0xd639f75f,
0x3fe90000,0x95a03c80,0x4bdb43e0,
0x3fea0000,0xbf204476,0xb0133646,
0xbff00000,0xa3c99e44,0x23457f66,
0xbff60000,0xa39facd9,0x117ead85,
0xbff70000,0xa85caafb,0xc1acde1a,
0x3fff0000,0xffe76c8b,0x43958106,
0xbff80000,0xfdf3b645,0xa1cac083,
0xc0070000,0xea0be76c,0x8b439581,
0x400f0000,0xa4ceb958,0x10624dd3,
};
#define epsA ((long double *)epsAi)
#define epsN ((long double *)epsNi)
#define epsL ((long double *)epsLi)
#endif
#if UNK
long double epsA[] = {
 1.81300000000000000000E-3L,
-5.90000000000000000000E-4L,
-4.68150000000000000000E1L,
 8.43814480000000000000E4L,
};
long double epsN[] = {
 1.81100000000000000000E-3L,
-5.90000000000000000000E-4L,
-4.68450000000000000000E1L,
 8.44282584000000000000E4L,
};
long double epsL[] = {
 2.45000000000000000000E-10L,
 5.79000000000000000000E-9L,
 2.78700000000000000000E-7L,
 7.12000000000000000000E-7L,
-3.90500000000000000000E-5L,
-2.49670000000000000000E-3L,
-5.13800000000000000000E-3L,
 1.99925000000000000000E0L,
-1.55000000000000000000E-2L,
-4.68093000000000000000E2L,
 8.43814480000000000000E4L,
};
#endif
#else /* not LDOUBLE: */
#if DEC
static short epsAi[4*4] = {
0035755,0121057,0065120,0153262,
0135432,0125072,0150615,0022756,
0141473,0041217,0056050,0172703,
0044244,0147271,0054020,0061116,
};
static short epsNi[4*4] = {
0035755,0057423,0105715,0177374,
0135432,0125072,0150615,0022756,
0141473,0060507,0127024,0075341,
0044244,0163041,0011500,0047245,
};
static short epsLi[11*4] = {
0030206,0130257,0042407,0104023,
0031306,0170540,0165726,0034767,
0032625,0120074,0100113,0155504,
0033077,0020104,0073260,0011466,
0134443,0144636,0042043,0042577,
0136043,0117654,0154421,0077256,
0136250,0056252,0175701,0126336,
0040377,0163554,0105503,0112601,
0136575,0171666,0042641,0145301,
0142352,0005747,0066213,0041626,
0044244,0147271,0054020,0061116,
};
#endif /* double DEC */
#if IBMPC
static short epsAi[4*4] = {
0x1ad6,0xed4a,0xb445,0x3f5d,
0xa4be,0x5a31,0x5547,0xbf43,
0x1eb8,0xeb85,0x6851,0xc047,
0x0c4a,0x2b02,0x99d7,0x40f4,
};
static short epsNi[4*4] = {
0xbfdf,0x7179,0xabe2,0x3f5d,
0xa4be,0x5a31,0x5547,0xbf43,
0x8f5c,0xf5c2,0x6c28,0xc047,
0x09d5,0x2268,0x9cc4,0x40f4,
};
static short epsLi[11*4] = {
0xf102,0xe8a0,0xd615,0x3df0,
0xc73f,0x1d7a,0xde2c,0x3e38,
0x7b68,0x9009,0xb407,0x3e92,
0x0267,0x8ed6,0xe408,0x3ea7,
0x68b0,0xc884,0x7933,0xbf04,
0x2fd6,0x9b22,0x73f5,0xbf64,
0x359c,0x5f78,0x0b95,0xbf75,
0x72b0,0x9168,0xfced,0x3fff,
0x3958,0xc8b4,0xbe76,0xbf8f,
0x6873,0xed91,0x417c,0xc07d,
0x0c4a,0x2b02,0x99d7,0x40f4,
};
#endif /* double IBMPC */
#if MIEEE
static short epsAi[4*4] = {
0x3f5d,0xb445,0xed4a,0x1ad6,
0xbf43,0x5547,0x5a31,0xa4be,
0xc047,0x6851,0xeb85,0x1eb8,
0x40f4,0x99d7,0x2b02,0x0c4a,
};
static short epsNi[4*4] = {
0x3f5d,0xabe2,0x7179,0xbfdf,
0xbf43,0x5547,0x5a31,0xa4be,
0xc047,0x6c28,0xf5c2,0x8f5c,
0x40f4,0x9cc4,0x2268,0x09d5,
};
static short epsLi[11*4] = {
0x3df0,0xd615,0xe8a0,0xf102,
0x3e38,0xde2c,0x1d7a,0xc73f,
0x3e92,0xb407,0x9009,0x7b68,
0x3ea7,0xe408,0x8ed6,0x0267,
0xbf04,0x7933,0xc884,0x68b0,
0xbf64,0x73f5,0x9b22,0x2fd6,
0xbf75,0x0b95,0x5f78,0x359c,
0x3fff,0xfced,0x9168,0x72b0,
0xbf8f,0xbe76,0xc8b4,0x3958,
0xc07d,0x417c,0xed91,0x6873,
0x40f4,0x99d7,0x2b02,0x0c4a,
};
#endif /* double MIEEE */
#if UNK
DOUBLE epsA[] = {
 1.81300000000000000000E-3,
-5.90000000000000000000E-4,
-4.68150000000000000000E1,
 8.43814480000000000000E4,
};
DOUBLE epsN[] = {
 1.81100000000000000000E-3,
-5.90000000000000000000E-4,
-4.68450000000000000000E1,
 8.44282584000000000000E4,
};
DOUBLE epsL[] = {
 2.45000000000000000000E-10,
 5.79000000000000000000E-9,
 2.78700000000000000000E-7,
 7.12000000000000000000E-7,
-3.90500000000000000000E-5,
-2.49670000000000000000E-3,
-5.13800000000000000000E-3,
 1.99925000000000000000E0,
-1.55000000000000000000E-2,
-4.68093000000000000000E2,
 8.43814480000000000000E4,
};
#else /* not UNK: */
#define epsA ((double *)epsAi)
#define epsN ((double *)epsNi)
#define epsL ((double *)epsLi)
#endif /* UNK */
#endif /* not LDOUBLE */

/* The results of the program are returned in these
 * global variables:
 */
DOUBLE jdeps = 0.0; /* Date for which obliquity was last computed */
DOUBLE eps = 0.0; /* The computed obliquity in radians */
DOUBLE coseps = 0.0; /* Cosine of the obliquity */
DOUBLE sineps = 0.0; /* Sine of the obliquity */
extern DOUBLE eps, coseps, sineps, STR, J1900, J2000, Jcentury, Jmillenium;
DOUBLE SIN(), COS(), FABS();

epsiln(J)
DOUBLE J; /* Julian date input */
{
DOUBLE T;

if( J == jdeps )
	goto done;

#if 0
T = (J - J2000)/Jcentury;

/* This expansion is from the AA.
 * Note the official 1976 IAU number is 23d 26' 21.448", but
 * the JPL numerical integration found 21.4119".
 */
if( FABS(T) < Two )
	eps = (((epsA[0]*T + epsA[1])*T +epsA[2])*T + epsA[3])*STR;
#endif

#if 0
T = (J - J1900)/Jcentury;
if( FABS(T) < Two )
	eps = (((epsN[0]*T + epsN[1])*T + epsN[2])*T + epsN[3])*STR;
/* This expansion is from Laskar, cited above.
 * Bretagnon and Simon say, in Planetary Programs and Tables, that it
 * is accurate to 0.1" over a span of 6000 years. Laskar estimates the
 * precision to be 0.01" after 1000 years and a few seconds of arc
 * after 10000 years.
 */
else
#endif

#if 1
	{
	T = (J - J2000)/Jmillenium;
	eps = ((((((((( epsL[0]*T + epsL[1])*T + epsL[2])*T
	+ epsL[3])*T + epsL[4])*T + epsL[5])*T + epsL[6])*T
	+ epsL[7])*T + epsL[8])*T + epsL[9])*T + epsL[10];
	eps *= STR;
	}
#endif

coseps = COS( eps );
sineps = SIN( eps );
jdeps = J;
done: ;
}
