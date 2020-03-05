/* Precession of the equinox and ecliptic
 * from epoch Julian date J to or from J2000.0
 *
 * Program by Steve Moshier.
 *
 *
 * The IAU formula is not used in this version.
 * IAU Coefficients are from:
 * J. H. Lieske, T. Lederle, W. Fricke, and B. Morando,
 * "Expressions for the Precession Quantities Based upon the IAU
 * (1976) System of Astronomical Constants,"  Astronomy and
 * Astrophysics 58, 1-16 (1977).
 *
 * Newer formulas that cover a much longer time span are from:
 * J. Laskar, "Secular terms of classical planetary theories
 * using the results of general theory," Astronomy and Astrophysics
 * 157, 59070 (1986).
 *
 * See also:
 * P. Bretagnon and G. Francou, "Planetary theories in rectangular
 * and spherical variables. VSOP87 solutions," Astronomy and
 * Astrophysics 202, 309-315 (1988).
 *
 * Laskar's expansions are said by Bretagnon and Francou
 * to have "a precision of about 1" over 10000 years before
 * and after J2000.0 in so far as the precession constants p^0_A
 * and epsilon^0_A are perfectly known."
 *
 * Bretagnon and Francou's expansions for the node and inclination
 * of the ecliptic were derived from Laskar's data but were truncated
 * after the term in T**6. I have recomputed these expansions from
 * Laskar's data, retaining powers up to T**10 in the result.
 *
 * The following table indicates the differences between the result
 * of the IAU formula and Laskar's formula using four different test
 * vectors, checking at J2000 plus and minus the indicated number
 * of years.
 *
 *   Years       Arc
 * from J2000  Seconds
 * ----------  -------
 *        0	  0
 *      100	.006	
 *      200     .006
 *      500     .015
 *     1000     .28
 *     2000    6.4
 *     3000   38.
 *    10000 9400.
 */

#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#endif
#if LDOUBLE
#if UNK
/* Precession coefficients taken from Laskar's paper: */
static DOUBLE pAcof[10] = {
-8.66000000000000000000E-10L,
-4.75900000000000000000E-8L,
 2.42400000000000000000E-7L,
 1.30950000000000000000E-5L,
 1.74510000000000000000E-4L,
-1.80550000000000000000E-3L,
-2.35316000000000000000E-1L,
 7.73200000000000000000E-2L,
 1.11197100000000000000E2L,
 5.02909660000000000000E4L
};
/* Node and inclination of the earth's orbit computed from
 * Laskar's data as done in Bretagnon and Francou's paper:
 */
static DOUBLE nodecof[11] = {
 6.64020000000000000000E-16L,
-2.69151000000000000000E-15L,
-1.54702100000000000000E-12L,
 7.52131300000000000000E-12L,
 6.31901310000000000000E-10L,
-3.48388152000000000000E-9L,
-1.81306589600000000000E-7L,
 2.75036225000000000000E-8L,
 7.43945314260000000000E-5L,
-4.20786043170000000000E-2L,
 3.05211265497500000000E0L
};
static DOUBLE inclcof[11] = {
 1.21470000000000000000E-16L,
 7.37590000000000000000E-17L,
-8.26287000000000000000E-14L,
 2.50341000000000000000E-13L,
 2.46508390000000000000E-11L,
-5.40004410000000000000E-11L,
 1.32115526000000000000E-9L,
-5.99873702700000000000E-7L,
-1.62427970910000000000E-5L,
 2.27849553700000000000E-3L,
 0.00000000000000000000E0L
};
#endif
#if IBMPC
static short pAcof[] = {
0x3761,0xf547,0x551b,0xee0b,0xbfe0, XPD
0x6a6d,0x43d6,0xc224,0xcc65,0xbfe6, XPD
0x3ffe,0x5966,0x33cb,0x8223,0x3fe9, XPD
0x9413,0x06aa,0x98c4,0xdbb2,0x3fee, XPD
0xdda4,0x9c6c,0xabe2,0xb6fc,0x3ff2, XPD
0xc6e3,0xe62d,0x86e7,0xeca6,0xbff5, XPD
0xe8c0,0xe6f2,0xad70,0xf0f6,0xbffc, XPD
0x6018,0x9d1f,0xf2ba,0x9e59,0x3ffb, XPD
0x4c98,0x8c15,0xea4a,0xde64,0x4005, XPD
0xef9e,0xc6a7,0xf74b,0xc472,0x400e, XPD
};
static short nodecof[] = {
0x38c8,0xf647,0x072a,0xbf64,0x3fcc, XPD
0x9b76,0xc03f,0x989c,0xc1f1,0xbfce, XPD
0x88e7,0x5905,0x4e3b,0xd9b9,0xbfd7, XPD
0xa57c,0x25f2,0xfb80,0x8450,0x3fda, XPD
0xf7b9,0x5a5a,0x1a04,0xadb2,0x3fe0, XPD
0x41a8,0xe90f,0x1783,0xef69,0xbfe2, XPD
0x8c7c,0x7312,0x2d05,0xc2ad,0xbfe8, XPD
0xe099,0x5ad6,0x1b01,0xec41,0x3fe5, XPD
0x8ff0,0x1112,0x428b,0x9c04,0x3ff1, XPD
0xbe8d,0x7207,0x9d56,0xac5a,0xbffa, XPD
0xa4e6,0x34d2,0xd051,0xc355,0x4000, XPD
};
static short inclcof[] = {
0x0681,0xeffb,0x9db4,0x8c0b,0x3fca, XPD
0xf1e2,0xed33,0xa0f0,0xaa13,0x3fc9, XPD
0xcd60,0x3987,0x33db,0xba10,0xbfd3, XPD
0xb113,0x604a,0xf0b7,0x8ced,0x3fd5, XPD
0x206c,0xe1ba,0xc131,0xd8d4,0x3fdb, XPD
0xb3b2,0xfa51,0x176b,0xed7f,0xbfdc, XPD
0xd03a,0x5b5d,0x04ac,0xb594,0x3fe1, XPD
0xe322,0xf2f6,0x01c7,0xa107,0xbfea, XPD
0x7f4e,0x73db,0x2422,0x8841,0xbfef, XPD
0x9c64,0xc468,0xcfd0,0x9552,0x3ff6, XPD
0x0000,0x0000,0x0000,0x0000,0x0000, XPD
};
#endif
#if MIEEE
static long pAcof[30] = {
0xbfe00000,0xee0b551b,0xf5473761,
0xbfe60000,0xcc65c224,0x43d66a6d,
0x3fe90000,0x822333cb,0x59663ffe,
0x3fee0000,0xdbb298c4,0x06aa9413,
0x3ff20000,0xb6fcabe2,0x9c6cdda4,
0xbff50000,0xeca686e7,0xe62dc6e3,
0xbffc0000,0xf0f6ad70,0xe6f2e8c0,
0x3ffb0000,0x9e59f2ba,0x9d1f6018,
0x40050000,0xde64ea4a,0x8c154c98,
0x400e0000,0xc472f74b,0xc6a7ef9e
};
static long nodecof[33] = {
0x3fcc0000,0xbf64072a,0xf64738c8,
0xbfce0000,0xc1f1989c,0xc03f9b76,
0xbfd70000,0xd9b94e3b,0x590588e7,
0x3fda0000,0x8450fb80,0x25f2a57c,
0x3fe00000,0xadb21a04,0x5a5af7b9,
0xbfe20000,0xef691783,0xe90f41a8,
0xbfe80000,0xc2ad2d05,0x73128c7c,
0x3fe50000,0xec411b01,0x5ad6e099,
0x3ff10000,0x9c04428b,0x11128ff0,
0xbffa0000,0xac5a9d56,0x7207be8d,
0x40000000,0xc355d051,0x34d2a4e6
};
static long inclcof[33] = {
0x3fca0000,0x8c0b9db4,0xeffb0681,
0x3fc90000,0xaa13a0f0,0xed33f1e2,
0xbfd30000,0xba1033db,0x3987cd60,
0x3fd50000,0x8cedf0b7,0x604ab113,
0x3fdb0000,0xd8d4c131,0xe1ba206c,
0xbfdc0000,0xed7f176b,0xfa51b3b2,
0x3fe10000,0xb59404ac,0x5b5dd03a,
0xbfea0000,0xa10701c7,0xf2f6e322,
0xbfef0000,0x88412422,0x73db7f4e,
0x3ff60000,0x9552cfd0,0xc4689c64,
0x00000000,0x00000000,0x00000000
};
#endif

#else /* regular double */
#if UNK
static DOUBLE pAcof[10] = {
-8.66000000000000000000E-10,
-4.75900000000000000000E-8,
 2.42400000000000000000E-7,
 1.30950000000000000000E-5,
 1.74510000000000000000E-4,
-1.80550000000000000000E-3,
-2.35316000000000000000E-1,
 7.73200000000000000000E-2,
 1.11197100000000000000E2,
 5.02909660000000000000E4
};
static DOUBLE nodecof[11] = {
 6.64020000000000000000E-16,
-2.69151000000000000000E-15,
-1.54702100000000000000E-12,
 7.52131300000000000000E-12,
 6.31901310000000000000E-10,
-3.48388152000000000000E-9,
-1.81306589600000000000E-7,
 2.75036225000000000000E-8,
 7.43945314260000000000E-5,
-4.20786043170000000000E-2,
 3.05211265497500000000E0
};
static DOUBLE inclcof[11] = {
 1.21470000000000000000E-16,
 7.37590000000000000000E-17,
-8.26287000000000000000E-14,
 2.50341000000000000000E-13,
 2.46508390000000000000E-11,
-5.40004410000000000000E-11,
 1.32115526000000000000E-9,
-5.99873702700000000000E-7,
-1.62427970910000000000E-5,
 2.27849553700000000000E-3,
 0.00000000000000000000E0
};
#endif
#if DEC
static short pAcof[40] = {
0130556,0005525,0015765,0043467,
0132114,0062702,0022103,0153152,
0032602,0021463,0145531,0063100,
0034133,0131230,0142006,0125224,
0035066,0176253,0161234,0066336,
0135754,0123206,0163746,0026707,
0137560,0173255,0070346,0171351,
0037236,0054762,0135235,0017540,
0041736,0062352,0045214,0012515,
0044104,0071367,0045706,0123760
};
static short nodecof[44] = {
0023477,0062007,0025366,0043471,
0124101,0170630,0116300,0037633,
0126331,0134516,0035531,0002611,
0027004,0050373,0100045,0171245,
0030455,0131032,0002132,0055370,
0131157,0064427,0101751,0007502,
0132502,0126455,0002563,0011214,
0031754,0040433,0000532,0153341,
0034634,0002102,0105421,0011220,
0137054,0055235,0053162,0003677,
0040503,0052720,0050464,0151245
};
static short inclcof[44] = {
0023014,0005635,0132357,0175406,
0022652,0011640,0170355,0031762,
0125272,0010063,0155471,0103715,
0025614,0166760,0133540,0045261,
0027330,0152301,0030741,0135040,
0127555,0077427,0065772,0050664,
0030665,0112004,0126133,0056720,
0133041,0003401,0143762,0173343,
0134210,0040444,0021163,0155577,
0036025,0051317,0150304,0064234,
0000000,0000000,0000000,0000000,
};
#endif
#if IBMPC
static short pAcof[40] = {
0xa8e7,0xa37e,0xc16a,0xbe0d,
0x7acd,0x4488,0x8cb8,0xbe69,
0x2cc8,0x796b,0x4466,0x3e90,
0xd553,0x1880,0x7653,0x3eeb,
0x8d9c,0x7c53,0xdf95,0x3f26,
0xc5b9,0xdcfc,0x94d0,0xbf5d,
0xde5d,0xae1c,0x1ed5,0xbfce,
0xa3ec,0x5753,0xcb3e,0x3fb3,
0x82aa,0x4951,0xcc9d,0x405b,
0xd4fe,0xe978,0x8e5e,0x40e8
};
static short nodecof[44] = {
0xc8e7,0xe55e,0xec80,0x3cc7,
0x07f3,0x1398,0x3e33,0xbce8,
0x20b1,0xc76b,0x3729,0xbd7b,
0xbe55,0x7004,0x8a1f,0x3da0,
0x4b5f,0x408b,0xb643,0x3e05,
0x21e8,0xf07d,0xed22,0xbe2d,
0x6252,0xa0ae,0x55a5,0xbe88,
0x5adc,0x602b,0x8823,0x3e5d,
0x2252,0x5162,0x8088,0x3f13,
0x40f8,0xaace,0x8b53,0xbfa5,
0x9a55,0x0a26,0x6aba,0x4008
};
static short inclcof[44] = {
0xff61,0xb69d,0x8173,0x3ca1,
0xa67e,0x1e1d,0x4274,0x3c95,
0x30fa,0x7b67,0x4206,0xbd37,
0x0956,0x16ec,0x9dbe,0x3d51,
0x3744,0x263c,0x1a98,0x3dbb,
0x4a36,0xed7f,0xafe2,0xbdcd,
0x6bba,0x958b,0xb280,0x3e16,
0x5edc,0x38fe,0x20e0,0xbea4,
0x7b70,0x844e,0x0824,0xbef1,
0x8d14,0xfa18,0xaa59,0x3f62,
0x0000,0x0000,0x0000,0x0000,
};
#endif
#if MIEEE
static short pAcof[40] = {
0xbe0d,0xc16a,0xa37e,0xa8e7,
0xbe69,0x8cb8,0x4488,0x7acd,
0x3e90,0x4466,0x796b,0x2cc8,
0x3eeb,0x7653,0x1880,0xd553,
0x3f26,0xdf95,0x7c53,0x8d9c,
0xbf5d,0x94d0,0xdcfc,0xc5b9,
0xbfce,0x1ed5,0xae1c,0xde5d,
0x3fb3,0xcb3e,0x5753,0xa3ec,
0x405b,0xcc9d,0x4951,0x82aa,
0x40e8,0x8e5e,0xe978,0xd4fe
};
static short nodecof[44] = {
0x3cc7,0xec80,0xe55e,0xc8e7,
0xbce8,0x3e33,0x1398,0x07f3,
0xbd7b,0x3729,0xc76b,0x20b1,
0x3da0,0x8a1f,0x7004,0xbe55,
0x3e05,0xb643,0x408b,0x4b5f,
0xbe2d,0xed22,0xf07d,0x21e8,
0xbe88,0x55a5,0xa0ae,0x6252,
0x3e5d,0x8823,0x602b,0x5adc,
0x3f13,0x8088,0x5162,0x2252,
0xbfa5,0x8b53,0xaace,0x40f8,
0x4008,0x6aba,0x0a26,0x9a55
};
static short inclcof[44] = {
0x3ca1,0x8173,0xb69d,0xff61,
0x3c95,0x4274,0x1e1d,0xa67e,
0xbd37,0x4206,0x7b67,0x30fa,
0x3d51,0x9dbe,0x16ec,0x0956,
0x3dbb,0x1a98,0x263c,0x3744,
0xbdcd,0xafe2,0xed7f,0x4a36,
0x3e16,0xb280,0x958b,0x6bba,
0xbea4,0x20e0,0x38fe,0x5edc,
0xbef1,0x0824,0x844e,0x7b70,
0x3f62,0xaa59,0xfa18,0x8d14,
0x0000,0x0000,0x0000,0x0000
};
#endif

#endif

extern DOUBLE J2000; /* = 2451545.0, 2000 January 1.5 */
extern DOUBLE B1950;
extern DOUBLE STR; /* = 4.8481368110953599359e-6 radians per arc second */
extern DOUBLE Jmillenium; /* 365250.0 days */
extern DOUBLE coseps, sineps; /* see epsiln.c */

/* Subroutine arguments:
 *
 * R = rectangular equatorial coordinate vector to be precessed.
 *     The result is written back into the input vector.
 * J = Julian date
 * direction =
 *      Precess from J to J2000: direction = 1
 *      Precess from J2000 to J: direction = -1
 * Note that if you want to precess from J1 to J2, you would
 * first go from J1 to J2000, then call the program again
 * to go from J2000 to J2.
 */

precess( R, J, direction )
DOUBLE R[], J;
int direction;
{
DOUBLE A, B, T, z, pA, W;
DOUBLE x[3];
DOUBLE *p;
int i;
DOUBLE SIN(), COS(), FABS();

if( J == B1950 )
	goto done;

x[0] = R[0];
x[1] = R[1];
x[2] = R[2];

/* if going from B1950 to date, first rotate to J2000. */
if( direction == -1 )
	r118_200(x);

T = (J - J2000)/Jmillenium;

/* Implementation by elementary rotations using Laskar's expansions.
 * First rotate about the x axis from the initial equator
 * to the ecliptic. (The input is equatorial.)
 */
if( direction == 1 )
	epsiln( J ); /* To J2000 */
else
	epsiln( J2000 ); /* From J2000 */
z = coseps*x[1] + sineps*x[2];
x[2] = -sineps*x[1] + coseps*x[2];
x[1] = z;

/* Precession in longitude
 */
/* T /= 10.0; */ /* thousands of years */
p = (DOUBLE *)pAcof;
pA = *p++;
for( i=0; i<9; i++ )
	pA = pA * T + *p++;
pA *= STR * T;

/* Node of the moving ecliptic on the J2000 ecliptic.
 */
p = (DOUBLE *)nodecof;
W = *p++;
for( i=0; i<10; i++ )
	W = W * T + *p++;

/* Rotate about z axis to the node.
 */
if( direction == 1 )
	z = W + pA;
else
	z = W;
B = COS(z);
A = SIN(z);
z = B * x[0] + A * x[1];
x[1] = -A * x[0] + B * x[1];
x[0] = z;

/* Rotate about new x axis by the inclination of the moving
 * ecliptic on the J2000 ecliptic.
 */
p = (DOUBLE *)inclcof;
z = *p++;
for( i=0; i<10; i++ )
	z = z * T + *p++;
if( direction == 1 )
	z = -z;
B = COS(z);
A = SIN(z);
z = B * x[1] + A * x[2];
x[2] = -A * x[1] + B * x[2];
x[1] = z;

/* Rotate about new z axis back from the node.
 */
if( direction == 1 )
	z = -W;
else
	z = -W - pA;
B = COS(z);
A = SIN(z);
z = B * x[0] + A * x[1];
x[1] = -A * x[0] + B * x[1];
x[0] = z;

/* Rotate about x axis to final equator.
 */
if( direction == 1 )
	epsiln( J2000 );
else
	epsiln( J );
z = coseps * x[1] - sineps * x[2];
x[2] = sineps * x[1] + coseps * x[2];
x[1] = z;

/* if going from date to 1950, rotate J2000 answer to 1950. */
if( direction == 1 )
	r200_118(x);

for( i=0; i<3; i++ )
	R[i] = x[i];
done: ;
}
