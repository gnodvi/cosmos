/* One-term nutation program uses only the 18.6 year term.
 *
 *
 * References:
 * "Summary of 1980 IAU Theory of Nutation (Final Report of the
 * IAU Working Group on Nutation)", P. K. Seidelmann et al., in
 * Transactions of the IAU Vol. XVIII A, Reports on Astronomy,
 * P. A. Wayman, ed.; D. Reidel Pub. Co., 1982.
 *
 * "Nutation and the Earth's Rotation",
 * I.A.U. Symposium No. 78, May, 1977, page 256.
 * I.A.U., 1980.
 *
 * Woolard, E.W., "A redevelopment of the theory of nutation",
 * The Astronomical Journal, 58, 1-3 (1953).
 */

#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#include "ssystem.h"
#include "const.h"
#endif

#if LDOUBLE
#if IBMPC
/* 1980 series */
static short nuAi[] = {
0xd17f,0x30cc,0x4f41,0x94fb,0x3fec, XPD
0xc469,0x5b0e,0x9761,0x8834,0x3ff6, XPD
0x41eb,0x9bf7,0x5c3f,0xf1c4,0xc009, XPD
0x5656,0x21e4,0xcb9d,0xfa16,0x4005, XPD
0x14a5,0x7b74,0x6349,0x8eb4,0xbff9, XPD
0xb780,0x8240,0xc7e2,0x8998,0xc003, XPD
0x73c1,0xe1ef,0xe392,0xe94e,0x3ff4, XPD
0x3d71,0xd70a,0x70a3,0x933d,0x4002, XPD
};
short C3600i[] = {0x0000,0x0000,0x0000,0xe100,0x400a, XPD};
short C360i[] = {0x0000,0x0000,0x0000,0xb400,0x4007, XPD};
#define nuA ((DOUBLE *)nuAi)
#define C3600 (*(DOUBLE *)C3600i)
#define C360 (*(DOUBLE *)C360i)
#endif /* long double IBMPC */
#if MIEEE
/* 1980 series */
long nuAi[24] = {
0x3fec0000,0x94fb4f41,0x30ccd17f,
0x3ff60000,0x88349761,0x5b0ec469,
0xc0090000,0xf1c45c3f,0x9bf741eb,
0x40050000,0xfa16cb9d,0x21e45656,
0xbff90000,0x8eb46349,0x7b7414a5,
0xc0030000,0x8998c7e2,0x8240b780,
0x3ff40000,0xe94ee392,0xe1ef73c1,
0x40020000,0x933d70a3,0xd70a3d71,
};
long C3600i[] = {0x400a0000,0xe1000000,0x00000000};
long C360i[] = {0x40070000,0xb4000000,0x00000000};
#define nuA ((DOUBLE *)nuAi)
#define C3600 (*(DOUBLE *)C3600i)
#define C360 (*(DOUBLE *)C360i)
#endif /* long double MIEEE */
#if UNK
static long double nuA[8] = {
 2.22000000000000000000E-6L,
 2.07833000000000000000E-3L,
-1.93413626080000000000E3L,
 1.25044522200000000000E2L,
-1.74200000000000000000E-2L,
-1.71996000000000000000E1L,
 8.90000000000000000000E-4L,
 9.20250000000000000000E0L,
};
long double C3600 = 3600.0L;
long double C360 = 360.0L;
#endif
#else /* not LDOUBLE: */
/* 1977 series */
/*
DOUBLE nuA[8] = {
 8.00000000000000000000E-3,
 7.48000000000000000000E0,
-6.96291123000000000000E6,
 9.33059790000000000000E5,
-1.73700000000000000000E-2,
-1.71860000000000000000E1,
 9.10000000000000000000E-4,
 9.20159000000000000000E0,
};
*/
/* 1980 series */
#if UNK
static double nuA[8] = {
 2.22000000000000000000E-6,
 2.07833000000000000000E-3,
-1.93413626080000000000E3,
 1.25044522200000000000E2,
-1.74200000000000000000E-2,
-1.71996000000000000000E1,
 8.90000000000000000000E-4,
 9.20250000000000000000E0,
};
double C3600 = 3600.0;
double C360 = 360.0;
#endif /* double UNK */
#if DEC
static short nuAi[32] = {
0033424,0175517,0040460,0146321,
0036010,0032227,0060533,0007304,
0142761,0142134,0037633,0173502,
0041772,0013313,0116441,0162126,
0136616,0132143,0044573,0072025,
0141211,0114307,0161202,0040270,
0035551,0047343,0111341,0167564,
0041023,0036560,0121727,0005075,
};
short C3600i[] = {0043141,0000000,0000000,0000000};
short C360i[] = {0042264,0000000,0000000,0000000};
#define nuA ((DOUBLE *)nuAi)
#define C3600 (*(DOUBLE *)C3600i)
#define C360 (*(DOUBLE *)C360i)
#endif /* double DEC */
#if IBMPC
static short nuAi[40] = {
0x199a,0xe826,0x9f69,0x3ec2,
0x61d9,0xec2b,0x0692,0x3f61,
0x7ee8,0x87f3,0x388b,0xc09e,
0x3c8b,0x73a4,0x42d9,0x405f,
0x6e83,0x692f,0xd68c,0xbf91,
0x4817,0xfc50,0x3318,0xc031,
0x3dee,0x725c,0x29dc,0x3f4d,
0xe148,0x147a,0x67ae,0x4022,
};
short C3600i[] = {0x0000,0x0000,0x2000,0x40ac};
short C360i[] = {0x0000,0x0000,0x8000,0x4076};
#define nuA ((DOUBLE *)nuAi)
#define C3600 (*(DOUBLE *)C3600i)
#define C360 (*(DOUBLE *)C360i)
#endif /* double IBMPC */
#if MIEEE
static short nuAi[40] = {
0x3ec2,0x9f69,0xe826,0x199a,
0x3f61,0x0692,0xec2b,0x61d9,
0xc09e,0x388b,0x87f3,0x7ee8,
0x405f,0x42d9,0x73a4,0x3c8b,
0xbf91,0xd68c,0x692f,0x6e83,
0xc031,0x3318,0xfc50,0x4817,
0x3f4d,0x29dc,0x725c,0x3dee,
0x4022,0x67ae,0x147a,0xe148,
};
short C3600i[] = {0x40ac,0x2000,0x0000,0x0000};
short C360i[] = {0x4076,0x8000,0x0000,0x0000};
#define nuA ((DOUBLE *)nuAi)
#define C3600 (*(DOUBLE *)C3600i)
#define C360 (*(DOUBLE *)C360i)
#endif /* double MIEEE */
#endif /* not LDOUBLE */

/* The answers are posted here by nutlo():
 */
DOUBLE jdnut = 0.0;	/* time to which the nutation applies */
DOUBLE nutl = 0.0;	/* nutation in longitude (radians) */
DOUBLE nuto = 0.0;	/* nutation in obliquity (radians) */
extern DOUBLE eps, coseps, sineps, DTR, STR, J1900, J2000, Jcentury;
DOUBLE FLOOR();
#define mod360(x) (x - C360*FLOOR(x/C360))

nutlo(J)
DOUBLE J;
{
DOUBLE T, OM, C, D;
DOUBLE SIN(), COS();

if( jdnut == J )
	return(0);
jdnut = J;

/* 1977 IAU series */
/*
 * T = (J-J1900)/Jcentury;
 * OM = (((nuA[0]*T + nuA[1])*T + nuA[2])*T + nuA[3])/C3600;
 */

/* 1980 IAU series */
T = (J-J2000)/Jcentury;
OM = ((nuA[0]*T + nuA[1])*T + nuA[2])*T + nuA[3];
OM = DTR * mod360(OM);

/* Woolard series
 * C = (-0.01737*T - 17.2327)*SIN(OM);
 * D = ( 0.00091*T +  9.2100)*COS(OM);
 */

/*
 * C = (nuA[4]*T + nuA[5])*SIN(OM);
 * D = ( nuA[6]*T + nuA[7])*COS(OM);
 */

C = nuA[5]*SIN(OM);
D = nuA[7]*COS(OM);

/* Save answers, expressed in radians */
nutl = STR * C;
nuto = STR * D;
return(1);
}



/* Nutation -- AA page B20
 * using nutation in longitude and obliquity from nutlo()
 * and obliquity of the ecliptic from epsiln()
 * both calculated for Julian date J.
 *
 * p[] = equatorial rectangular position vector of object for
 * mean ecliptic and equinox of date.
 */
DOUBLE nuN00, nuN01, nuN02, nuN10, nuN11, nuN12, nuN20, nuN21, nuN22;

nutate( J, p )
DOUBLE J;
DOUBLE p[];
{
DOUBLE ce, se, cl, sl, sino, f;
DOUBLE p1[3];
int i;
DOUBLE SIN(), COS();

if( J == jdnut )
	goto matmpy;
nutlo(J); /* be sure we calculated nutl and nuto */
epsiln(J); /* and also the obliquity of date */

f = eps + nuto;
ce = COS( f );
se = SIN( f );
sino = SIN(nuto);
cl = COS( nutl );
sl = SIN( nutl );

/* Apply adjustment
 * to equatorial rectangular coordinates of object.
 *
 * This is a composite of three rotations: rotate about x axis
 * to ecliptic of date; rotate about new z axis by the nutation
 * in longitude; rotate about new x axis back to equator of date
 * plus nutation in obliquity.
 */
nuN00 =   cl;
nuN01 =  -sl*coseps;
nuN02 =  -sl*sineps;

nuN10 = sl*ce;
nuN11 = cl*coseps*ce + sineps*se;
nuN12 = -( sino + (One-cl)*sineps*se );

nuN20 = sl*se;
nuN21 = sino + (cl-One)*se*coseps;
nuN22 = cl*sineps*se + coseps*ce;

matmpy:

p1[0] = nuN00*p[0] + nuN01*p[1] + nuN02*p[2];
p1[1] = nuN10*p[0] + nuN11*p[1] + nuN12*p[2];
p1[2] = nuN20*p[0] + nuN21*p[1] + nuN22*p[2];

for( i=0; i<3; i++ )
	p[i] = p1[i];
}


/* Inverse transformation using previously
 * computed matrix elements.
 */
invnut( p1 )
DOUBLE p1[];
{
DOUBLE p[3];
int i;

p[0] = nuN00*p1[0] + nuN10*p1[1] + nuN20*p1[2];
p[1] = nuN01*p1[0] + nuN11*p1[1] + nuN21*p1[2];
p[2] = nuN02*p1[0] + nuN12*p1[1] + nuN22*p1[2];

for( i=0; i<3; i++ )
	p1[i] = p[i];
}
