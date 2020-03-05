/* oblate.c
 */

#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#include "ssystem.h"
#include "const.h"
#endif

#ifdef DEBUG
#undef DEBUG
#endif
#define DEBUG 0

extern DOUBLE GMs[];
extern DOUBLE Rij[NTOTAL][NTOTAL];
extern DOUBLE Rij3[NTOTAL][NTOTAL];

extern DOUBLE AM, AE, LOVENO, PHASE, EMRAT;
extern DOUBLE PSLP1, PSLPI, PSLPIA, PSLPIB, JDEPOCH;
extern DOUBLE LBET, LGAM, LGAMBET;
extern DOUBLE Jm[], Je[], Cnm[], Snm[];

DOUBLE ofdate[3];
DOUBLE phi, th, psi, phi1, th1, psi1, Nx, Ny, Nz;
DOUBLE sinpsi, cospsi, sinth, costh, sinphi, cosphi;
DOUBLE X00, X01, X02, X10, X11, X12, X20, X21, X22;
extern DOUBLE B1950;
DOUBLE SQRT(), SIN(), COS(), FLOOR();

#if LDOUBLE
#if IBMPC
short Cr375i[]  = {0x0000,0x0000,0x0000,0xc000,0x3ffd, XPD};
short C2r5i[]   = {0x0000,0x0000,0x0000,0xa000,0x4000, XPD};
short C3r75i[]  = {0x0000,0x0000,0x0000,0xf000,0x4000, XPD};
short C4r375i[] = {0x0000,0x0000,0x0000,0x8c00,0x4001, XPD};
short C7r5i[]   = {0x0000,0x0000,0x0000,0xf000,0x4001, XPD};
short C15i[]    = {0x0000,0x0000,0x0000,0xf000,0x4002, XPD};
short C17r5i[]  = {0x0000,0x0000,0x0000,0x8c00,0x4003, XPD};
short C52r5i[]  = {0x0000,0x0000,0x0000,0xd200,0x4004, XPD};
short C105i[]   = {0x0000,0x0000,0x0000,0xd200,0x4005, XPD};
#define Cr375 (*(long double *)Cr375i)
#define C2r5 (*(long double *)C2r5i)
#define C3r75 (*(long double *)C3r75i)
#define C4r375 (*(long double *)C4r375i)
#define C7r5 (*(long double *)C7r5i)
#define C15 (*(long double *)C15i)
#define C17r5 (*(long double *)C17r5i)
#define C52r5 (*(long double *)C52r5i)
#define C105 (*(long double *)C105i)
#endif
#if MIEEE
long Cr375i[]   = {0x3ffd0000,0xc0000000,0x00000000};
long C2r5i[]    = {0x40000000,0xa0000000,0x00000000};
long C3r75i[]   = {0x40000000,0xf0000000,0x00000000};
long C4r375i[]  = {0x40010000,0x8c000000,0x00000000};
long C7r5i[]    = {0x40010000,0xf0000000,0x00000000};
long C15i[]     = {0x40020000,0xf0000000,0x00000000};
long C17r5i[]   = {0x40030000,0x8c000000,0x00000000};
long C52r5i[]   = {0x40040000,0xd2000000,0x00000000};
long C105i[]    = {0x40050000,0xd2000000,0x00000000};
#define Cr375 (*(long double *)Cr375i)
#define C2r5 (*(long double *)C2r5i)
#define C3r75 (*(long double *)C3r75i)
#define C4r375 (*(long double *)C4r375i)
#define C7r5 (*(long double *)C7r5i)
#define C15 (*(long double *)C15i)
#define C17r5 (*(long double *)C17r5i)
#define C52r5 (*(long double *)C52r5i)
#define C105 (*(long double *)C105i)
#endif
#if UNK
long double Cr375 =  3.75000000000000000000E-1L;
long double C2r5 =  2.50000000000000000000E0L;
long double C3r75 =  3.75000000000000000000E0L;
long double C4r375 = 4.37500000000000000000E0L;
long double C7r5 = 7.50000000000000000000E0L;
long double C15 = 1.50000000000000000000E1L;
long double C17r5 = 1.75000000000000000000E1L;
long double C52r5 = 5.25000000000000000000E1L;
long double C105 = 1.05000000000000000000E2L;
#endif
#else /* LDOUBLE */
DOUBLE Cr375 =  3.75000000000000000000E-1;
DOUBLE C2r5 =  2.50000000000000000000E0;
DOUBLE C3r75 =  3.75000000000000000000E0;
DOUBLE C4r375 = 4.37500000000000000000E0;
DOUBLE C7r5 = 7.50000000000000000000E0;
DOUBLE C15 = 1.50000000000000000000E1;
DOUBLE C17r5 = 1.75000000000000000000E1;
DOUBLE C52r5 = 5.25000000000000000000E1;
DOUBLE C105 = 1.05000000000000000000E2;
#endif /* not LDOUBLE */


oblate( JD, yp, v, ymoon )
DOUBLE JD;
DOUBLE yp[], v[], ymoon[];
{
DOUBLE f, a, b, rem, res;
DOUBLE xyz[3], c, s;
int i;

/* x, y, z are the the (equatorial J2000)
 * solar system barycentric Cartesian coordinates
 * of the external point mass relative to the center of
 * the extended body.
 * Bx, By, Bz are the principal axes of inertia of the extended body.
 * The Euler angles are
 * theta, the angle from the z axis to the Bz axis
 * phi, measured in the x-y plane from the x axis to the line of nodes
 * psi, measured in the body's equatorial (Bx-By) plane
 *      from the line of nodes to the Bx axis.
 */
phi1 = yp[0];
phi = yp[1];
th1 = yp[2];
th = yp[3];
psi1 = yp[4];
v[1] = phi1; /* Copy the first derivative vector */
v[3] = th1;
v[5] = psi1;
/* Add slope representing the Moon's mean rotation */
psi1 += PSLP1;
/* Calculate psi = yp[5] + PSLP1*(JD - JDEPOCH).
 * Since only the sine and cosine of psi are needed, psi may
 * here be reduced accurately modulo the mean rotation period
 * of 2 pi / PSLP1 days.
 */
b = JD - JDEPOCH; /* elapsed days assumed arithmetically exact */
/* Subtract nearest integer number of rotation periods from the time */
f = FLOOR( b/PSLPI ); /* integer number of periods elapsed */
a = b - f * PSLPIA; /* arithmetic is exact while  f < 2^26  periods */
b = a - f * PSLPIB; /* note, PSLPIA + PSLPIB = 1 period */
psi = yp[5] + PSLP1*b; /* yp[5] is the integrator output for psi */
#if DEBUG
printf( "phi %.5e, th %.5e, psi %.5e\n", (double) phi, (double) th, (double) psi );
#endif


/* Construct matrix elements that will be used by mfigure()
 * to convert the extended-body-to-point-mass vector
 * from space (xyz) coordinates to body (Bxyz) coordinates.
 * This is a composite of three rotations:
 *   first by an angle phi about the z axis,
 *   second by an angle theta about the new x axis,
 *   third by an angle psi about the final z axis.
 * Each rotation is counter-clockwise, looking toward the
 * origin along the indicated axis.  See H. Goldstein,
 * _Classsical Mechanics_, Addison-Wesley, 1950, pp 107-109.
 */
sinpsi = SIN(psi);
cospsi = COS(psi);
sinth = SIN(th);
costh = COS(th);
sinphi = SIN(phi);
cosphi = COS(phi);
a = costh*sinphi;
b = costh*cosphi;
X00 = cospsi*cosphi - a*sinpsi;
X01 = cospsi*sinphi + b*sinpsi;
X02 = sinpsi*sinth;
X10 = -sinpsi*cosphi - a*cospsi;
X11 = -sinpsi*sinphi + b*cospsi;
X12 = cospsi*sinth;
X20 = sinth*sinphi;
X21 = -sinth*cosphi;
X22 = costh;

/* Initialize the torque accumulators.
 */
Nx = Zero;
Ny = Zero;
Nz = Zero;

/* Accelerations due to the figure of the Moon.
 */
/* effect of the point-mass Earth : 2e-12 au/d^2 */
xyz[0] = -ymoon[1];
xyz[1] = -ymoon[3];
xyz[2] = -ymoon[5];
mfigure( IEARTH, xyz, v );

/* effect of the point-mass Sun: 1e-17 au/d^2 */
xyz[0] = yp[6*ISUN+1] - yp[6*IMOON+1];
xyz[1] = yp[6*ISUN+3] - yp[6*IMOON+3];
xyz[2] = yp[6*ISUN+5] - yp[6*IMOON+5];
#if DEBUG
printf( "Moon->Sun %.5e %.5e %.5e\n", (double) xyz[0], (double) xyz[1], (double) xyz[2] );
#endif
/* a = Rij[ISUN][IMOON]; */
mfigure( ISUN, xyz, v );

/* Update accelerations of libration angles */
librate( v );


/* Oblateness of the earth
 */

/* effect of the point-mass Moon : 7e-11 au/d^2
 */
for( i=0; i<3; i++ )
	ofdate[i] = ymoon[2*i+1];
/* Precess from basic epoch to date */
precess( ofdate, JD, -1 );
nutate( JD, ofdate );
rem = Rij[IEARTH][IMOON];
legendre( ofdate, AE, rem, Je, xyz, 0 );
invnut( xyz );
precess( xyz, JD, 1 );
a = -GMs[IMOON];
i = 6 * IEARTH;
v[i] += a * xyz[0];
v[i+2] += a * xyz[1];
v[i+4] += a * xyz[2];

a = GMs[IEARTH];
#if DEBUG
printf( "E,M %.5e %.5e %.5e\n",
   (double )(a*xyz[0]), (double )(a*xyz[1]), (double )(a*xyz[2]) );
#endif
i = 6 * IMOON;
v[i] += a * xyz[0];
v[i+2] += a * xyz[1];
v[i+4] += a * xyz[2];

/* tides: 1e-15 au/d^2
 */
a = AE * rem;
b = a*a;
a = b*b*a*(One + One/EMRAT);
a = -Three*LOVENO*a*GMs[IMOON] * Rij3[IEARTH][IMOON];
xyz[0] = a * (ofdate[0] + PHASE * ofdate[1]);
xyz[1] = a * (ofdate[1] - PHASE * ofdate[0]);
xyz[2] = a * ofdate[2];
invnut( xyz );
precess( xyz, JD, 1 );
c = One/(One+EMRAT);
s = EMRAT*c;
i = 6 * IEARTH;
v[i] -= c * xyz[0];
v[i+2] -= c * xyz[1];
v[i+4] -= c * xyz[2];
#if DEBUG
printf( "T,M %.5e %.5e %.5e\n",
    (double )(s*xyz[0]), (double )(s*xyz[1]), (double )(s*xyz[2]) );
#endif
i = 6 * IMOON;
v[i] += s * xyz[0];
v[i+2] += s * xyz[1];
v[i+4] += s * xyz[2];


/* effect of the Sun: 7e-16 au/d^2
 */
res = Rij[ISUN][IEARTH];
ofdate[0] = yp[6*ISUN+1] - yp[6*IEARTH+1];
ofdate[1] = yp[6*ISUN+3] - yp[6*IEARTH+3];
ofdate[2] = yp[6*ISUN+5] - yp[6*IEARTH+5];
precess( ofdate, JD, -1 );
nutate( JD, ofdate );
legendre( ofdate, AE, res, Je, xyz, 0 );
invnut( xyz );
precess( xyz, JD, 1 );
a = -GMs[ISUN];
#if DEBUG
printf( "E,S %.5e %.5e %.5e\n",
 (double )(a*xyz[0]), (double )(a*xyz[1]), (double )(a*xyz[2]) );
#endif
i = 6 * IEARTH;
v[i] += a * xyz[0];
v[i+2] += a * xyz[1];
v[i+4] += a * xyz[2];
a = GMs[IEARTH];
i = 6 * ISUN;
v[i] += a * xyz[0];
v[i+2] += a * xyz[1];
v[i+4] += a * xyz[2];
}



mfigure( iobj, xyz, v )
int iobj; /* Index of distance from center of extended body to point mass */
DOUBLE xyz[]; /* space coordinates of point mass */
DOUBLE v[]; /* Velocity and acceleration state vector */
{
DOUBLE a, rem;
DOUBLE Bxyz[3], Fxyz[3];
int i;

/* Apply the transformation from space to body coordinates.
 */
Bxyz[0] = X00*xyz[0] + X01*xyz[1] + X02*xyz[2];
Bxyz[1] = X10*xyz[0] + X11*xyz[1] + X12*xyz[2];
Bxyz[2] = X20*xyz[0] + X21*xyz[1] + X22*xyz[2];

rem = Rij[iobj][IMOON];

/* Compute acceleration due to Moon's oblateness
 */
legendre( Bxyz, AM, rem, Jm, Fxyz, 1 );

/* The torque on the extended body, in body coordinates is r X F.
 * -GMiobj Fxyz is the acceleration of the extended body,
 * so the torque = -GMiobj r X Fxyz M
 * where M is the mass of the extended body.  M is not actually
 * included, since it will cancel out later.
 */
a = -GMs[iobj];
Nx += a*(Bxyz[1] * Fxyz[2] - Bxyz[2] * Fxyz[1]);
Ny += a*(Bxyz[2] * Fxyz[0] - Bxyz[0] * Fxyz[2]);
Nz += a*(Bxyz[0] * Fxyz[1] - Bxyz[1] * Fxyz[0]);

/* Transform from body coordinates to space coordinates
 * by the inverse (i.e., the transpose) of the previous transformation.
 */
xyz[0] = X00*Fxyz[0] + X10*Fxyz[1] + X20*Fxyz[2];
xyz[1] = X01*Fxyz[0] + X11*Fxyz[1] + X21*Fxyz[2];
xyz[2] = X02*Fxyz[0] + X12*Fxyz[1] + X22*Fxyz[2];

#if DEBUG
printf( "%d: %.5e %.5e %.5e\n",
 iobj, (double )(a*xyz[0]), (double )(a*xyz[1]), (double )(a*xyz[2]) );
#endif
i = 6 * IMOON;
v[i] += a * xyz[0];
v[i+2] += a * xyz[1];
v[i+4] += a * xyz[2];
a = GMs[IMOON];
i = 6 * iobj;
v[i] += a * xyz[0];
v[i+2] += a * xyz[1];
v[i+4] += a * xyz[2];
}




/* Expansion in Legendre functions.
 * The potential function is
 *                inf.  n    n
 *      GM  [      -    -   a      m           m           m        ]
 *  U = --- [ 1 -  >    >   --  ( C  cos mL + S  sin mL ) P (sin d) ]
 *       r  [      -    -    n     n           n           n        ]
 *                n=2  m=0  r
 *
 * where
 * a = mean radius
 * r = distance
 * L = longitude re principal axes
 * d = latitude  re principal axes
 * P = Legendre function
 *
 * The part GM/r has already been evaluated by the main program.
 * The part under the summations is converted to an acceleration
 * by taking its gradient in spherical coordinates and transforming
 * to Cartesian coordinates xi, eta, zeta.  The program then transforms
 * this into the principal axis body coordinate system.
 */
legendre( xyz, arad, r, Jn, Bxyz, assoc )
DOUBLE xyz[];  /* body coordinates of the point mass */
DOUBLE arad;   /* mean radius of extended body */
DOUBLE r;      /* center-to-center distance = |x,y,z| */
DOUBLE Jn[];   /* zonal harmonic coefficients */
DOUBLE Bxyz[]; /* output acceleration in body coordinates */
int assoc;     /* 1 = do tesseral harmonics */
{
DOUBLE P[3], dP[3], Pnm[8], dPnm[8], cosm[4], sinm[4];
DOUBLE coslat, sinlat, coslon, sinlon;
DOUBLE x, x2, x2p, sx2p, a, b, c, s, aor, aorn, np1;
DOUBLE xez0, xez1, xez2, t0, t1, t2;
int i, k, m;

/* Latitude of the point mass, in body coordinates:
 * Bz = r sin(lat)
 */
sinlat = xyz[2]*r;
x = sinlat;
x2 = sinlat * sinlat;
x2p = One - x2;
sx2p = SQRT( x2p );
coslat = sx2p;

/* Longitude of the point mass expressed in body coordinates:
 * By = r_xy sin(lon)
 * Bx = r_xy cos(lon)
 */
a = coslat/r; /* r is 1/distance */
coslon = xyz[0]/a;
sinlon = xyz[1]/a;

#if DEBUG
printf( "sinlat %.5e, sinlon %.5e\n", (double )sinlat, (double )sinlon );
#endif
/* Legendre polynomials.
 * P0(x) = 1
 * P1(x) = x
 * (n+1) Pn+1(x) = (2n+1) x Pn(x) - n Pn-1(x)
 */
P[0] = OneandaHalf*x2 - Half;          /* P2(x) */
P[1] = (C2r5 * x2 - OneandaHalf)*x;    /* P3(x) */
P[2] = (C4r375*x2 - C3r75)*x2 + Cr375; /* P4(x) */

/* Derivatives of Legendre polynomials
 * P'0(x) = 0
 * P'1(x) = 1
 * (x^2 - 1) P'n(x) = n[ x Pn(x) - Pn-1(x) ]
 */
dP[0] = Three * x;             /* dP2(x)/dx */
dP[1] = C7r5*x2 - OneandaHalf; /* dP3(x)/dx */
dP[2] = (C17r5*x2 - C7r5)*x;   /* dP4(x)/dx */

if( assoc )
	{
/*
 * Associated Legendre functions
 * for integer u:
 *
 *   u                  u/2   u        u
 *  P  (x)  =  (1 - x^2)     d Pv(x)/dx
 *   v
 *                                                    u
 * (Caution: some math textbooks multiply this by (-1)  ).
 * For v also an integer,
 *   u
 *  P  (x)  = 0   if u > v.
 *   v
 *
 */

	Pnm[0] = Three * x2p;               /* P^2_2(x) */
	Pnm[1] = sx2p * dP[1];              /* P^1_3(x) */
	Pnm[2] = C15 * x * x2p;             /* P^2_3(x) */
	b = x2p * sx2p;
	Pnm[3] = C15 * b;                   /* P^3_3(x) */
	Pnm[4] = sx2p * dP[2];              /* P^1_4(x) */
	Pnm[5] = x2p * (C52r5 * x2 - C7r5); /* P^2_4(x) */
	Pnm[6] = C105 * x * b;              /* P^3_4(x) */
	Pnm[7] = C105 * x2p * x2p;          /* P^4_4(x) */

/* Derivatives of associated Legendre functions
 */
	a = -x/x2p;
	b = One/sx2p;
	dPnm[0] = -Six*x;                    /* dP^2_2(x)/dx */
	dPnm[1] =     a*Pnm[1] + b*Pnm[2];   /* dP^1_3(x)/dx */
	dPnm[2] = Two*a*Pnm[2] + b*Pnm[3];   /* dP^2_3(x)/dx */
	dPnm[3] = Three*a*Pnm[3];            /* dP^3_3(x)/dx */
	dPnm[4] =     a*Pnm[4] + b*Pnm[5];   /* dP^1_4(x)/dx */
	dPnm[5] = Two*a*Pnm[5] + b*Pnm[6];   /* dP^2_4(x)/dx */
	dPnm[6] = Three*a*Pnm[6] + b*Pnm[7]; /* dP^3_4(x)/dx */
	dPnm[7] = Four*a*Pnm[7];             /* dP^4_4(x)/dx */


/* sines and cosines of m * longitude, m = 1, 2, 3, 4
 */
	c = coslon;
	s = sinlon;
	cosm[0] = c;
	sinm[0] = s;
	for( i=1; i<4; i++ )
		{
		b = sinlon * c + coslon * s;
		c = coslon * c - sinlon * s;
		s = b;
		cosm[i] = c;
		sinm[i] = s;
		}
	}

xez0 = Zero; /* xi */
xez1 = Zero; /* eta */
xez2 = Zero; /* zeta */
np1 = Two;
aor = arad * r;
aorn = aor;
k = 0;
for( i=0; i<3; i++ )
	{
	aorn *= aor;
	np1 += One;
/* zonal harmonics */
	a = Jn[i];
	t0 = a * np1 * P[i];
	t1 = Zero;
	t2 = -coslat * a * dP[i];

	if( assoc )
		{
/* tesseral harmonics */
		for( m=0; m<=(i+1); m++ )
			{
			if( (i == 0) && (m == 0) )
				continue;
			a = Pnm[k];
			b = Cnm[k]*cosm[m] + Snm[k]*sinm[m];
			t0 -= np1*a*b;
			t1 += ((m+1)*a*(-Cnm[k]*sinm[m]+Snm[k]*cosm[m]))/coslat;
			t2 += coslat*dPnm[k]*b;
			k += 1;
			}
		}
	xez0 += aorn * t0;
	xez1 += aorn * t1;
	xez2 += aorn * t2;
	}
r = r*r; /* 1/r^2 */
xez0 *= r;
if( assoc )
	xez1 *= r;
xez2 *= r;
/* Rotate xi,eta,zeta into the Bx,By,Bz body coordinate system.
 */
/* First rotate counterclockwise about the eta axis by the latitude. */
Bxyz[0] = coslat*xez0 - sinlat*xez2;
Bxyz[1] = xez1;
Bxyz[2] = sinlat*xez0 + coslat*xez2;
/* Then rotate clockwise about the new z axis by the longitude. */
b  = coslon*Bxyz[0] - sinlon*Bxyz[1];
Bxyz[1] = sinlon*Bxyz[0] + coslon*Bxyz[1];
Bxyz[0] = b;
}


/* Differential equation for Lunar librations
 */
extern DOUBLE AMR2, BMR2, CMR2;

librate( v )
DOUBLE v[];
{
DOUBLE omBx, omBy, omBz, omB1x, omB1y, omB1z;
DOUBLE phi2, th2, psi2, a;

/* The angular velocity of the extended body,
 * expressed in body coordinates.  The suffix 1 indicates
 * differentiation with respect to time.
 */
a = phi1 * sinth;
omBx = a*sinpsi + th1*cospsi;
omBy = a*cospsi - th1*sinpsi;
omBz = phi1*costh + psi1;

/*
 * The ratio of the torque N to the principal moment of
 * inertia C is calculated by
 *
 *    N          r X F          r X (M a)         r X a
 *   ----  =  -----------  =  ------------  =  ------------
 *    C       C/MR^2 MR^2     C/MR^2 MR^2       C/MR^2 R^2
 *
 * where a is the acceleration due to oblateness,
 * M is the mass of the extended body,
 * R is the mean radius of the extended body,
 * r is the vector from the center of the body to the point mass.
 * r X a = (Nx,Ny,Nz), calculated earlier.
 */

/* Euler's equations for the change in angular velocity.
 */
omB1x = LGAMBET*omBy*omBz + Nx/AMR2;
omB1y = LBET*omBz*omBx + Ny/BMR2;
omB1z = -LGAM*omBx*omBy + Nz/CMR2;

/* Differential equations for the Euler angles
 */
phi2 = (omB1x*sinpsi + omB1y*cospsi + th1*(psi1 - phi1*costh))/sinth;
th2 = omB1x*cospsi - omB1y*sinpsi - a*psi1;
psi2 = omB1z - phi2*costh + a*th1;

v[0] = phi2;
v[2] = th2;
v[4] = psi2;
#if DEBUG
printf( "phi2 %.12e, th2 %.12e, psi2 %.12e\n",
   (double) phi2, (double) th2, (double) psi2 );
#endif
}
