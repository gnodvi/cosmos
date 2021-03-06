#ifndef PROTOS_H
#define PROTOS_H
/*
 *   This file was automatically generated by version 1.7 of cextract.
 *   Manual editing not recommended.
 *
 *   Created: Wed May 12 21:05:54 2004
 */

#include "prec.h"

/* Adams-Bashforth-Moulton integration formulas. */
extern int adstart ( DOUBLE h, DOUBLE yn[], DOUBLE work[], int neq, int ord, DOUBLE t );
extern DOUBLE adstep ( DOUBLE *t, DOUBLE yn[], int nequat );
extern DOUBLE difpol ( DOUBLE *diffn, int k, DOUBLE t );

extern int divdif ( DOUBLE vec[], int k, DOUBLE *diffn );

extern int dupdate ( DOUBLE *diffn, int k, DOUBLE f );
extern DOUBLE intpol ( DOUBLE *diffn, DOUBLE *coeffs, int k );
extern int rinternal ( int f, DOUBLE *work );
extern int setptrs ( DOUBLE *w );
extern int sinternal ( int f );
extern int vcopy ( DOUBLE *p, DOUBLE *q, int n );

/*							asinl.c */
extern long double acosl ( long double x );
extern long double asinl ( long double x );

/*							atanl.c
 *	Inverse circular tangent, long double precision */
extern long double atan2l ( long double x, long double y );
extern long double atanl ( long double x );

/* Obliquity of the ecliptic at Julian date J */
extern int epsiln ( DOUBLE J );

/* Include file for extended precision arithmetic programs.  */
extern void dectoe ( unsigned short *d, unsigned short *e );
extern void etodec ( unsigned short *x, unsigned short *d );
extern void todec ( unsigned short *x, unsigned short *y );

/*							mconf.h
 *	Common include file for math routines */
extern int findcenter ( DOUBLE t, DOUBLE yin[] );

/*							ceill()
 *							floorl()
 *							frexpl()
 *							ldexpl()
 *							fabsl()
 *							signbitl()
 *							isnanl()
 *							isfinitel()
 */
extern long double ceill ( long double x );
extern long double fabsl ( long double x );
extern long double floorl ( long double x );
extern long double frexpl ( long double x, int *pw2 );
extern long double ldexpl ( long double x, int pw2 );

/*							ieee.c
 *
 *    Extended precision IEEE binary floating point arithmetic routines */
extern void asctoe ( char *s, unsigned short *y );
extern void asctoe24 ( char *s, unsigned short *y );
extern void asctoe53 ( char *s, unsigned short *y );
extern void asctoe64 ( char *s, unsigned short *y );
extern void asctoeg ( char *ss, unsigned short *y, int oprec );
extern void e24toasc ( unsigned short x[], char *string, int ndigs );
extern void e24toe ( unsigned short *e, unsigned short *y );
extern void e53toasc ( unsigned short x[], char *string, int ndigs );
extern void e53toe ( unsigned short *e, unsigned short *y );
extern void e64toasc ( unsigned short x[], char *string, int ndigs );
extern void e64toe ( unsigned short *e, unsigned short *y );
extern void eabs ( unsigned short x[] );
extern void eadd ( unsigned short *a, unsigned short *b, unsigned short *c );
extern void eadd1 ( unsigned short *a, unsigned short *b, unsigned short *c );
extern void eaddm ( unsigned short *x, unsigned short *y );
extern void eclear ( unsigned short *x );
extern void ecleaz ( unsigned short *xi );
extern void ecleazs ( unsigned short *xi );
extern int ecmp ( unsigned short *a, unsigned short *b );
extern int ecmpm ( unsigned short *a, unsigned short *b );
extern void ediv ( unsigned short *a, unsigned short *b, unsigned short *c );
extern int edivm ( unsigned short den[], unsigned short num[] );
extern void efloor ( unsigned short x[], unsigned short y[] );
extern void efrexp ( unsigned short x[], long *exp, unsigned short s[] );
extern void eifrac ( unsigned short *x, long *i, unsigned short *frac );
extern void einfin ( unsigned short *x );
extern void einit ( void );
extern void eiremain ( unsigned short den[], unsigned short num[] );
extern int eisinf ( unsigned short x[] );
extern int eisneg ( unsigned short x[] );
extern void eldexp ( unsigned short x[], long pwr2, unsigned short y[] );
extern void emdnorm ( unsigned short s[], int lost, int subflg, long exp, int rcntrl );
extern void emov ( unsigned short *a, unsigned short *b );
extern void emovi ( unsigned short *a, unsigned short *b );
extern void emovo ( unsigned short *a, unsigned short *b );
extern void emovz ( unsigned short *a, unsigned short *b );
extern void emul ( unsigned short *a, unsigned short *b, unsigned short *c );
extern int emulm ( unsigned short a[], unsigned short b[] );
extern void eneg ( unsigned short x[] );
extern int enormlz ( unsigned short x[] );
extern void eremain ( unsigned short a[], unsigned short b[], unsigned short c[] );
extern void eround ( unsigned short *x, unsigned short *y );
extern void eshdn1 ( unsigned short *x );
extern void eshdn6 ( unsigned short *x );
extern void eshdn8 ( unsigned short *x );
extern int eshift ( unsigned short *x, int sc );
extern void eshup1 ( unsigned short *x );
extern void eshup6 ( unsigned short *x );
extern void eshup8 ( unsigned short *x );
extern void esub ( unsigned short *a, unsigned short *b, unsigned short *c );
extern void esubm ( unsigned short *x, unsigned short *y );
extern void etoasc ( unsigned short x[], char *string, int ndigs );
extern void etoe24 ( unsigned short *x, unsigned short *e );
extern void etoe53 ( unsigned short *x, unsigned short *e );
extern void etoe64 ( unsigned short *x, unsigned short *e );
extern void euifrac ( unsigned short *x, long *i, unsigned short *frac );
extern void ltoe ( long *lp, unsigned short *y );
extern void ultoe ( unsigned long *lp, unsigned short *y );

/* Orbits of the five minor planets as implemented in the
 * JPL DE102 and DE118/DE200 numerical integrations. */
extern int aroids ( DOUBLE JD, DOUBLE y[] );
extern int r102_200 ( DOUBLE vec[] );
extern int r118_200 ( DOUBLE vec[] );
extern int r200_102 ( DOUBLE vec[] );
extern int r200_118 ( DOUBLE vec[] );
extern int mtherr ( char *str, int code );

/* One-term nutation program uses only the 18.6 year term. */
extern int invnut ( DOUBLE p1[] );
extern int nutate ( DOUBLE J, DOUBLE p[] );
extern int nutlo ( DOUBLE J );

/* oblate.c */
extern int legendre ( DOUBLE xyz[], DOUBLE arad, DOUBLE r, DOUBLE Jn[], DOUBLE Bxyz[], int assoc );
extern int librate ( DOUBLE v[] );
extern int mfigure ( int iobj, DOUBLE xyz[], DOUBLE v[] );
extern int oblate ( DOUBLE JD, DOUBLE yp[], DOUBLE v[], DOUBLE ymoon[] );

/*							mconf.h
 *
 *	Common include file for math routines */
extern long double p1evll ( long double x, long double *P, int n );
extern long double polevll ( long double x, long double *P, int n );

/* Precession of the equinox and ecliptic
 * from epoch Julian date J to or from J2000.0 */
extern int precess ( DOUBLE R[], DOUBLE J, int direction );
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extern int rdnums ( void );
void set_const_by_hand (void );
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extern int reltiv ( DOUBLE y[], DOUBLE v[] );

/* Runge-Kutta numerical integration
 * of a system of ordinary differential equations. */
extern int rkstart ( int neq, DOUBLE work[] );
extern int rungek ( int neqn, DOUBLE x, DOUBLE yold[], DOUBLE h, DOUBLE ynew[], DOUBLE delta[] );

/*							sinl.c
 *
 *	Circular sine, long double precision */
extern long double cosl ( long double x );
extern long double sinl ( long double x );

/*							sqrtl.c
 *
 *	Square root, long double precision */
extern long double sqrtl ( long double x );

/* Main program for numerical integration of the solar system. */
extern int distances ( DOUBLE y[] );
extern int fixsun ( DOUBLE y[], DOUBLE v[] );
extern int fromemb ( DOUBLE yinp[], DOUBLE y[] );
extern int func ( DOUBLE t, DOUBLE yin[], DOUBLE v[] );
extern void resstate ( void );
extern void savstate ( void );
extern int wfile (int fd,  DOUBLE t, DOUBLE y[] );

/*							tanl.c
 *
 *	Circular tangent, long double precision */
extern long double cotl ( long double x );
extern long double tancotl ( long double xx, int cotflg );
extern long double tanl ( long double x );

/*							atan2()
 *
 *	Quadrant correct inverse circular tangent */
extern DOUBLE zatan2 ( DOUBLE x, DOUBLE y );
#endif
