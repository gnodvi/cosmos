//------------------------------------------------------------------------------
/* Adams-Bashforth-Moulton integration formulas.
 *
 * Reference:
 * Shampine, L. F. and M. K. Gordon, _Computer Solution of
 * Ordinary Differential Equations_, W. H. Freeman, 1975.
 *
 * Program by Steve Moshier.
 */
//------------------------------------------------------------------------------

#ifndef NOINCS
#include "mconf.h"
#if MSC
#include <io.h>
#endif
#if UNIX
#include <unistd.h>
#endif
#include <stdlib.h>
#include "prec.h"
#include "int.h"
#include "const.h"
#endif

/* Divided differences */
#define N 19


/*Predictor coefficients:
 1.0
 1.0 / 2.0
 5.0 / 12.0
 3.0 / 8.0
 251.0 / 720.0
 95.0 / 288.0
 19087.0 / 60480.0
 5257.0 / 17280.0
 1070017.0 / 3628800.0
 25713.0 / 89600.0
 26842253.0 / 95800320.0
 4777223.0 / 17418240.0
 703604254357.0 / 2615348736000.0
 106364763817.0 / 402361344000.0
 1166309819657.0 / 4483454976000.0
 2.5221445e7 / 9.8402304e7
 8.092989203533249e15 / 3.201186852864e16
 8.5455477715379e13 / 3.4237292544e14
 1.2600467236042756559e19 / 5.109094217170944e19
*/

/* Corrector coefficients:
   1.0
  -1.0 / 2.0
  -1.0 / 12.0
  -1.0 / 24.0
  -19.0 / 720.0
  -3.0 / 160.0
  -863.0 / 60480.0
  -275.0 / 24192.0
  -33953.0 / 3628800.0
  -8183.0 / 1036800.0
  -3250433.0 / 479001600.0
  -4671.0 / 788480.0
  -13695779093.0 / 2615348736000.0
  -2224234463.0 / 475517952000.0
  -132282840127.0 / 31384184832000.0
  -2639651053.0 / 689762304000.0
  1.11956703448001e14 / 3.201186852864e16
  5.0188465e7 / 1.5613165568e10
  2.334028946344463e15 / 7.86014494949376e17
*/
#if LDOUBLE
#if IBMPC
static short precof[] = {
0x0000,0x0000,0x0000,0x8000,0x3fff, XPD
0x0000,0x0000,0x0000,0x8000,0x3ffe, XPD
0x5555,0x5555,0x5555,0xd555,0x3ffd, XPD
0x0000,0x0000,0x0000,0xc000,0x3ffd, XPD
0xd27d,0x7d27,0x27d2,0xb27d,0x3ffd, XPD
0x38e4,0xe38e,0x8e38,0xa8e3,0x3ffd, XPD
0x5440,0xea99,0x43fe,0xa195,0x3ffd, XPD
0x3519,0x6dfc,0x518a,0x9bc3,0x3ffd, XPD
0xb126,0x0fb1,0xf045,0x96f8,0x3ffd, XPD
0x80bb,0x54d8,0x721a,0x92ee,0x3ffd, XPD
0x6246,0xcff2,0x02c2,0x8f75,0x3ffd, XPD
0xac15,0xb603,0x8869,0x8c6c,0x3ffd, XPD
0xdc52,0x2596,0x2625,0x89be,0x3ffd, XPD
0x2820,0xc6b2,0x0f57,0x8759,0x3ffd, XPD
0xf898,0xbc14,0x9903,0x8530,0x3ffd, XPD
0x4d10,0xe1e8,0xff92,0x833a,0x3ffd, XPD
0x883c,0x772e,0x97fc,0x8170,0x3ffd, XPD
0xe995,0xbf9f,0x86c4,0xff96,0x3ffc, XPD
0x78b9,0x72bf,0x1a81,0xfc8c,0x3ffc, XPD
};
static short corcof[] = {
0x0000,0x0000,0x0000,0x8000,0x3fff, XPD
0x0000,0x0000,0x0000,0x8000,0xbffe, XPD
0xaaab,0xaaaa,0xaaaa,0xaaaa,0xbffb, XPD
0xaaab,0xaaaa,0xaaaa,0xaaaa,0xbffa, XPD
0xd82e,0x2d82,0x82d8,0xd82d,0xbff9, XPD
0x999a,0x9999,0x9999,0x9999,0xbff9, XPD
0x9474,0x1e9c,0x473f,0xe9c9,0xbff8, XPD
0xe4e9,0x93a3,0x4e8f,0xba3e,0xbff8, XPD
0x7e46,0xc950,0x28ab,0x994c,0xbff8, XPD
0x0d67,0x5b26,0xc557,0x814f,0xbff8, XPD
0x9d4e,0x3987,0xd5e1,0xde5b,0xbff7, XPD
0x8c4d,0x7bad,0x9646,0xc21e,0xbff7, XPD
0xf0ac,0x1b33,0x9124,0xab98,0xbff7, XPD
0x0c7c,0xb92d,0xb357,0x9945,0xbff7, XPD
0xe201,0xa74b,0x9502,0x8a1d,0xbff7, XPD
0xc3dc,0x1655,0xb86d,0xfacc,0xbff6, XPD
0x6a55,0x5ce2,0xcb35,0xe533,0x3ff6, XPD
0xb891,0xaf49,0x4d0b,0xd2aa,0x3ff6, XPD
0x370d,0x381c,0x10d3,0xc29b,0x3ff6, XPD
};
#endif
#if MIEEE
static long precof[3*19] = {
0x3fff0000,0x80000000,0x00000000,
0x3ffe0000,0x80000000,0x00000000,
0x3ffd0000,0xd5555555,0x55555555,
0x3ffd0000,0xc0000000,0x00000000,
0x3ffd0000,0xb27d27d2,0x7d27d27d,
0x3ffd0000,0xa8e38e38,0xe38e38e4,
0x3ffd0000,0xa19543fe,0xea995440,
0x3ffd0000,0x9bc3518a,0x6dfc3519,
0x3ffd0000,0x96f8f045,0x0fb1b126,
0x3ffd0000,0x92ee721a,0x54d880bb,
0x3ffd0000,0x8f7502c2,0xcff26246,
0x3ffd0000,0x8c6c8869,0xb603ac15,
0x3ffd0000,0x89be2625,0x2596dc52,
0x3ffd0000,0x87590f57,0xc6b22820,
0x3ffd0000,0x85309903,0xbc14f898,
0x3ffd0000,0x833aff92,0xe1e84d10,
0x3ffd0000,0x817097fc,0x772e883c,
0x3ffc0000,0xff9686c4,0xbf9fe995,
0x3ffc0000,0xfc8c1a81,0x72bf78b9,
};
static long corcof[3*19] = {
0x3fff0000,0x80000000,0x00000000,
0xbffe0000,0x80000000,0x00000000,
0xbffb0000,0xaaaaaaaa,0xaaaaaaab,
0xbffa0000,0xaaaaaaaa,0xaaaaaaab,
0xbff90000,0xd82d82d8,0x2d82d82e,
0xbff90000,0x99999999,0x9999999a,
0xbff80000,0xe9c9473f,0x1e9c9474,
0xbff80000,0xba3e4e8f,0x93a3e4e9,
0xbff80000,0x994c28ab,0xc9507e46,
0xbff80000,0x814fc557,0x5b260d67,
0xbff70000,0xde5bd5e1,0x39879d4e,
0xbff70000,0xc21e9646,0x7bad8c4d,
0xbff70000,0xab989124,0x1b33f0ac,
0xbff70000,0x9945b357,0xb92d0c7c,
0xbff70000,0x8a1d9502,0xa74be201,
0xbff60000,0xfaccb86d,0x1655c3dc,
0x3ff60000,0xe533cb35,0x5ce26a55,
0x3ff60000,0xd2aa4d0b,0xaf49b891,
0x3ff60000,0xc29b10d3,0x381c370d,
};
#endif /* ldouble MIEEE */
#if UNK
static long double precof[19] = {
  1.0000000000000000000000000E0L,
  5.0000000000000000000000000E-1L,
  4.1666666666666666666666667E-1L,
  3.7500000000000000000000000E-1L,
  3.4861111111111111111111111E-1L,
  3.2986111111111111111111111E-1L,
  3.1559193121693121693121693E-1L,
  3.0422453703703703703703704E-1L,
  2.9486800044091710758377425E-1L,
  2.8697544642857142857142857E-1L,
  2.8018959644393672171449949E-1L,
  2.7426554003159905937683715E-1L,
  2.6902884677364877431014997E-1L,
  2.6435134836660650979434048E-1L,
  2.6013639612760103693745669E-1L,
  2.5630949657438915251415251E-1L,
  2.5281214672903923486004312E-1L,
  2.4959765029771566740858701E-1L,
  2.4662820258225745779739230E-1L,
};
static long double corcof[19] = {
   1.0000000000000000000000000E0L,
  -5.0000000000000000000000000E-1L,
  -8.3333333333333333333333333E-2L,
  -4.1666666666666666666666667E-2L,
  -2.6388888888888888888888889E-2L,
  -1.8750000000000000000000000E-2L,
  -1.4269179894179894179894180E-2L,
  -1.1367394179894179894179894E-2L,
  -9.3565365961199294532627866E-3L,
  -7.8925540123456790123456790E-3L,
  -6.7858499846347068569290792E-3L,
  -5.9240564123376623376623377E-3L,
  -5.2366932579502850666871831E-3L,
  -4.6774984070422645158094894E-3L,
  -4.2149522390054728568837916E-3L,
  -3.8268995532118844233041764E-3L,
  3.4973498453499176541093925E-3L,
  3.2144964313235674514561069E-3L,
  2.9694477154582096111947101E-3L,
};
#endif
#else /* LDOUBLE */
#if UNK
static DOUBLE precof[19] = {
 1.00000000000000000000E0,
 5.00000000000000000000E-1,
 4.16666666666666666667E-1,
 3.75000000000000000000E-1,
 3.48611111111111111111E-1,
 3.29861111111111111111E-1,
 3.15591931216931216931E-1,
 3.04224537037037037037E-1,
 2.94868000440917107584E-1,
 2.86975446428571428571E-1,
 2.80189596443936721714E-1,
 2.74265540031599059377E-1,
 2.69028846773648774310E-1,
 2.64351348366606509794E-1,
 2.60136396127601036937E-1,
 2.56309496574389152514E-1,
 2.52812146729039234860E-1,
 2.49597650297715667409E-1,
 2.46628202582257457797E-1,
};
static DOUBLE corcof[19] = {
 1.00000000000000000000E0,
-5.00000000000000000000E-1,
-8.33333333333333333333E-2,
-4.16666666666666666667E-2,
-2.63888888888888888889E-2,
-1.87500000000000000000E-2,
-1.42691798941798941799E-2,
-1.13673941798941798942E-2,
-9.35653659611992945326E-3,
-7.89255401234567901235E-3,
-6.78584998463470685693E-3,
-5.92405641233766233766E-3,
-5.23669325795028506669E-3,
-4.67749840704226451581E-3,
-4.21495223900547285688E-3,
-3.82689955321188442330E-3,
 3.49734984534991765411E-3,
 3.21449643132356745146E-3,
 2.96944771545820961119E-3,
};
#endif /* double UNK */
#if DEC
static short precof[19*4] = {
0040200,0000000,0000000,0000000,
0040000,0000000,0000000,0000000,
0037725,0052525,0052525,0052525,
0037700,0000000,0000000,0000000,
0037662,0076447,0151175,0023722,
0037650,0161616,0034343,0107071,
0037641,0112503,0177352,0114524,
0037633,0141521,0105155,0176065,
0037626,0174360,0042417,0130661,
0037622,0167162,0015124,0154201,
0037617,0072402,0141317,0171142,
0037614,0066210,0064666,0001654,
0037611,0137046,0022445,0113334,
0037607,0054417,0053706,0131050,
0037605,0030231,0001674,0012371,
0037603,0035377,0111341,0164115,
0037601,0070227,0176167,0027210,
0037577,0113206,0142277,0117752,
0037574,0106032,0100562,0137571,
};
static short corcof[19*4] = {
0040200,0000000,0000000,0000000,
0140000,0000000,0000000,0000000,
0137252,0125252,0125252,0125253,
0137052,0125252,0125252,0125253,
0136730,0026602,0154055,0101330,
0136631,0114631,0114631,0114632,
0136551,0144507,0037436,0116224,
0136472,0037116,0107623,0121745,
0136431,0046050,0125711,0050176,
0136401,0047705,0053533,0023015,
0136336,0055725,0160471,0103635,
0136302,0017226,0043173,0126614,
0136253,0114221,0022033,0031761,
0136231,0042663,0053671,0026414,
0136212,0016625,0001247,0045742,
0136172,0146270,0066426,0052704,
0036145,0031713,0032534,0161152,
0036122,0125115,0005657,0044671,
0036102,0115420,0151470,0016067,
};
#endif /* double DEC */
#if IBMPC
static short precof[19*4] = {
0x0000,0x0000,0x0000,0x3ff0,
0x0000,0x0000,0x0000,0x3fe0,
0xaaab,0xaaaa,0xaaaa,0x3fda,
0x0000,0x0000,0x0000,0x3fd8,
0xa4fa,0xfa4f,0x4fa4,0x3fd6,
0x71c7,0xc71c,0x1c71,0x3fd5,
0x532b,0x7fdd,0x32a8,0x3fd4,
0xbf87,0x314d,0x786a,0x3fd3,
0xf636,0x08a1,0xdf1e,0x3fd2,
0x9b10,0x434a,0x5dce,0x3fd2,
0xfe4c,0x5859,0xeea0,0x3fd1,
0xc076,0x0d36,0x8d91,0x3fd1,
0xb2dc,0xc4a4,0x37c4,0x3fd1,
0xd645,0xeaf8,0xeb21,0x3fd0,
0x829f,0x2077,0xa613,0x3fd0,
0x3d0a,0xf25c,0x675f,0x3fd0,
0xe5d1,0xff8e,0x2e12,0x3fd0,
0xf3fd,0xd897,0xf2d0,0x3fcf,
0x57ef,0x502e,0x9183,0x3fcf,
};
static short corcof[19*4] = {
0x0000,0x0000,0x0000,0x3ff0,
0x0000,0x0000,0x0000,0xbfe0,
0x5555,0x5555,0x5555,0xbfb5,
0x5555,0x5555,0x5555,0xbfa5,
0xb05b,0x5b05,0x05b0,0xbf9b,
0x3333,0x3333,0x3333,0xbf93,
0xd393,0xe7e3,0x3928,0xbf8d,
0x747d,0xd1f2,0x47c9,0xbf87,
0x2a10,0x1579,0x2985,0xbf83,
0x64c2,0xaaeb,0x29f8,0xbf80,
0x30f4,0xbc27,0xcb7a,0xbf7b,
0x75b2,0xc8cf,0x43d2,0xbf78,
0x667e,0x2483,0x7312,0xbf75,
0x25a2,0x6af7,0x28b6,0xbf73,
0xe97c,0xa054,0x43b2,0xbf71,
0xcab8,0x0da2,0x5997,0xbf6f,
0x9c4d,0x66ab,0xa679,0x3f6c,
0xe937,0xa175,0x5549,0x3f6a,
0x0387,0x1a67,0x5362,0x3f68,
};
#endif /* double IBMPC */
#if MIEEE
static short precof[19*4] = {
0x3ff0,0x0000,0x0000,0x0000,
0x3fe0,0x0000,0x0000,0x0000,
0x3fda,0xaaaa,0xaaaa,0xaaab,
0x3fd8,0x0000,0x0000,0x0000,
0x3fd6,0x4fa4,0xfa4f,0xa4fa,
0x3fd5,0x1c71,0xc71c,0x71c7,
0x3fd4,0x32a8,0x7fdd,0x532b,
0x3fd3,0x786a,0x314d,0xbf87,
0x3fd2,0xdf1e,0x08a1,0xf636,
0x3fd2,0x5dce,0x434a,0x9b10,
0x3fd1,0xeea0,0x5859,0xfe4c,
0x3fd1,0x8d91,0x0d36,0xc076,
0x3fd1,0x37c4,0xc4a4,0xb2dc,
0x3fd0,0xeb21,0xeaf8,0xd645,
0x3fd0,0xa613,0x2077,0x829f,
0x3fd0,0x675f,0xf25c,0x3d0a,
0x3fd0,0x2e12,0xff8e,0xe5d1,
0x3fcf,0xf2d0,0xd897,0xf3fd,
0x3fcf,0x9183,0x502e,0x57ef,
};
static short corcof[19*4] = {
0x3ff0,0x0000,0x0000,0x0000,
0xbfe0,0x0000,0x0000,0x0000,
0xbfb5,0x5555,0x5555,0x5555,
0xbfa5,0x5555,0x5555,0x5555,
0xbf9b,0x05b0,0x5b05,0xb05b,
0xbf93,0x3333,0x3333,0x3333,
0xbf8d,0x3928,0xe7e3,0xd393,
0xbf87,0x47c9,0xd1f2,0x747d,
0xbf83,0x2985,0x1579,0x2a10,
0xbf80,0x29f8,0xaaeb,0x64c2,
0xbf7b,0xcb7a,0xbc27,0x30f4,
0xbf78,0x43d2,0xc8cf,0x75b2,
0xbf75,0x7312,0x2483,0x667e,
0xbf73,0x28b6,0x6af7,0x25a2,
0xbf71,0x43b2,0xa054,0xe97c,
0xbf6f,0x5997,0x0da2,0xcab8,
0x3f6c,0xa679,0x66ab,0x9c4d,
0x3f6a,0x5549,0xa175,0xe937,
0x3f68,0x5362,0x1a67,0x0387,
};
#endif /* double MIEEE */
#endif /* not LDOUBLE */

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
/* Compute zeroth through kth backward differences
 * of the data in the input array
 */
//------------------------------------------------------------------------------
divdif (
        DOUBLE vec[], // input array of k+1 data items 
        int k,
        DOUBLE *diffn // output array of ith differences 
        )
{
  DOUBLE diftbl[N];
  DOUBLE *p, *q;
  DOUBLE y;
  int i, o;

  /* Copy the given data (zeroth difference) into temp array
    */
  p = diftbl;
  q = vec;
  for( i=0; i<=k; i++ )
    *p++ = *q++;

  /* On the first outer loop, k-1 first differences are calculated.
    * These overwrite the original data in the temp array.
    */
  o = k;
  for( o=k; o>0; o-- )
  {
    p = diftbl;
    q = p;
    for( i=0; i<o; i++ )
    {
      y = *p++;
      *q++ = *p - y;
    }
    *diffn++ = *p; /* copy out the last (undifferenced) item */
#if DEBUG
    printf( "%.5e ", *p );
#endif
  }

#if DEBUG
  printf( "%.5e\n", *(q-1) );
#endif
  *diffn++ = *(q-1);
}

//------------------------------------------------------------------------------
/* Update array of differences, given new data value.
 * diffn is an array of k+1 differences, starting with the
 * zeroth difference (the previous original data value).
 */
//------------------------------------------------------------------------------
dupdate( diffn, k, f )
  register DOUBLE *diffn;  /* input and output array of differences */
int k; /* max order of differences */
DOUBLE f; /* new data point (zeroth difference) */
{

  DOUBLE new, old;
  int i;

  new = f;
  for( i=0; i<k; i++ )
  {
    old = *diffn;
    *diffn++ = new;
    #if DEBUG
    printf( "%.5e ", new );
    #endif
    new = new - old;
  }

#if DEBUG
  printf( "%.5e\n", new );
#endif

  *diffn++ = new;
}

//------------------------------------------------------------------------------

/* Evaluate the interpolating polynomial
 *
 *              (x - x )
 *                    n    1
 * P(x) = f  +  --------  D f   +  ...
 *         n       h         n
 *
 *     (x - x )(x - x   )...(x - x     )
 *           n       n-1          n+2-k    k-1
 *  +  ---------------------------------  D    f
 *                   k-1                        n
 *                  h     (k-1)!
 *
 *
 *         j
 *  where D denotes the jth backward difference, see dupdate(), and
 *
 *  f   =   f( x , y(x ) )  is the interpolated derivative y'(x ) .
 *   n          n     n                                        n
 *
 * The subroutine argument t is linearly scaled so that t = 1.0
 * will evaluate the polynomial at x = x_n + h,
 * t = 0.0 corresponds to x = x_n, etc.
 */
//------------------------------------------------------------------------------

DOUBLE difpol ( diffn, k, t )
  DOUBLE *diffn;
    int k; /* differences go up to order k-1 */
    DOUBLE t; /* scaled argument */
{

  DOUBLE f, fac, s, u;
  int i;

  f = *diffn++; /* the zeroth difference = nth data point */
  u = One;
  /*s = x/h - n;*/
  s = One;	/* to evaluate the polynomial at x = xn + h */
  fac = One;

  for( i=1; i<k; i++ )
  {
    if( s == 0 )
      break;
    u *= s / fac;
    f += u  *  (*diffn++);
    fac += One;
    s += One;
  }

  return( f );

}
//------------------------------------------------------------------------------
/* Integrate the interpolating polynomial from x[n] to x[n+1]
 * to obtain the change in the integrated function from y[n] to y[n+1]
 * given by:
 *
 *                     k-1
 *                      -       i
 *  y    =  y   +   h   >   c  D  f    .
 *   n+1     n          -    i     n
 *                     i=0
 *
 *  This subroutine returns the summation term, not multiplied by h.
 *  The coefficients c_i of the integration formula are either
 *  precof[] or corcof[], given above.
 */
//------------------------------------------------------------------------------

DOUBLE intpol ( diffn, coeffs, k )

  DOUBLE *diffn; /* array of backward differences */
  DOUBLE *coeffs; /* coefficients of integration formula */
  int k; /* differences used go up to order k-1 */
{

  DOUBLE s;
  int i;

  s = Zero;
  coeffs += k;
  diffn += k;

  for( i=0; i<k; i++ )
  {
    s += (*--coeffs) * (*--diffn);
  }

  return( s );
}
//------------------------------------------------------------------------------
/* Copy array of n elements from p to q.
 */
//------------------------------------------------------------------------------
vcopy ( p, q, n )

  register DOUBLE *p, *q;
  register int n;
{
  do
    *q++ = *p++;
  while( --n );

}
//------------------------------------------------------------------------------
/* Adams initialization program.
 */
//------------------------------------------------------------------------------

/* Addresses within the work array */
DOUBLE *dv;
DOUBLE *dvp;
DOUBLE *vp;
DOUBLE *yp;
DOUBLE *vn;
DOUBLE *delta;
DOUBLE *sdelta;
DOUBLE *adamy0;

DOUBLE ccor;
DOUBLE hstep; /* step size (constant) */
int order; /* Order of the prediction formula */
int ordp1;
int asiz;
int dsiz;
int ineq;
long jstep; /* counts steps taken */

//------------------------------------------------------------------------------

/* Initialize pointers in work array, compute derivatives at
 * initial position, and start the difference tables.

 * neq >= nequat in adstep() below.  If neq > nequat, unused space
 * is left for extra variables that are not actually integrated.
 */

//------------------------------------------------------------------------------

adstart ( h, yn, work, neq, ord, t )

  DOUBLE h;
  DOUBLE yn[], work[];
  int    neq, ord;
  DOUBLE t; 
{

  DOUBLE *p;
  int j;

  hstep = h;
  ineq  = neq;
  /*ccor = hstep * ((DOUBLE *)&precof) [ord];*/
  ccor = hstep * ((DOUBLE *)precof) [ord];
  jstep = 0;
  dsiz = ord + 2;
  asiz = neq * dsiz;
  order = ord;
  ordp1 = ord + 1;

  setptrs (work);

  func (t, yn, vn);

  p = dv;

  for (j = 0; j < neq; j++)
  {
    dupdate (p, 0, vn[j]);
    p += dsiz;
  }

}
//------------------------------------------------------------------------------
/* Calculate and initialize array pointers
 */
//------------------------------------------------------------------------------
setptrs ( w )

  DOUBLE *w;
{

  DOUBLE *p;

  p = w;
  dv = p;
  p += asiz;
  dvp = p;
  p += dsiz;
  vp = p;
  p += ineq;
  yp = p;
  p += ineq;
  vn = p;
  p += ineq;
  delta = p;
  p += ineq;
  sdelta = p;
  p += ineq;
  adamy0 = p;

}
//------------------------------------------------------------------------------
/* Subroutines to save and restore the internal
 * state of the integrator.
 */
//------------------------------------------------------------------------------

sinternal(f) /* save state */

  register int f; /* file */
{
  write( f, &ineq, sizeof(int) );
  write( f, &order, sizeof(int) );
  write( f, &ordp1, sizeof(int) );
  write( f, &asiz, sizeof(int) );
  write( f, &dsiz, sizeof(int) );
  write( f, &jstep, sizeof(long) );
  write( f, &ccor, sizeof(DOUBLE) );
  write( f, &hstep, sizeof(DOUBLE) );
}

rinternal(f, work) /* restore state */
  register int f;
    DOUBLE *work;
{
  read( f, &ineq, sizeof(int) );
  read( f, &order, sizeof(int) );
  read( f, &ordp1, sizeof(int) );
  read( f, &asiz, sizeof(int) );
  read( f, &dsiz, sizeof(int) );
  read( f, &jstep, sizeof(long) );
  read( f, &ccor, sizeof(DOUBLE) );
  read( f, &hstep, sizeof(DOUBLE) );
  setptrs( work );
}

//------------------------------------------------------------------------------
//
/* Adams-Bashforth-Moulton step.
 */
//------------------------------------------------------------------------------

DOUBLE adstep ( t, yn, nequat )

  DOUBLE *t;
  DOUBLE yn[];
  int nequat;
{

  DOUBLE e, e0, time;
  DOUBLE *pdv;
  int i, j;
  DOUBLE intpol();

  time = *t;
  jstep += 1;

  /* Do Runge-Kutta for the first ord + 1 steps.
    */

  if( jstep <= ordp1 )
  {
    /* Take two half steps */
    e0 = Half * hstep;
    rungek(  nequat, time, yn, e0, yn, delta );

    e = time + e0;
    func( e, yn, vn );
    rungek(  nequat, e, yn, e0, yn, delta );
    func( time+hstep, yn, vn );

    pdv = dv;

    for( j=0; j<nequat; j++ )
    {
      dupdate( pdv, (int )jstep, vn[j] );
      pdv += dsiz;
    }

    e = Zero;
    goto done;
  }



  /* Predict the next position
    * based on current y' difference table dv[].
    */
  pdv = dv;
  for (i=0; i<nequat; i++)
  {
    yp[i] = yn[i] + hstep * intpol( pdv, precof, order );
    pdv += dsiz;
  }

  /* Evaluate derivatives at the predicted position
    */
  func (time+hstep, yp, vp);

  /* Correct the predicted position and velocity using the derivatives
    * evalutated at yp.
    */

  pdv = dv;
  for( i=0; i<nequat; i++ )
  {
    vcopy( pdv, dvp, dsiz ); /* y' difference table */
    dupdate( dvp, ordp1, vp[i] );
    /* Note, the following line is equivalent to:
      *		yn[i] = yn[i] + hstep * intpol( dap, corcof, ordp1 );
      */
    yn[i] = yp[i] + ccor * dvp[order];
    pdv += dsiz;
  }

  /* Evaluate derivative at the final position
    */
  func( time+hstep, yn, vn );

  e = Zero;
  pdv = dv;
  for( i=0; i<nequat; i++ )
  {
    dupdate( pdv, ordp1, vn[i] );
    /* Note, the following line is equivalent to:
      *		yn[i] = yn[i] + hstep * intpol( dap, corcof, ordp1 );
      */
    yn[i] = yp[i] + ccor * pdv[order];
    /* Estimate the error
      */
    /*		e0 = hstep * pdv[order] * ((DOUBLE *)&corcof) [order];*/
    e0 = hstep * pdv[order] * ((DOUBLE *)corcof) [order];
    if( e0 < 0 )
      e0 = -e0;
    if( e0 > e )
      e = e0;
    pdv += dsiz;
  }

done:
  time += hstep;
  *t = time;

  return ( e );
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
