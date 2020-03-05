#ifndef NOINCS
#include "mconf.h"
#endif

#if UNK
/* almost 2^16384 */
long double MAXNUML = 1.189731495357231765021263853E4932L;
/* 2^-64 */
long double MACHEPL = 5.42101086242752217003726400434970855712890625E-20L;
/* log( MAXNUML ) */
long double MAXLOGL =  1.1356523406294143949492E4L;
/* log( underflow threshold = 2^(-16382) ) */
long double MINLOGL = -1.1355137111933024058873E4L;
long double LOGE2L  = 6.9314718055994530941723E-1L;
long double LOG2EL  = 1.4426950408889634073599E0L;
long double PIL     = 3.1415926535897932384626L;
long double PIO2L   = 1.5707963267948966192313L;
long double PIO4L   = 7.8539816339744830961566E-1L;
long double Zero = 0.0L;
long double Fourth = 0.25L;
long double Half = 0.5L;
long double One = 1.0L;
long double OneandaHalf = 1.5L;
long double Two = 2.0L;
long double Three = 3.0L;
long double ThreeandaHalf = 3.5L;
long double Four = 4.0L;
long double Five = 5.0L;
long double Six = 6.0L;
long double Ten = 10.0L;
#endif
#if IBMPC
short MAXNUML[] = {0xffff,0xffff,0xffff,0xffff,0x7ffe, XPD};
short MAXLOGL[] = {0x79ab,0xd1cf,0x17f7,0xb172,0x400c, XPD};
short MINLOGL[] = {0xeb2f,0x1210,0x8c67,0xb16c,0xc00c, XPD};
short MACHEPL[] = {0x0000,0x0000,0x0000,0x8000,0x3fbf, XPD};
short LOGE2L[]  = {0x79ac,0xd1cf,0x17f7,0xb172,0x3ffe, XPD};
short LOG2EL[]  = {0xf0bc,0x5c17,0x3b29,0xb8aa,0x3fff, XPD};
short PIL[]     = {0xc235,0x2168,0xdaa2,0xc90f,0x4000, XPD};
short PIO2L[]   = {0xc235,0x2168,0xdaa2,0xc90f,0x3fff, XPD};
short PIO4L[]   = {0xc235,0x2168,0xdaa2,0xc90f,0x3ffe, XPD};
short Zero[]    = {0x0000,0x0000,0x0000,0x0000,0x0000, XPD};
short Fourth[]  = {0x0000,0x0000,0x0000,0x8000,0x3ffd, XPD};
short Half[]    = {0x0000,0x0000,0x0000,0x8000,0x3ffe, XPD};
short One[]     = {0x0000,0x0000,0x0000,0x8000,0x3fff, XPD};
short OneandaHalf[] = {0x0000,0x0000,0x0000,0xc000,0x3fff, XPD};
short Two[]     = {0x0000,0x0000,0x0000,0x8000,0x4000, XPD};
short Three[]   = {0x0000,0x0000,0x0000,0xc000,0x4000, XPD};
short ThreeandaHalf[]   = {0x0000,0x0000,0x0000,0xe000,0x4000, XPD};
short Four[]    = {0x0000,0x0000,0x0000,0x8000,0x4001, XPD};
short Five[]    = {0x0000,0x0000,0x0000,0xa000,0x4001, XPD};
short Six[]     = {0x0000,0x0000,0x0000,0xc000,0x4001, XPD};
short Ten[]     = {0x0000,0x0000,0x0000,0xa000,0x4002, XPD};
#endif
#if MIEEE
long MAXNUML[] = {0x7ffe0000,0xffffffff,0xffffffff};
long MAXLOGL[] = {0x400c0000,0xb17217f7,0xd1cf79ab};
long MINLOGL[] = {0xc00c0000,0xb16c8c67,0x1210eb2f};
long MACHEPL[] = {0x3fbf0000,0x00000000,0x00000000};
long LOGE2L[]  = {0x3ffe0000,0xb17217f7,0xd1cf79ac};
long LOG2EL[]  = {0x3fff0000,0xb8aa3b29,0x5c17f0bc};
long PIL[]     = {0x40000000,0xc90fdaa2,0x2168c235};
long PIO2L[]   = {0x3fff0000,0xc90fdaa2,0x2168c235};
long PIO4L[]   = {0x3ffe0000,0xc90fdaa2,0x2168c235};
long Zero[]    = {0x00000000,0x00000000,0x00000000};
long Fourth[]  = {0x3ffd0000,0x80000000,0x00000000};
long Half[]    = {0x3ffe0000,0x80000000,0x00000000};
long One[]     = {0x3fff0000,0x80000000,0x00000000};
long OneandaHalf[] = {0x3fff0000,0xc0000000,0x00000000};
long Two[]     = {0x40000000,0x80000000,0x00000000};
long Three[]   = {0x40000000,0xc0000000,0x00000000};
long ThreeandaHalf[] = {0x40000000,0xe0000000,0x00000000};
long Four[]    = {0x40010000,0x80000000,0x00000000};
long Five[]    = {0x40010000,0xa0000000,0x00000000};
long Six[]     = {0x40010000,0xc0000000,0x00000000};
long Ten[]     = {0x40020000,0xa0000000,0x00000000};
#endif



/* Polynomial evaluator:
 *  P[0] x^n  +  P[1] x^(n-1)  +  ...  +  P[n]
 */
long double polevll( x, P, n )
long double x;
long double *P;
int n;
{
register long double y;

y = *P++;
do
	{
	y = y * x + *P++;
	}
while( --n );
return(y);
}



/* Polynomial evaluator:
 *  x^n  +  P[0] x^(n-1)  +  P[1] x^(n-2)  +  ...  +  P[n]
 */
long double p1evll( x, P, n )
long double x;
long double *P;
int n;
{
register long double y;

n -= 1;
y = x + *P++;
do
	{
	y = y * x + *P++;
	}
while( --n );
return( y );
}
