/*							ieee.c
 *
 *    Extended precision IEEE binary floating point arithmetic routines
 *
 * Numbers are stored in C language as arrays of 16-bit unsigned
 * short integers.  The arguments of the routines are pointers to
 * the arrays.
 *
 *
 * External e type data structure, simulates Intel 8087 chip
 * temporary real format but possibly with a larger significand:
 *
 *	NE-1 significand words	(least significant word first,
 *				 most significant bit is normally set)
 *	exponent		(value = EXONE for 1.0,
 *				top bit is the sign)
 *
 *
 * Internal data structure of a number (a "word" is 16 bits):
 *
 * ei[0]	sign word	(0 for positive, 0xffff for negative)
 * ei[1]	biased exponent	(value = EXONE for the number 1.0)
 * ei[2]	high guard word	(always zero after normalization)
 * ei[3]
 * to ei[NI-2]	significand	(NI-4 significand words,
 *				 most significant word first,
 *				 most significant bit is set)
 * ei[NI-1]	low guard word	(0x8000 bit is rounding place)
 *
 *
 *
 *		Routines for external format numbers
 *
 *	asctoe( string, e )	ASCII string to extended double e type
 *	asctoe64( string, &d )	ASCII string to long double
 *	asctoe53( string, &d )	ASCII string to double
 *	asctoe24( string, &f )	ASCII string to single
 *	asctoeg( string, e, prec ) ASCII string to specified precision
 *	e24toe( &f, e )		IEEE single precision to e type
 *	e53toe( &d, e )		IEEE double precision to e type
 *	e64toe( &d, e )		IEEE long double precision to e type
 *	eabs(e)			absolute value
 *	eadd( a, b, c )		c = b + a
 *	eclear(e)		e = 0
 *	ecmp( a, b )		compare a to b, return 1, 0, or -1
 *	ediv( a, b, c )		c = b / a
 *	efloor( a, b )		truncate to integer, toward -infinity
 *	efrexp( a, exp, s )	extract exponent and significand
 *	eifrac( e, &l, frac )   e to long integer and e type fraction
 *	euifrac( e, &l, frac )  e to unsigned long integer and e type fraction
 *	einfin( e )		set e to infinity, leaving its sign alone
 *	eldexp( a, n, b )	multiply by 2**n
 *	emov( a, b )		b = a
 *	emul( a, b, c )		c = b * a
 *	eneg(e)			e = -e
 *	eround( a, b )		b = nearest integer value to a
 *	esub( a, b, c )		c = b - a
 *	e24toasc( &f, str, n )	single to ASCII string, n digits after decimal
 *	e53toasc( &d, str, n )	double to ASCII string, n digits after decimal
 *	e64toasc( &d, str, n )	long double to ASCII string
 *	etoasc( e, str, n )	e to ASCII string, n digits after decimal
 *	etoe24( e, &f )		convert e type to IEEE single precision
 *	etoe53( e, &d )		convert e type to IEEE double precision
 *	etoe64( e, &d )		convert e type to IEEE long double precision
 *	ltoe( &l, e )		long (32 bit) integer to e type
 *	ultoe( &l, e )		unsigned long (32 bit) integer to e type
 *      eisneg( e )             1 if sign bit of e != 0, else 0
 *      eisinf( e )             1 if e has maximum exponent
 *
 *
 *		Routines for internal format numbers
 *
 *	eaddm( ai, bi )		add significands, bi = bi + ai
 *	ecleaz(ei)		ei = 0
 *	ecleazs(ei)		set ei = 0 but leave its sign alone
 *	ecmpm( ai, bi )		compare significands, return 1, 0, or -1
 *	edivm( ai, bi )		divide  significands, bi = bi / ai
 *	emdnorm(ai,l,s,exp)	normalize and round off
 *	emovi( a, ai )		convert external a to internal ai
 *	emovo( ai, a )		convert internal ai to external a
 *	emovz( ai, bi )		bi = ai, low guard word of bi = 0
 *	emulm( ai, bi )		multiply significands, bi = bi * ai
 *	enormlz(ei)		left-justify the significand
 *	eshdn1( ai )		shift significand and guards down 1 bit
 *	eshdn8( ai )		shift down 8 bits
 *	eshdn6( ai )		shift down 16 bits
 *	eshift( ai, n )		shift ai n bits up (or down if n < 0)
 *	eshup1( ai )		shift significand and guards up 1 bit
 *	eshup8( ai )		shift up 8 bits
 *	eshup6( ai )		shift up 16 bits
 *	esubm( ai, bi )		subtract significands, bi = bi - ai
 *
 *
 * The result is always normalized and rounded to NI-4 word precision
 * after each arithmetic operation.
 *
 * Exception flags and NaNs are NOT fully supported.
 * This arithmetic should never produce a NaN output, but it might
 * be confused by a NaN input.
 * Define INFINITY in mconf.h for support of infinity; otherwise a
 * saturation arithmetic is implemented.
 * Denormals are always supported here where appropriate (e.g., not
 * for conversion to DEC numbers).
 */

/*
 * Revision history:
 *
 *  5 Jan 84	PDP-11 assembly language version
 *  2 Mar 86	fixed bug in asctoq()
 *  6 Dec 86	C language version
 * 30 Aug 88	100 digit version, improved rounding
 * 15 May 92    80-bit long double support
 *
 * Author:  S. L. Moshier.
 */

#include <stdio.h>
#include "ehead.h"
#include "mconf.h"

/* Turn UNK into something else */
#if UNK
#if BIGENDIAN
#undef MIEEE
#define MIEEE 1
#else
#undef IBMPC
#define IBMPC 1
#endif
#undef UNK
#endif


/* Control register for rounding precision.
 * This can be set to 80 (if NE=6), 64, 56, 53, or 24 bits.
 */
int rndprc = NBITS;
extern int rndprc;

void eaddm(), esubm(), emdnorm(), asctoeg();
static void toe24(), toe53(), toe64();
void eremain(), einit(), eiremain();
int ecmpm(), edivm(), emulm(), eisneg(), eisinf();
void emovi(), emovo(), emovz(), ecleaz(), eadd1();
void etodec(), todec(), dectoe();




void einit()
{
}

/*
; Clear out entire external format number.
;
; unsigned short x[];
; eclear( x );
*/

void eclear( x )
register unsigned short *x;
{
register int i;

for( i=0; i<NE; i++ )
	*x++ = 0;
}



/* Move external format number from a to b.
 *
 * emov( a, b );
 */

void emov( a, b )
register unsigned short *a, *b;
{
register int i;

for( i=0; i<NE; i++ )
	*b++ = *a++;
}


/*
;	Absolute value of external format number
;
;	short x[NE];
;	eabs( x );
*/

void eabs(x)
unsigned short x[];	/* x is the memory address of a short */
{

x[NE-1] &= 0x7fff; /* sign is top bit of last word of external format */
}




/*
;	Negate external format number
;
;	unsigned short x[NE];
;	eneg( x );
*/

void eneg(x)
unsigned short x[];
{

x[NE-1] ^= 0x8000; /* Toggle the sign bit */
}



/* Return 1 if external format number is negative,
 * else return zero.
 */
int eisneg(x)
unsigned short x[];
{

if( x[NE-1] & 0x8000 )
	return( 1 );
else
	return( 0 );
}


/* Return 1 if external format number has maximum possible exponent,
 * else return zero.
 */
int eisinf(x)
unsigned short x[];
{

if( (x[NE-1] & 0x7fff) == 0x7fff )
	return( 1 );
else
	return( 0 );
}


/*
; Fill entire number, including exponent and significand, with
; largest possible number.  These programs implement a saturation
; value that is an ordinary, legal number.  A special value
; "infinity" may also be implemented; this would require tests
; for that value and implementation of special rules for arithmetic
; operations involving inifinity.
*/

void einfin(x)
register unsigned short *x;
{
register int i;

#ifdef INFINITY
for( i=0; i<NE-1; i++ )
	*x++ = 0;
*x |= 32767;
#else
for( i=0; i<NE-1; i++ )
	*x++ = 0xffff;
*x |= 32766;
if( rndprc < NBITS )
	{
	if( rndprc == 64 )
		{
		*(x-5) = 0;
		}
	if( rndprc == 53 )
		{
		*(x-4) = 0xf800;
		}
	else
		{
		*(x-4) = 0;
		*(x-3) = 0;
		*(x-2) = 0xff00;
		}
	}
#endif
}



/* Move in external format number,
 * converting it to internal format.
 */
void emovi( a, b )
unsigned short *a, *b;
{
register unsigned short *p, *q;
int i;

q = b;
p = a + (NE-1);	/* point to last word of external number */
/* get the sign bit */
if( *p & 0x8000 )
	*q++ = 0xffff;
else
	*q++ = 0;
/* get the exponent */
*q = *p--;
*q++ &= 0x7fff;	/* delete the sign bit */
#ifdef INFINITY
if( (*(q-1) & 0x7fff) == 0x7fff )
	{
	for( i=2; i<NI; i++ )
		*q++ = 0;
	return;
	}
#endif
/* clear high guard word */
*q++ = 0;
/* move in the significand */
for( i=0; i<NE-1; i++ )
	*q++ = *p--;
/* clear low guard word */
*q = 0;
}


/* Move internal format number out,
 * converting it to external format.
 */
void emovo( a, b )
unsigned short *a, *b;
{
register unsigned short *p, *q;
unsigned short i;

p = a;
q = b + (NE-1); /* point to output exponent */
/* combine sign and exponent */
i = *p++;
if( i )
	*q-- = *p++ | 0x8000;
else
	*q-- = *p++;
#ifdef INFINITY
if( *(p-1) == 0x7fff )
	{
	einfin(b);
	return;
	}
#endif
/* skip over guard word */
++p;
/* move the significand */
for( i=0; i<NE-1; i++ )
	*q-- = *p++;
}




/* Clear out internal format number.
 */

void ecleaz( xi )
register unsigned short *xi;
{
register int i;

for( i=0; i<NI; i++ )
	*xi++ = 0;
}

/* same, but don't touch the sign. */

void ecleazs( xi )
register unsigned short *xi;
{
register int i;

++xi;
for(i=0; i<NI-1; i++)
	*xi++ = 0;
}




/* Move internal format number from a to b.
 */
void emovz( a, b )
register unsigned short *a, *b;
{
register int i;

for( i=0; i<NI-1; i++ )
	*b++ = *a++;
/* clear low guard word */
*b = 0;
}


/*
;	Compare significands of numbers in internal format.
;	Guard words are included in the comparison.
;
;	unsigned short a[NI], b[NI];
;	cmpm( a, b );
;
;	for the significands:
;	returns	+1 if a > b
;		 0 if a == b
;		-1 if a < b
*/
int ecmpm( a, b )
register unsigned short *a, *b;
{
int i;

a += M; /* skip up to significand area */
b += M;
for( i=M; i<NI; i++ )
	{
	if( *a++ != *b++ )
		goto difrnt;
	}
return(0);

difrnt:
if( *(--a) > *(--b) )
	return(1);
else
	return(-1);
}


/*
;	Shift significand down by 1 bit
*/

void eshdn1(x)
register unsigned short *x;
{
register unsigned short bits;
int i;

x += M;	/* point to significand area */

bits = 0;
for( i=M; i<NI; i++ )
	{
	if( *x & 1 )
		bits |= 1;
	*x >>= 1;
	if( bits & 2 )
		*x |= 0x8000;
	bits <<= 1;
	++x;
	}	
}



/*
;	Shift significand up by 1 bit
*/

void eshup1(x)
register unsigned short *x;
{
register unsigned short bits;
int i;

x += NI-1;
bits = 0;

for( i=M; i<NI; i++ )
	{
	if( *x & 0x8000 )
		bits |= 1;
	*x <<= 1;
	if( bits & 2 )
		*x |= 1;
	bits <<= 1;
	--x;
	}
}



/*
;	Shift significand down by 8 bits
*/

void eshdn8(x)
register unsigned short *x;
{
register unsigned short newbyt, oldbyt;
int i;

x += M;
oldbyt = 0;
for( i=M; i<NI; i++ )
	{
	newbyt = *x << 8;
	*x >>= 8;
	*x |= oldbyt;
	oldbyt = newbyt;
	++x;
	}
}

/*
;	Shift significand up by 8 bits
*/

void eshup8(x)
register unsigned short *x;
{
int i;
register unsigned short newbyt, oldbyt;

x += NI-1;
oldbyt = 0;

for( i=M; i<NI; i++ )
	{
	newbyt = *x >> 8;
	*x <<= 8;
	*x |= oldbyt;
	oldbyt = newbyt;
	--x;
	}
}

/*
;	Shift significand up by 16 bits
*/

void eshup6(x)
register unsigned short *x;
{
int i;
register unsigned short *p;

p = x + M;
x += M + 1;

for( i=M; i<NI-1; i++ )
	*p++ = *x++;

*p = 0;
}

/*
;	Shift significand down by 16 bits
*/

void eshdn6(x)
register unsigned short *x;
{
int i;
register unsigned short *p;

x += NI-1;
p = x + 1;

for( i=M; i<NI-1; i++ )
	*(--p) = *(--x);

*(--p) = 0;
}

/*
;	Add significands
;	x + y replaces y
*/

void eaddm( x, y )
unsigned short *x, *y;
{
register unsigned long a;
int i;
unsigned int carry;

x += NI-1;
y += NI-1;
carry = 0;
for( i=M; i<NI; i++ )
	{
	a = (unsigned long )(*x) + (unsigned long )(*y) + carry;
	if( a & 0x10000 )
		carry = 1;
	else
		carry = 0;
	*y = (unsigned short )a;
	--x;
	--y;
	}
}

/*
;	Subtract significands
;	y - x replaces y
*/

void esubm( x, y )
unsigned short *x, *y;
{
unsigned long a;
int i;
unsigned int carry;

x += NI-1;
y += NI-1;
carry = 0;
for( i=M; i<NI; i++ )
	{
	a = (unsigned long )(*y) - (unsigned long )(*x) - carry;
	if( a & 0x10000 )
		carry = 1;
	else
		carry = 0;
	*y = (unsigned short )a;
	--x;
	--y;
	}
}


/* Divide significands */

static unsigned short equot[NI] = {0};

int edivm( den, num )
unsigned short den[], num[];
{
int i;
register unsigned short *p, *q;
unsigned short j;

p = &equot[0];
*p++ = num[0];
*p++ = num[1];

for( i=M; i<NI; i++ )
	{
	*p++ = 0;
	}

/* Use faster compare and subtraction if denominator
 * has only 15 bits of significance.
 */
p = &den[M+2];
if( *p++ == 0 )
	{
	for( i=M+3; i<NI; i++ )
		{
		if( *p++ != 0 )
			goto fulldiv;
		}
	if( (den[M+1] & 1) != 0 )
		goto fulldiv;
	eshdn1(num);
	eshdn1(den);

	p = &den[M+1];
	q = &num[M+1];

	for( i=0; i<NBITS+2; i++ )
		{
		if( *p <= *q )
			{
			*q -= *p;
			j = 1;
			}
		else
			{
			j = 0;
			}
		eshup1(equot);
		equot[NI-2] |= j;
		eshup1(num);
		}
	goto divdon;
	}

/* The number of quotient bits to calculate is
 * NBITS + 1 scaling guard bit + 1 roundoff bit.
 */
fulldiv:

p = &equot[NI-2];
for( i=0; i<NBITS+2; i++ )
	{
	if( ecmpm(den,num) <= 0 )
		{
		esubm(den, num);
		j = 1;	/* quotient bit = 1 */
		}
	else
		j = 0;
	eshup1(equot);
	*p |= j;
	eshup1(num);
	}

divdon:

eshdn1( equot );
eshdn1( equot );

/* test for nonzero remainder after roundoff bit */
p = &num[M];
j = 0;
for( i=M; i<NI; i++ )
	{
	j |= *p++;
	}
if( j )
	j = 1;


for( i=0; i<NI; i++ )
	num[i] = equot[i];
return( (int )j );
}


/* Multiply significands */
int emulm( a, b )
unsigned short a[], b[];
{
unsigned short *p, *q;
int i, j, k;

equot[0] = b[0];
equot[1] = b[1];
for( i=M; i<NI; i++ )
	equot[i] = 0;

p = &a[NI-2];
k = NBITS;
while( *p == 0 ) /* significand is not supposed to be all zero */
	{
	eshdn6(a);
	k -= 16;
	}
if( (*p & 0xff) == 0 )
	{
	eshdn8(a);
	k -= 8;
	}

q = &equot[NI-1];
j = 0;
for( i=0; i<k; i++ )
	{
	if( *p & 1 )
		eaddm(b, equot);
/* remember if there were any nonzero bits shifted out */
	if( *q & 1 )
		j |= 1;
	eshdn1(a);
	eshdn1(equot);
	}

for( i=0; i<NI; i++ )
	b[i] = equot[i];

/* return flag for lost nonzero bits */
return(j);
}



/*
 * Normalize and round off.
 *
 * The internal format number to be rounded is "s".
 * Input "lost" indicates whether or not the number is exact.
 * This is the so-called sticky bit.
 *
 * Input "subflg" indicates whether the number was obtained
 * by a subtraction operation.  In that case if lost is nonzero
 * then the number is slightly smaller than indicated.
 *
 * Input "exp" is the biased exponent, which may be negative.
 * the exponent field of "s" is ignored but is replaced by
 * "exp" as adjusted by normalization and rounding.
 *
 * Input "rcntrl" is the rounding control.
 */

static int rlast = -1;
static int rw = 0;
static unsigned short rmsk = 0;
static unsigned short rmbit = 0;
static unsigned short rebit = 0;
static int re = 0;
static unsigned short rbit[NI] = {0,0,0,0,0,0,0,0};

void emdnorm( s, lost, subflg, exp, rcntrl )
unsigned short s[];
int lost;
int subflg;
long exp;
int rcntrl;
{
int i, j;
unsigned short r;

/* Normalize */
j = enormlz( s );

/* a blank significand could mean either zero or infinity. */
#ifndef INFINITY
if( j > NBITS )
	{
	ecleazs( s );
	return;
	}
#endif
exp -= j;
#ifndef INFINITY
if( exp >= 32767L )
	goto overf;
#endif
if( exp < 0L )
	{
	if( exp > (long )(-NBITS-1) )
		{
		j = (int )exp;
		i = eshift( s, j );
		if( i )
			lost = 1;
		}
	else
		{
		ecleazs( s );
		return;
		}
	}
/* Round off, unless told not to by rcntrl. */
if( rcntrl == 0 )
	goto mdfin;
/* Set up rounding parameters if the control register changed. */
if( rndprc != rlast )
	{
	ecleaz( rbit );
	switch( rndprc )
		{
		default:
		case NBITS:
			rw = NI-1; /* low guard word */
			rmsk = 0xffff;
			rmbit = 0x8000;
			rbit[rw-1] = 1;
			re = NI-2;
			rebit = 1;
			break;
		case 64:
			rw = 7;
			rmsk = 0xffff;
			rmbit = 0x8000;
			rbit[rw-1] = 1;
			re = rw-1;
			rebit = 1;
			break;
/* For DEC arithmetic */
		case 56:
			rw = 6;
			rmsk = 0xff;
			rmbit = 0x80;
			rbit[rw] = 0x100;
			re = rw;
			rebit = 0x100;
			break;
		case 53:
			rw = 6;
			rmsk = 0x7ff;
			rmbit = 0x0400;
			rbit[rw] = 0x800;
			re = rw;
			rebit = 0x800;
			break;
		case 24:
			rw = 4;
			rmsk = 0xff;
			rmbit = 0x80;
			rbit[rw] = 0x100;
			re = rw;
			rebit = 0x100;
			break;
		}
	rlast = rndprc;
	}

if( rndprc >= 64 )
	{
	r = s[rw] & rmsk;
	if( rndprc == 64 )
		{
		i = rw+1;
		while( i < NI )
			{
			if( s[i] )
				r |= 1;
			s[i] = 0;
			++i;
			}
		}
	}
else
	{
	if( exp <= 0 )
		eshdn1(s);
	r = s[rw] & rmsk;
/* These tests assume NI = 8 */
	i = rw+1;
	while( i < NI )
		{
		if( s[i] )
			r |= 1;
		s[i] = 0;
		++i;
		}
/*
	if( rndprc == 24 )
		{
		if( s[5] || s[6] )
			r |= 1;
		s[5] = 0;
		s[6] = 0;
		}
*/
	s[rw] &= ~rmsk;
	}
if( (r & rmbit) != 0 )
	{
	if( r == rmbit)
		{
		if( lost == 0 )
			{ /* round to even */
			if( (s[re] & rebit) == 0 )
				goto mddone;
			}
		else
			{
			if( subflg != 0 )
				goto mddone;
			}
		}
	eaddm( rbit, s );
	}
mddone:
if( (rndprc < 64) && (exp <= 0) )
	{
	eshup1(s);
	}
if( s[2] != 0 )
	{ /* overflow on roundoff */
	eshdn1(s);
	exp += 1;
	}
mdfin:
s[NI-1] = 0;
if( exp >= 32767L )
	{
#ifndef INFINITY
overf:
#endif
#ifdef INFINITY
	s[1] = 32767;
	for( i=2; i<NI-1; i++ )
		s[i] = 0;
#else
	s[1] = 32766;
	s[2] = 0;
	for( i=M+1; i<NI-1; i++ )
		s[i] = 0xffff;
	s[NI-1] = 0;
	if( rndprc < 64 )
		{
		s[rw] &= ~rmsk;
		if( rndprc == 24 )
			{
			s[5] = 0;
			s[6] = 0;
			}
		}
#endif
	return;
	}
if( exp < 0 )
	s[1] = 0;
else
	s[1] = (unsigned short )exp;
}



/*
;	Subtract external format numbers.
;
;	unsigned short a[NE], b[NE], c[NE];
;	esub( a, b, c );	 c = b - a
*/

static int subflg = 0;

void esub( a, b, c )
unsigned short *a, *b, *c;
{

subflg = 1;
eadd1( a, b, c );
}


/*
;	Add.
;
;	unsigned short a[NE], b[NE], c[NE];
;	eadd( a, b, c );	 c = b + a
*/
void eadd( a, b, c )
unsigned short *a, *b, *c;
{

subflg = 0;
eadd1( a, b, c );
}

void eadd1( a, b, c )
unsigned short *a, *b, *c;
{
unsigned short ai[NI], bi[NI], ci[NI];
int i, lost, j, k;
long lt, lta, ltb;

#ifdef INFINITY
if( eisinf(a) )
	{
	emov(a,c);
	if( subflg )
		eneg(c);
	return;
	}
if( eisinf(b) )
	{
	emov(b,c);
	return;
	}
#endif
emovi( a, ai );
emovi( b, bi );
if( subflg )
	ai[0] = ~ai[0];

/* compare exponents */
lta = ai[E];
ltb = bi[E];
lt = lta - ltb;
if( lt > 0L )
	{	/* put the larger number in bi */
	emovz( bi, ci );
	emovz( ai, bi );
	emovz( ci, ai );
	ltb = bi[E];
	lt = -lt;
	}
lost = 0;
if( lt != 0L )
	{
	if( lt < (long )(-NBITS-1) )
		goto done;	/* answer same as larger addend */
	k = (int )lt;
	lost = eshift( ai, k ); /* shift the smaller number down */
	}
else
	{
/* exponents were the same, so must compare significands */
	i = ecmpm( ai, bi );
	if( i == 0 )
		{ /* the numbers are identical in magnitude */
		/* if different signs, result is zero */
		if( ai[0] != bi[0] )
			{
			eclear(c);
			return;
			}
		/* if same sign, result is double */
		/* double denomalized tiny number */
		if( (bi[E] == 0) && ((bi[3] & 0x8000) == 0) )
			{
			eshup1( bi );
			goto done;
			}
		/* add 1 to exponent unless both are zero! */
		for( j=1; j<NI-1; j++ )
			{
			if( bi[j] != 0 )
				{
/* This could overflow, but let emovo take care of that. */
				ltb += 1;
				break;
				}
			}
		bi[E] = (unsigned short )ltb;
		goto done;
		}
	if( i > 0 )
		{	/* put the larger number in bi */
		emovz( bi, ci );
		emovz( ai, bi );
		emovz( ci, ai );
		}
	}
if( ai[0] == bi[0] )
	{
	eaddm( ai, bi );
	subflg = 0;
	}
else
	{
	esubm( ai, bi );
	subflg = 1;
	}
emdnorm( bi, lost, subflg, ltb, 64 );

done:
emovo( bi, c );
}



/*
;	Divide.
;
;	unsigned short a[NE], b[NE], c[NE];
;	ediv( a, b, c );	c = b / a
*/
void ediv( a, b, c )
unsigned short *a, *b, *c;
{
unsigned short ai[NI], bi[NI];
int i;
long lt, lta, ltb;

#ifdef INFINITY
if( eisinf(b) )
	{
	if( eisneg(a) ^ eisneg(b) )
		*(c+(NE-1)) = 0x8000;
	else
		*(c+(NE-1)) = 0;
	einfin(c);
	return;
	}
if( eisinf(a) )
	{
	eclear(c);
	return;
	}
#endif
emovi( a, ai );
emovi( b, bi );
lta = ai[E];
ltb = bi[E];
if( bi[E] == 0 )
	{ /* See if numerator is zero. */
	for( i=1; i<NI-1; i++ )
		{
		if( bi[i] != 0 )
			{
			ltb -= enormlz( bi );
			goto dnzro1;
			}
		}
	eclear(c);
	return;
	}
dnzro1:

if( ai[E] == 0 )
	{	/* possible divide by zero */
	for( i=1; i<NI-1; i++ )
		{
		if( ai[i] != 0 )
			{
			lta -= enormlz( ai );
			goto dnzro2;
			}
		}
	if( ai[0] == bi[0] )
		*(c+(NE-1)) = 0;
	else
		*(c+(NE-1)) = 0x8000;
	einfin(c);
	mtherr( "ediv", SING );
	return;
	}
dnzro2:

i = edivm( ai, bi );
/* calculate exponent */
lt = ltb - lta + EXONE;
emdnorm( bi, i, 0, lt, 64 );
/* set the sign */
if( ai[0] == bi[0] )
	bi[0] = 0;
else
	bi[0] = 0Xffff;
emovo( bi, c );
}



/*
;	Multiply.
;
;	unsigned short a[NE], b[NE], c[NE];
;	emul( a, b, c );	c = b * a
*/
void emul( a, b, c )
unsigned short *a, *b, *c;
{
unsigned short ai[NI], bi[NI];
int i, j;
long lt, lta, ltb;

#ifdef INFINITY
if( eisinf(a) || eisinf(b) )
	{
	if( eisneg(a) ^ eisneg(b) )
		*(c+(NE-1)) = 0x8000;
	else
		*(c+(NE-1)) = 0;
	einfin(c);
	return;
	}
#endif
emovi( a, ai );
emovi( b, bi );
lta = ai[E];
ltb = bi[E];
if( ai[E] == 0 )
	{
	for( i=1; i<NI-1; i++ )
		{
		if( ai[i] != 0 )
			{
			lta -= enormlz( ai );
			goto mnzer1;
			}
		}
	eclear(c);
	return;
	}
mnzer1:

if( bi[E] == 0 )
	{
	for( i=1; i<NI-1; i++ )
		{
		if( bi[i] != 0 )
			{
			ltb -= enormlz( bi );
			goto mnzer2;
			}
		}
	eclear(c);
	return;
	}
mnzer2:

/* Multiply significands */
j = emulm( ai, bi );
/* calculate exponent */
lt = lta + ltb - (EXONE - 1);
emdnorm( bi, j, 0, lt, 64 );
/* calculate sign of product */
if( ai[0] == bi[0] )
	bi[0] = 0;
else
	bi[0] = 0xffff;
emovo( bi, c );
}




/*
; Convert IEEE double precision to e type
;	double d;
;	unsigned short x[N+2];
;	e53toe( &d, x );
*/
void e53toe( e, y )
unsigned short *e, *y;
{
#ifdef DEC

dectoe( e, y ); /* see etodec.c */

#else

register unsigned short r;
register unsigned short *p;
unsigned short yy[NI];
int denorm, k;

denorm = 0;	/* flag if denormalized number */
ecleaz(yy);
#ifdef IBMPC
e += 3;
#endif
r = *e;
yy[0] = 0;
if( r & 0x8000 )
	yy[0] = 0xffff;
yy[M] = (r & 0x0f) | 0x10;
r &= ~0x800f;	/* strip sign and 4 significand bits */
#ifdef INFINITY
if( r == 0x7ff0 )
	{
	einfin( y );
	if( r & 0x8000 )
		eneg(y);
	return;
	}
#endif
r >>= 4;
/* If zero exponent, then the significand is denormalized.
 * So, take back the understood high significand bit. */ 
if( r == 0 )
	{
	denorm = 1;
	yy[M] &= ~0x10;
	}
r += EXONE - 01777;
yy[E] = r;
p = &yy[M+1];
#ifdef IBMPC
*p++ = *(--e);
*p++ = *(--e);
*p++ = *(--e);
#endif
#ifdef MIEEE
++e;
*p++ = *e++;
*p++ = *e++;
*p++ = *e++;
#endif
(void )eshift( yy, -5 );
if( denorm )
	{ /* if zero exponent, then normalize the significand */
	if( (k = enormlz(yy)) > NBITS )
		ecleazs(yy);
	else
		yy[E] -= (unsigned short )(k-1);
	}
emovo( yy, y );
#endif /* not DEC */
}

void e64toe( e, y )
unsigned short *e, *y;
{
unsigned short yy[NI];
unsigned short *p, *q;
int i;

p = yy;
for( i=0; i<NE-5; i++ )
	*p++ = 0;
#ifdef IBMPC
for( i=0; i<5; i++ )
	*p++ = *e++;
#endif
#ifdef MIEEE
p = &yy[0] + (NE-1);
*p-- = *e++;
++e;
for( i=0; i<4; i++ )
	*p-- = *e++;
#endif
p = yy;
q = y;
#ifdef INFINITY
if( *p == 0x7fff )
	{
	einfin( y );
	if( *p & 0x8000 )
		eneg(y);
	return;
	}
#endif
for( i=0; i<NE; i++ )
	*q++ = *p++;
}


/*
; Convert IEEE single precision to e type
;	float d;
;	unsigned short x[N+2];
;	dtox( &d, x );
*/
void e24toe( e, y )
unsigned short *e, *y;
{
register unsigned short r;
register unsigned short *p;
unsigned short yy[NI];
int denorm, k;

denorm = 0;	/* flag if denormalized number */
ecleaz(yy);
#ifdef IBMPC
e += 1;
#endif
#ifdef DEC
e += 1;
#endif
r = *e;
yy[0] = 0;
if( r & 0x8000 )
	yy[0] = 0xffff;
yy[M] = (r & 0x7f) | 0200;
r &= ~0x807f;	/* strip sign and 7 significand bits */
#ifdef INFINITY
if( r == 0x7f80 )
	{
	einfin( y );
	if( r & 0x8000 )
		eneg(y);
	return;
	}
#endif
r >>= 7;
/* If zero exponent, then the significand is denormalized.
 * So, take back the understood high significand bit. */ 
if( r == 0 )
	{
	denorm = 1;
	yy[M] &= ~0200;
	}
r += EXONE - 0177;
yy[E] = r;
p = &yy[M+1];
#ifdef IBMPC
*p++ = *(--e);
#endif
#ifdef DEC
*p++ = *(--e);
#endif
#ifdef MIEEE
++e;
*p++ = *e++;
#endif
(void )eshift( yy, -8 );
if( denorm )
	{ /* if zero exponent, then normalize the significand */
	if( (k = enormlz(yy)) > NBITS )
		ecleazs(yy);
	else
		yy[E] -= (unsigned short )(k-1);
	}
emovo( yy, y );
}


void etoe64( x, e )
unsigned short *x, *e;
{
unsigned short xi[NI];
long exp;
int rndsav;

emovi( x, xi );
exp = (long )xi[E]; /* adjust exponent for offset */
#ifdef INFINITY
if( eisinf(x) )
	goto nonorm;
#endif
/* round off to nearest or even */
rndsav = rndprc;
rndprc = 64;
emdnorm( xi, 0, 0, exp, 64 );
rndprc = rndsav;
nonorm:
toe64( xi, e );
}

/* move out internal format to ieee long double */
static void toe64( a, b )
unsigned short *a, *b;
{
register unsigned short *p, *q;
unsigned short i;

p = a;
#ifdef MIEEE
q = b;
#else
q = b + 4; /* point to output exponent */
#endif

/* combine sign and exponent */
i = *p++;
#ifdef MIEEE
if( i )
	*q++ = *p++ | 0x8000;
else
	*q++ = *p++;
*q++ = 0;
#else
if( i )
	*q-- = *p++ | 0x8000;
else
	*q-- = *p++;
#endif
/* skip over guard word */
++p;
/* move the significand */
#ifdef MIEEE
for( i=0; i<4; i++ )
	*q++ = *p++;
#else
for( i=0; i<4; i++ )
	*q-- = *p++;
#endif
}


/*
; e type to IEEE double precision
;	double d;
;	unsigned short x[NE];
;	etoe53( x, &d );
*/

#ifdef DEC

void etoe53( x, e )
unsigned short *x, *e;
{
etodec( x, e ); /* see etodec.c */
}

static void toe53( x, y )
unsigned short *x, *y;
{
todec( x, y );
}

#else

void etoe53( x, e )
unsigned short *x, *e;
{
unsigned short xi[NI];
long exp;
int rndsav;

emovi( x, xi );
exp = (long )xi[E] - (EXONE - 0x3ff); /* adjust exponent for offsets */
#ifdef INFINITY
if( eisinf(x) )
	goto nonorm;
#endif
/* round off to nearest or even */
rndsav = rndprc;
rndprc = 53;
emdnorm( xi, 0, 0, exp, 64 );
rndprc = rndsav;
nonorm:
toe53( xi, e );
}


static void toe53( x, y )
unsigned short *x, *y;
{
unsigned short i;
unsigned short *p;


p = &x[0];
#ifdef IBMPC
y += 3;
#endif
*y = 0;	/* output high order */
if( *p++ )
	*y = 0x8000;	/* output sign bit */

i = *p++;
if( i >= (unsigned int )2047 )
	{	/* Saturate at largest number less than infinity. */
#ifdef INFINITY
	*y |= 0x7ff0;
#ifdef IBMPC
	*(--y) = 0;
	*(--y) = 0;
	*(--y) = 0;
#endif
#ifdef MIEEE
	++y;
	*y++ = 0;
	*y++ = 0;
	*y++ = 0;
#endif
#else
	*y |= (unsigned short )0x7fef;
#ifdef IBMPC
	*(--y) = 0xffff;
	*(--y) = 0xffff;
	*(--y) = 0xffff;
#endif
#ifdef MIEEE
	++y;
	*y++ = 0xffff;
	*y++ = 0xffff;
	*y++ = 0xffff;
#endif
#endif
	return;
	}
if( i == 0 )
	{
	(void )eshift( x, 4 );
	}
else
	{
	i <<= 4;
	(void )eshift( x, 5 );
	}
i |= *p++ & (unsigned short )0x0f;	/* *p = xi[M] */
*y |= (unsigned short )i; /* high order output already has sign bit set */
#ifdef IBMPC
*(--y) = *p++;
*(--y) = *p++;
*(--y) = *p;
#endif
#ifdef MIEEE
++y;
*y++ = *p++;
*y++ = *p++;
*y++ = *p++;
#endif
}

#endif /* not DEC */



/*
; e type to IEEE single precision
;	float d;
;	unsigned short x[N+2];
;	xtod( x, &d );
*/
void etoe24( x, e )
unsigned short *x, *e;
{
long exp;
unsigned short xi[NI];
int rndsav;

emovi( x, xi );
exp = (long )xi[E] - (EXONE - 0177); /* adjust exponent for offsets */
#ifdef INFINITY
if( eisinf(x) )
	goto nonorm;
#endif
/* round off to nearest or even */
rndsav = rndprc;
rndprc = 24;
emdnorm( xi, 0, 0, exp, 64 );
rndprc = rndsav;
nonorm:
toe24( xi, e );
}

static void toe24( x, y )
unsigned short *x, *y;
{
unsigned short i;
unsigned short *p;

p = &x[0];
#ifdef IBMPC
y += 1;
#endif
#ifdef DEC
y += 1;
#endif
*y = 0;	/* output high order */
if( *p++ )
	*y = 0x8000;	/* output sign bit */

i = *p++;
if( i >= 255 )
	{	/* Saturate at largest number less than infinity. */
#ifdef INFINITY
	*y |= (unsigned short )0x7f80;
#ifdef IBMPC
	*(--y) = 0;
#endif
#ifdef DEC
	*(--y) = 0;
#endif
#ifdef MIEEE
	++y;
	*y = 0;
#endif
#else
	*y |= (unsigned short )0x7f7f;
#ifdef IBMPC
	*(--y) = 0xffff;
#endif
#ifdef DEC
	*(--y) = 0xffff;
#endif
#ifdef MIEEE
	++y;
	*y = 0xffff;
#endif
#endif
	return;
	}
if( i == 0 )
	{
	(void )eshift( x, 7 );
	}
else
	{
	i <<= 7;
	(void )eshift( x, 8 );
	}
i |= *p++ & (unsigned short )0x7f;	/* *p = xi[M] */
*y |= i;	/* high order output already has sign bit set */
#ifdef IBMPC
*(--y) = *p;
#endif
#ifdef DEC
*(--y) = *p;
#endif
#ifdef MIEEE
++y;
*y = *p;
#endif
}


/* Compare two e type numbers.
 *
 * unsigned short a[NE], b[NE];
 * ecmp( a, b );
 *
 *  returns +1 if a > b
 *           0 if a == b
 *          -1 if a < b
 */
int ecmp( a, b )
unsigned short *a, *b;
{
unsigned short ai[NI], bi[NI];
register unsigned short *p, *q;
register int i;
int msign;

emovi( a, ai );
p = ai;
emovi( b, bi );
q = bi;

if( *p != *q )
	{ /* the signs are different */
/* -0 equals + 0 */
	for( i=1; i<NI-1; i++ )
		{
		if( ai[i] != 0 )
			goto nzro;
		if( bi[i] != 0 )
			goto nzro;
		}
	return(0);
nzro:
	if( *p == 0 )
		return( 1 );
	else
		return( -1 );
	}
/* both are the same sign */
if( *p == 0 )
	msign = 1;
else
	msign = -1;
i = NI-1;
do
	{
	if( *p++ != *q++ )
		{
		goto diff;
		}
	}
while( --i > 0 );

return(0);	/* equality */



diff:

if( *(--p) > *(--q) )
	return( msign );		/* p is bigger */
else
	return( -msign );	/* p is littler */
}




/* Find nearest integer to x = floor( x + 0.5 )
 *
 * unsigned short x[NE], y[NE]
 * eround( x, y );
 */
void eround( x, y )
unsigned short *x, *y;
{

eadd( ehalf, x, y );
efloor( y, y );
}




/*
; convert long (32-bit) integer to e type
;
;	long l;
;	unsigned short x[NE];
;	ltoe( &l, x );
; note &l is the memory address of l
*/
void ltoe( lp, y )
long *lp;	/* lp is the memory address of a long integer */
unsigned short *y;	/* y is the address of a short */
{
unsigned short yi[NI];
unsigned long ll;
int k;

ecleaz( yi );
if( *lp < 0 )
	{
	ll =  (unsigned long )( -(*lp) ); /* make it positive */
	yi[0] = 0xffff; /* put correct sign in the e type number */
	}
else
	{
	ll = (unsigned long )( *lp );
	}
/* move the long integer to yi significand area */
yi[M] = (unsigned short )(ll >> 16); 
yi[M+1] = (unsigned short )ll;

yi[E] = EXONE + 15;	/* exponent if normalize shift count were 0 */
if( (k = enormlz( yi )) > NBITS ) /* normalize the significand */
	ecleaz( yi );	/* it was zero */
else
	yi[E] -= (unsigned short )k; /* subtract shift count from exponent */
emovo( yi, y );	/* output the answer */
}

/*
; convert unsigned long (32-bit) integer to e type
;
;	unsigned long l;
;	unsigned short x[NE];
;	ltox( &l, x );
; note &l is the memory address of l
*/
void ultoe( lp, y )
unsigned long *lp; /* lp is the memory address of a long integer */
unsigned short *y;	/* y is the address of a short */
{
unsigned short yi[NI];
unsigned long ll;
int k;

ecleaz( yi );
ll = *lp;

/* move the long integer to ayi significand area */
yi[M] = (unsigned short )(ll >> 16); 
yi[M+1] = (unsigned short )ll;

yi[E] = EXONE + 15;	/* exponent if normalize shift count were 0 */
if( (k = enormlz( yi )) > NBITS ) /* normalize the significand */
	ecleaz( yi );	/* it was zero */
else
	yi[E] -= (unsigned short )k; /* subtract shift count from exponent */
emovo( yi, y );	/* output the answer */
}


/*
;	Find long integer and fractional parts

;	long i;
;	unsigned short x[NE], frac[NE];
;	xifrac( x, &i, frac );
 
  The integer output has the sign of the input.  The fraction is
  the positive fractional part of abs(x).
*/
void eifrac( x, i, frac )
unsigned short *x;
long *i;
unsigned short *frac;
{
unsigned short xi[NI];
int k;

emovi( x, xi );
k = (int )xi[E] - (EXONE - 1);
if( k <= 0 )
	{
/* if exponent <= 0, integer = 0 and real output is fraction */
	*i = 0L;
	emovo( xi, frac );
	return;
	}
if( k > 31 )
	{
/*
;	long integer overflow: output large integer
;	and correct fraction
*/
	if( xi[0] )
		*i = 0x80000000;
	else
		*i = 0x7fffffff;
	(void )eshift( xi, k );
	goto lab11;
	}

if( k > 16 )
	{
/*
; shift more than 16 bits: shift up k-16, output the integer,
; then complete the shift to get the fraction.
*/
	k -= 16;
	(void )eshift( xi, k );

	*i = (long )(((unsigned long )xi[M] << 16) | xi[M+1]);
	eshup6( xi );
	goto lab10;
	}

/* shift not more than 16 bits */
(void )eshift( xi, k );
*i = (long )xi[M] & 0xffff;

lab10:

if( xi[0] )
	*i = -(*i);
lab11:

xi[0] = 0;
xi[E] = EXONE - 1;
xi[M] = 0;
if( (k = enormlz( xi )) > NBITS )
	ecleaz( xi );
else
	xi[E] -= (unsigned short )k;

emovo( xi, frac );
}


/*
;	Find unsigned long integer and fractional parts

;	unsigned long i;
;	unsigned short x[NE], frac[NE];
;	xifrac( x, &i, frac );

  A negative e type input yields integer output = 0
  but correct fraction.
*/
void euifrac( x, i, frac )
unsigned short *x;
long *i;
unsigned short *frac;
{
unsigned short xi[NI];
int k;

emovi( x, xi );
k = (int )xi[E] - (EXONE - 1);
if( k <= 0 )
	{
/* if exponent <= 0, integer = 0 and argument is fraction */
	*i = 0L;
	emovo( xi, frac );
	return;
	}
if( k > 32 )
	{
/*
;	long integer overflow: output large integer
;	and correct fraction
*/
	*i = 0xffffffff;
	(void )eshift( xi, k );
	goto lab10;
	}

if( k > 16 )
	{
/*
; shift more than 16 bits: shift up k-16, output the integer,
; then complete the shift to get the fraction.
*/
	k -= 16;
	(void )eshift( xi, k );

	*i = (long )(((unsigned long )xi[M] << 16) | xi[M+1]);
	eshup6( xi );
	goto lab10;
	}

/* shift not more than 16 bits */
(void )eshift( xi, k );
*i = (long )xi[M] & 0xffff;

lab10:

if( xi[0] )
	*i = 0L;

xi[0] = 0;
xi[E] = EXONE - 1;
xi[M] = 0;
if( (k = enormlz( xi )) > NBITS )
	ecleaz( xi );
else
	xi[E] -= (unsigned short )k;

emovo( xi, frac );
}



/*
;	Shift significand
;
;	Shifts significand area up or down by the number of bits
;	given by the variable sc.
*/
int eshift( x, sc )
unsigned short *x;
int sc;
{
unsigned short lost;
unsigned short *p;

if( sc == 0 )
	return( 0 );

lost = 0;
p = x + NI-1;

if( sc < 0 )
	{
	sc = -sc;
	while( sc >= 16 )
		{
		lost |= *p;	/* remember lost bits */
		eshdn6(x);
		sc -= 16;
		}

	while( sc >= 8 )
		{
		lost |= *p & 0xff;
		eshdn8(x);
		sc -= 8;
		}

	while( sc > 0 )
		{
		lost |= *p & 1;
		eshdn1(x);
		sc -= 1;
		}
	}
else
	{
	while( sc >= 16 )
		{
		eshup6(x);
		sc -= 16;
		}

	while( sc >= 8 )
		{
		eshup8(x);
		sc -= 8;
		}

	while( sc > 0 )
		{
		eshup1(x);
		sc -= 1;
		}
	}
if( lost )
	lost = 1;
return( (int )lost );
}



/*
;	normalize
;
; Shift normalizes the significand area pointed to by argument
; shift count (up = positive) is returned.
*/
int enormlz(x)
unsigned short x[];
{
register unsigned short *p;
int sc;

sc = 0;
p = &x[M];
if( *p != 0 )
	goto normdn;
++p;
if( *p & 0x8000 )
	return( 0 );	/* already normalized */
while( *p == 0 )
	{
	eshup6(x);
	sc += 16;
/* With guard word, there are NBITS+16 bits available.
 * return true if all are zero.
 */
	if( sc > NBITS )
		return( sc );
	}
/* see if high byte is zero */
while( (*p & 0xff00) == 0 )
	{
	eshup8(x);
	sc += 8;
	}
/* now shift 1 bit at a time */
while( (*p  & 0x8000) == 0)
	{
	eshup1(x);
	sc += 1;
	if( sc > NBITS )
		{
		mtherr( "enormlz", UNDERFLOW );
		return( sc );
		}
	}
return( sc );

/* Normalize by shifting down out of the high guard word
   of the significand */
normdn:

if( *p & 0xff00 )
	{
	eshdn8(x);
	sc -= 8;
	}
while( *p != 0 )
	{
	eshdn1(x);
	sc -= 1;

	if( sc < -NBITS )
		{
		mtherr( "enormlz", OVERFLOW );
		return( sc );
		}
	}
return( sc );
}




/* Convert e type number to decimal format ASCII string.
 * The constants are for 64 bit precision.
 */

#define NTEN 12
#define MAXP 4096

static unsigned short etens[NTEN+1][NE] = {
{0xc94c,0x979a,0x8a20,0x5202,0xc460,0x7525,},/* 10**4096 */
{0xa74d,0x5de4,0xc53d,0x3b5d,0x9e8b,0x5a92,},/* 10**2048 */
{0x650d,0x0c17,0x8175,0x7586,0xc976,0x4d48,},
{0xcc65,0x91c6,0xa60e,0xa0ae,0xe319,0x46a3,},
{0xddbc,0xde8d,0x9df9,0xebfb,0xaa7e,0x4351,},
{0xc66f,0x8cdf,0x80e9,0x47c9,0x93ba,0x41a8,},
{0x3cbf,0xa6d5,0xffcf,0x1f49,0xc278,0x40d3,},
{0xf020,0xb59d,0x2b70,0xada8,0x9dc5,0x4069,},
{0x0000,0x0000,0x0400,0xc9bf,0x8e1b,0x4034,},
{0x0000,0x0000,0x0000,0x2000,0xbebc,0x4019,},
{0x0000,0x0000,0x0000,0x0000,0x9c40,0x400c,},
{0x0000,0x0000,0x0000,0x0000,0xc800,0x4005,},
{0x0000,0x0000,0x0000,0x0000,0xa000,0x4002,}, /* 10**1 */
};

static unsigned short emtens[NTEN+1][NE] = {
{0x2de4,0x9fde,0xd2ce,0x04c8,0xa6dd,0x0ad8,}, /* 10**-4096 */
{0x4925,0x2de4,0x3436,0x534f,0xceae,0x256b,}, /* 10**-2048 */
{0x87a6,0xc0bd,0xda57,0x82a5,0xa2a6,0x32b5,},
{0x7133,0xd21c,0xdb23,0xee32,0x9049,0x395a,},
{0xfa91,0x1939,0x637a,0x4325,0xc031,0x3cac,},
{0xac7d,0xe4a0,0x64bc,0x467c,0xddd0,0x3e55,},
{0x3f24,0xe9a5,0xa539,0xea27,0xa87f,0x3f2a,},
{0x67de,0x94ba,0x4539,0x1ead,0xcfb1,0x3f94,},
{0x4c2f,0xe15b,0xc44d,0x94be,0xe695,0x3fc9,},
{0xfdc2,0xcefc,0x8461,0x7711,0xabcc,0x3fe4,},
{0xd3c3,0x652b,0xe219,0x1758,0xd1b7,0x3ff1,},
{0x3d71,0xd70a,0x70a3,0x0a3d,0xa3d7,0x3ff8,},
{0xcccd,0xcccc,0xcccc,0xcccc,0xcccc,0x3ffb,}, /* 10**-1 */
};

void e24toasc( x, string, ndigs )
unsigned short x[];
char *string;
int ndigs;
{
unsigned short w[NI];

#ifdef INFINITY
#ifdef IBMPC
if( (x[1] & 0x7f80) == 0x7f80 )
#else
if( (x[0] & 0x7f80) == 0x7f80 )
#endif
	{
#ifdef IBMPC
	if( x[1] & 0x8000 )
#else
	if( x[0] & 0x8000 )
#endif
		sprintf( string, " -Infinity " );
	else
		sprintf( string, " Infinity " );
	return;
	}
#endif
e24toe( x, w );
etoasc( w, string, ndigs );
}


void e53toasc( x, string, ndigs )
unsigned short x[];
char *string;
int ndigs;
{
unsigned short w[NI];

#ifdef INFINITY
#ifdef IBMPC
if( (x[3] & 0x7ff0) == 0x7ff0 )
#else
if( (x[0] & 0x7ff0) == 0x7ff0 )
#endif
	{
#ifdef IBMPC
	if( x[3] & 0x8000 )
#else
	if( x[0] & 0x8000 )
#endif
		sprintf( string, " -Infinity " );
	else
		sprintf( string, " Infinity " );
	return;
	}
#endif
e53toe( x, w );
etoasc( w, string, ndigs );
}


void e64toasc( x, string, ndigs )
unsigned short x[];
char *string;
int ndigs;
{
unsigned short w[NI];

#ifdef INFINITY
#ifdef IBMPC
if( (x[4] & 0x7fff) == 0x7fff )
#else
if( (x[0] & 0x7fff) == 0x7fff )
#endif
	{
#ifdef IBMPC
	if( x[4] & 0x8000 )
#else
	if( x[0] & 0x8000 )
#endif
		sprintf( string, " -Infinity " );
	else
		sprintf( string, " Infinity " );
	return;
	}
#endif
e64toe( x, w );
etoasc( w, string, ndigs );
}



void etoasc( x, string, ndigs )
unsigned short x[];
char *string;
int ndigs;
{
long digit;
unsigned short y[NI], t[NI], u[NI], w[NI];
unsigned short *p, *r, *ten;
unsigned short sign;
int i, j, k, expon, rndsav;
char *s, *ss;
unsigned short m;

rndsav = rndprc;
rndprc = NBITS;		/* set to full precision */
emov( x, y ); /* retain external format */
if( y[NE-1] & 0x8000 )
	{
	sign = 0xffff;
	y[NE-1] &= 0x7fff;
	}
else
	{
	sign = 0;
	}
expon = 0;
ten = &etens[NTEN][0];
emov( eone, t );
/* Test for zero exponent */
if( y[NE-1] == 0 )
	{
	for( k=0; k<NE-1; k++ )
		{
		if( y[k] != 0 )
			goto tnzro; /* denormalized number */
		}
	goto isone; /* legal all zeros */
	}
tnzro:

/* Test for infinity.  Don't bother with illegal infinities.
 */
if( y[NE-1] == 0x7fff )
	{
	if( sign )
		sprintf( string, " -Infinity " );
	else
		sprintf( string, " Infinity " );
	goto bxit;
	}

/* Test for exponent nonzero but significand denormalized.
 * This is an error condition.
 */
if( (y[NE-1] != 0) && ((y[NE-2] & 0x8000) == 0) )
	{
	mtherr( "etoasc", DOMAIN );
	sprintf( string, "NaN" );
	goto bxit;
	}

/* Compare to 1.0 */
i = ecmp( eone, y );
if( i == 0 )
	goto isone;

if( i < 0 )
	{ /* Number is greater than 1 */
/* Convert significand to an integer and strip trailing decimal zeros. */
	emov( y, u );
	u[NE-1] = EXONE + NBITS - 1;

	p = &etens[NTEN-4][0];
	m = 16;
do
	{
	ediv( p, u, t );
	efloor( t, w );
	for( j=0; j<NE-1; j++ )
		{
		if( t[j] != w[j] )
			goto noint;
		}
	emov( t, u );
	expon += (int )m;
noint:
	p += NE;
	m >>= 1;
	}
while( m != 0 );

/* Rescale from integer significand */
	u[NE-1] += y[NE-1] - (unsigned int )(EXONE + NBITS - 1);
	emov( u, y );
/* Find power of 10 */
	emov( eone, t );
	m = MAXP;
	p = &etens[0][0];
	while( ecmp( ten, u ) <= 0 )
		{
		if( ecmp( p, u ) <= 0 )
			{
			ediv( p, u, u );
			emul( p, t, t );
			expon += (int )m;
			}
		m >>= 1;
		if( m == 0 )
			break;
		p += NE;
		}
	}
else
	{ /* Number is less than 1.0 */
/* Pad significand with trailing decimal zeros. */
	if( y[NE-1] == 0 )
		{
		while( (y[NE-2] & 0x8000) == 0 )
			{
			emul( ten, y, y );
			expon -= 1;
			}
		}
	else
		{
		emovi( y, w );
		for( i=0; i<NDEC+1; i++ )
			{
			if( (w[NI-1] & 0x7) != 0 )
				break;
/* multiply by 10 */
			emovz( w, u );
			eshdn1( u );
			eshdn1( u );
			eaddm( w, u );
			u[1] += 3;
			while( u[2] != 0 )
				{
				eshdn1(u);
				u[1] += 1;
				}
			if( u[NI-1] != 0 )
				break;
			if( eone[NE-1] <= u[1] )
				break;
			emovz( u, w );
			expon -= 1;
			}
		emovo( w, y );
		}
	k = -MAXP;
	p = &emtens[0][0];
	r = &etens[0][0];
	emov( y, w );
	emov( eone, t );
	while( ecmp( eone, w ) > 0 )
		{
		if( ecmp( p, w ) >= 0 )
			{
			emul( r, w, w );
			emul( r, t, t );
			expon += k;
			}
		k /= 2;
		if( k == 0 )
			break;
		p += NE;
		r += NE;
		}
	ediv( t, eone, t );
	}
isone:
/* Find the first (leading) digit. */
emovi( t, w );
emovz( w, t );
emovi( y, w );
emovz( w, y );
eiremain( t, y );
digit = equot[NI-1];
while( (digit == 0) && (ecmp(y,ezero) != 0) )
	{
	eshup1( y );
	emovz( y, u );
	eshup1( u );
	eshup1( u );
	eaddm( u, y );
	eiremain( t, y );
	digit = equot[NI-1];
	expon -= 1;
	}
s = string;
if( sign )
	*s++ = '-';
else
	*s++ = ' ';
*s++ = (char )digit + '0';
*s++ = '.';
/* Examine number of digits requested by caller. */
if( ndigs < 0 )
	ndigs = 0;
if( ndigs > NDEC )
	ndigs = NDEC;
/* Generate digits after the decimal point. */
for( k=0; k<=ndigs; k++ )
	{
/* multiply current number by 10, without normalizing */
	eshup1( y );
	emovz( y, u );
	eshup1( u );
	eshup1( u );
	eaddm( u, y );
	eiremain( t, y );
	*s++ = (char )equot[NI-1] + '0';
	}
digit = equot[NI-1];
--s;
ss = s;
/* round off the ASCII string */
if( digit > 4 )
	{
/* Test for critical rounding case in ASCII output. */
	if( digit == 5 )
		{
		emovo( y, t );
		if( ecmp(t,ezero) != 0 )
			goto roun;	/* round to nearest */
		if( (*(s-1) & 1) == 0 )
			goto doexp;	/* round to even */
		}
/* Round up and propagate carry-outs */
roun:
	--s;
	k = *s & 0x7f;
/* Carry out to most significant digit? */
	if( k == '.' )
		{
		--s;
		k = *s;
		k += 1;
		*s = (char )k;
/* Most significant digit carries to 10? */
		if( k > '9' )
			{
			expon += 1;
			*s = '1';
			}
		goto doexp;
		}
/* Round up and carry out from less significant digits */
	k += 1;
	*s = (char )k;
	if( k > '9' )
		{
		*s = '0';
		goto roun;
		}
	}
doexp:
/*
if( expon >= 0 )
	sprintf( ss, "e+%d", expon );
else
	sprintf( ss, "e%d", expon );
*/
	sprintf( ss, "E%d", expon );
bxit:
rndprc = rndsav;
}




/*
;								ASCTOQ
;		ASCTOQ.MAC		LATEST REV: 11 JAN 84
;					SLM, 3 JAN 78
;
;	Convert ASCII string to quadruple precision floating point
;
;		Numeric input is free field decimal number
;		with max of 15 digits with or without 
;		decimal point entered as ASCII from teletype.
;	Entering E after the number followed by a second
;	number causes the second number to be interpreted
;	as a power of 10 to be multiplied by the first number
;	(i.e., "scientific" notation).
;
;	Usage:
;		asctoq( string, q );
*/

/* ASCII to single */
void asctoe24( s, y )
char *s;
unsigned short *y;
{
asctoeg( s, y, 24 );
}


/* ASCII to double */
void asctoe53( s, y )
char *s;
unsigned short *y;
{
#ifdef DEC
asctoeg( s, y, 56 );
#else
asctoeg( s, y, 53 );
#endif
}


/* ASCII to long double */
void asctoe64( s, y )
char *s;
unsigned short *y;
{
asctoeg( s, y, 64 );
}

/* ASCII to super double */
void asctoe( s, y )
char *s;
unsigned short *y;
{
asctoeg( s, y, NBITS );
}

/* Space to make a copy of the input string: */
static char lstr[82] = {0};

void asctoeg( ss, y, oprec )
char *ss;
unsigned short *y;
int oprec;
{
unsigned short yy[NI], xt[NI], tt[NI];
int esign, decflg, sgnflg, nexp, exp, prec, lost;
int k, trail, c, rndsav;
long lexp;
unsigned short nsign, *p;
char *sp, *s;

/* Copy the input string. */
s = ss;
while( *s == ' ' ) /* skip leading spaces */
	++s;
sp = lstr;
for( k=0; k<79; k++ )
	{
	if( (*sp++ = *s++) == '\0' )
		break;
	}
*sp = '\0';
s = lstr;

rndsav = rndprc;
rndprc = NBITS; /* Set to full precision */
lost = 0;
nsign = 0;
decflg = 0;
sgnflg = 0;
nexp = 0;
exp = 0;
prec = 0;
ecleaz( yy );
trail = 0;

nxtcom:
k = *s - '0';
if( (k >= 0) && (k <= 9) )
	{
/* Ignore leading zeros */
	if( (prec == 0) && (decflg == 0) && (k == 0) )
		goto donchr;
/* Identify and strip trailing zeros after the decimal point. */
	if( (trail == 0) && (decflg != 0) )
		{
		sp = s;
		while( (*sp >= '0') && (*sp <= '9') )
			++sp;
/* Check for syntax error */
		c = *sp & 0x7f;
		if( (c != 'e') && (c != 'E') && (c != '\0')
			&& (c != '\n') && (c != '\r') && (c != ' ')
			&& (c != ',') )
			goto error;
		--sp;
		while( *sp == '0' )
			*sp-- = 'z';
		trail = 1;
		if( *s == 'z' )
			goto donchr;
		}
/* If enough digits were given to more than fill up the yy register,
 * continuing until overflow into the high guard word yy[2]
 * guarantees that there will be a roundoff bit at the top
 * of the low guard word after normalization.
 */
	if( yy[2] == 0 )
		{
		if( decflg )
			nexp += 1; /* count digits after decimal point */
		eshup1( yy );	/* multiply current number by 10 */
		emovz( yy, xt );
		eshup1( xt );
		eshup1( xt );
		eaddm( xt, yy );
		ecleaz( xt );
		xt[NI-2] = (unsigned short )k;
		eaddm( xt, yy );
		}
	else
		{
		lost |= k;
		}
	prec += 1;
	goto donchr;
	}

switch( *s )
	{
	case 'z':
		break;
	case 'E':
	case 'e':
		goto expnt;
	case '.':	/* decimal point */
		if( decflg )
			goto error;
		++decflg;
		break;
	case '-':
		nsign = 0xffff;
		if( sgnflg )
			goto error;
		++sgnflg;
		break;
	case '+':
		if( sgnflg )
			goto error;
		++sgnflg;
		break;
	case ',':
	case ' ':
	case '\0':
	case '\n':
	case '\r':
		goto daldone;
	case 'i':
	case 'I':
		ecleaz(yy);
		yy[E] = 0x7fff;  /* infinity */
		goto aexit;
	default:
	error:
		mtherr( "asctoe", DOMAIN );
		eclear(y);
		goto aexit;
	}
donchr:
++s;
goto nxtcom;

/* Exponent interpretation */
expnt:

esign = 1;
exp = 0;
++s;
/* check for + or - */
if( *s == '-' )
	{
	esign = -1;
	++s;
	}
if( *s == '+' )
	++s;
while( (*s >= '0') && (*s <= '9') )
	{
	exp *= 10;
	exp += *s++ - '0';
	}
if( esign < 0 )
	exp = -exp;

daldone:
nexp = exp - nexp;
/* Pad trailing zeros to minimize power of 10, per IEEE spec. */
while( (nexp > 0) && (yy[2] == 0) )
	{
	emovz( yy, xt );
	eshup1( xt );
	eshup1( xt );
	eaddm( yy, xt );
	eshup1( xt );
	if( xt[2] != 0 )
		break;
	nexp -= 1;
	emovz( xt, yy );
	}
if( (k = enormlz(yy)) > NBITS )
	{
	ecleaz(yy);
	goto aexit;
	}
lexp = (EXONE - 1 + NBITS) - k;
emdnorm( yy, lost, 0, lexp, 64 );
/* convert to external format */


/* Multiply by 10**nexp.  If precision is 64 bits,
 * the maximum relative error incurred in forming 10**n
 * for 0 <= n <= 324 is 8.2e-20, at 10**180.
 * For 0 <= n <= 999, the peak relative error is 1.4e-19 at 10**947.
 * For 0 >= n >= -999, it is -1.55e-19 at 10**-435.
 */
lexp = yy[E];
if( nexp == 0 )
	{
	k = 0;
	goto expdon;
	}
esign = 1;
if( nexp < 0 )
	{
	nexp = -nexp;
	esign = -1;
	if( nexp > 4096 )
		{ /* Punt.  Can't handle this without 2 divides. */
		emovi( etens[0], tt );
		lexp -= tt[E];
		k = edivm( tt, yy );
		lexp += EXONE;
		nexp -= 4096;
		}
	}
p = &etens[NTEN][0];
emov( eone, xt );
exp = 1;
do
	{
	if( exp & nexp )
		emul( p, xt, xt );
	p -= NE;
	exp = exp + exp;
	}
while( exp <= MAXP );

emovi( xt, tt );
if( esign < 0 )
	{
	lexp -= tt[E];
	k = edivm( tt, yy );
	lexp += EXONE;
	}
else
	{
	lexp += tt[E];
	k = emulm( tt, yy );
	lexp -= EXONE - 1;
	}

expdon:

/* Round and convert directly to the destination type */
if( oprec == 53 )
	lexp -= EXONE - 0x3ff;
else if( oprec == 24 )
	lexp -= EXONE - 0177;
#ifdef DEC
else if( oprec == 56 )
	lexp -= EXONE - 0201;
#endif
rndprc = oprec;
emdnorm( yy, k, 0, lexp, 64 );

aexit:

rndprc = rndsav;
yy[0] = nsign;
switch( oprec )
	{
#ifdef DEC
	case 56:
		todec( yy, y ); /* see etodec.c */
		break;
#endif
	case 53:
		toe53( yy, y );
		break;
	case 24:
		toe24( yy, y );
		break;
	case 64:
		toe64( yy, y );
		break;
	case NBITS:
		emovo( yy, y );
		break;
	}
}


 
/* y = largest integer not greater than x
 * (truncated toward minus infinity)
 *
 * unsigned short x[NE], y[NE]
 *
 * efloor( x, y );
 */
static unsigned short bmask[] = {
0xffff,
0xfffe,
0xfffc,
0xfff8,
0xfff0,
0xffe0,
0xffc0,
0xff80,
0xff00,
0xfe00,
0xfc00,
0xf800,
0xf000,
0xe000,
0xc000,
0x8000,
0x0000,
};

void efloor( x, y )
unsigned short x[], y[];
{
register unsigned short *p;
int e, expon, i;
unsigned short f[NE];

emov( x, f ); /* leave in external format */
expon = (int )f[NE-1];
e = (expon & 0x7fff) - (EXONE - 1);
if( e <= 0 )
	{
	eclear(y);
	goto isitneg;
	}
/* number of bits to clear out */
e = NBITS - e;
emov( f, y );
if( e <= 0 )
	return;

p = &y[0];
while( e >= 16 )
	{
	*p++ = 0;
	e -= 16;
	}
/* clear the remaining bits */
*p &= bmask[e];
/* truncate negatives toward minus infinity */
isitneg:

if( (unsigned short )expon & (unsigned short )0x8000 )
	{
	for( i=0; i<NE-1; i++ )
		{
		if( f[i] != y[i] )
			{
			esub( eone, y, y );
			break;
			}
		}
	}
}


/* unsigned short x[], s[];
 * long *exp;
 *
 * efrexp( x, exp, s );
 *
 * Returns s and exp such that  s * 2**exp = x and .5 <= s < 1.
 * For example, 1.1 = 0.55 * 2**1
 * Handles denormalized numbers properly using long integer exp.
 */
void efrexp( x, exp, s )
unsigned short x[];
long *exp;
unsigned short s[];
{
unsigned short xi[NI];
long li;

emovi( x, xi );
li = (long )((short )xi[1]);

if( li == 0 )
	{
	li -= enormlz( xi );
	}
xi[1] = 0x3ffe;
emovo( xi, s );
*exp = li - 0x3ffe;
}



/* unsigned short x[], y[];
 * long pwr2;
 *
 * eldexp( x, pwr2, y );
 *
 * Returns y = x * 2**pwr2.
 */
void eldexp( x, pwr2, y )
unsigned short x[];
long pwr2;
unsigned short y[];
{
unsigned short xi[NI];
long li;
int i;

emovi( x, xi );
li = xi[1];
li += pwr2;
i = 0;
emdnorm( xi, i, i, li, 64 );
emovo( xi, y );
}


/* c = remainder after dividing b by a
 * Least significant integer quotient bits left in equot[].
 */
void eremain( a, b, c )
unsigned short a[], b[], c[];
{
unsigned short den[NI], num[NI];

if( ecmp(a,ezero) == 0 )
	{
	mtherr( "eremain", SING );
	eclear( c );
	return;
	}
emovi( a, den );
emovi( b, num );
eiremain( den, num );
/* Sign of remainder = sign of quotient */
if( a[0] == b[0] )
	num[0] = 0;
else
	num[0] = 0xffff;
emovo( num, c );
}

void eiremain( den, num )
unsigned short den[], num[];
{
long ld, ln;
unsigned short j;

ld = den[E];
ld -= enormlz( den );
ln = num[E];
ln -= enormlz( num );
ecleaz( equot );
while( ln >= ld )
	{
	if( ecmpm(den,num) <= 0 )
		{
		esubm(den, num);
		j = 1;
		}
	else
		{
		j = 0;
		}
	eshup1(equot);
	equot[NI-1] |= j;
	eshup1(num);
	ln -= 1;
	}
emdnorm( num, 0, 0, ln, 0 );
}

