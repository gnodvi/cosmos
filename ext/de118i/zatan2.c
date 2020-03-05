/*							atan2()
 *
 *	Quadrant correct inverse circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, z, atan2();
 *
 * z = atan2( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between 0 and +2pi whose tangent
 * is y/x.
 *
 *
 *
 * ACCURACY:
 *
 * See atan.c.
 *
 */


/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
Certain routines from the Library, including this one, may
be used and distributed freely provided this notice is retained
and source code is included with all distributions.
*/

#include "prec.h"
extern long double atanl (long double);

extern DOUBLE PI, Two, OneandaHalf, Half, Zero;

DOUBLE zatan2( x, y )
DOUBLE x, y;
{
DOUBLE z, w;
short code;
DOUBLE ATAN();


code = 0;

if( x < 0 )
	code = 2;
if( y < 0 )
	code |= 1;

if( x == 0 )
	{
	if( code & 1 )
		return( OneandaHalf*PI );
	if( y == Zero )
		return( Zero );
	return( Half*PI );
	}

if( y == 0 )
	{
	if( code & 2 )
		return( PI );
	return( Zero );
	}


if( code == 0 )
	w = Zero;
else if( code == 1 )
	w = Two*PI;
else
	w = PI;
/*
switch( code )
	{
	case 0: w = 0.0; break;
	case 1: w = 2.0 * PI; break;
	case 2:
	case 3: w = PI; break;
	}
*/
z = ATAN( y/x );

return( w + z );
}
