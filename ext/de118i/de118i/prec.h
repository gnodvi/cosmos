
#ifndef PREC_H
#define PREC_H
/* Definitions for long double precision */


#ifndef LDOUBLE
#define LDOUBLE 1
#endif


// т.е. здесь всегда LDOUBLE !!!!!!!!!!
// 

#if LDOUBLE

#define DOUBLE long double
#define SQRT sqrtl
#define SIN sinl
#define COS cosl
#define TAN tanl
#define LOG logl
#define ASIN asinl
#define ATAN atanl
#define FLOOR floorl
#define FABS fabsl

#else

#define DOUBLE double
#define SQRT sqrt
#define SIN sin
#define COS cos
#define TAN tan
#define LOG log
#define ASIN asin
#define ATAN atan
#define FLOOR floor
#define FABS fabs

#endif

DOUBLE SQRT(), SIN(), COS(), TAN(), LOG(), ASIN(), ATAN();
DOUBLE FLOOR(), FABS();

#endif
