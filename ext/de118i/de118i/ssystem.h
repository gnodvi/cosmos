//------------------------------------------------------------------------------
/* Configuration settings for ssystem.c
 */
//------------------------------------------------------------------------------

#define DE200 0
#define DE102 0
#define DE118 1

//------------------------------------------------------------------------------
//
// дальше начинаются извращения в попытке описать расположение данных для разных
// объектов в одном массиве !! (((
// 
//------------------------------------------------------------------------------


/* Define nonzero to add an extra body such as an asteroid.
   You must manually insert the desired initial state vectors
   just before the moon in aconst.h.
   Otherwise define this as 0. */
#define EXTRA_BODY 0

/* Compute relativity corrections, or not: */
#define DOREL 1

/* Include the Moon as a separate body, or not: */
#define MOON   1
#define LIBRAT 1

/* Include the five asteroids that are in JPL's model, or not: */
#define AROIDS 1

/* Index to first mass. The lunar libration angles
 * are the first items in the array.
 */
// т.е. это как фиктивное тело (масса)? но там 6 коэфф-тов илбрации 

#if LIBRAT
#define FMASS 1  
#else
#define FMASS 0
#endif

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if MOON

#define NBODY (11+EXTRA_BODY)
#define IMOON (9+FMASS+EXTRA_BODY)
#else

#undef  LIBRAT
#define LIBRAT 0

#undef  FMASS
#define FMASS 0

#define NBODY (10+EXTRA_BODY)

#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define IAROIDS (NBODY+FMASS)

#if AROIDS
#define NTOTAL (NBODY+FMASS+5)
#endif

#ifndef NTOTAL
#define NTOTAL (NBODY+FMASS)
#endif

#define IEARTH (2+FMASS)
#define ISUN (NBODY+FMASS-1)

#define NEQ (6*NTOTAL)

#define NEQC (6*(NBODY+FMASS))

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------

