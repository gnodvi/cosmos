// -*-   mode: c    coding: koi8   -*- -----------------------------------------

//------------------------------------------------------------------------------
/* Main program for numerical integration of the solar system.
 * See the readme file for discussion.
 *
 * Before compiling, set type of computer in mconf.h,
 * arithmetic precision in prec.h, and model parameters
 * in ssystem.h.
 *
 * If DE200 compatible outputs are desired, include the
 * rotation r118_200() in the output routine.
 */
//------------------------------------------------------------------------------


/* install trap handler for arithmetic debugging */
#define FPESHOW 0

#ifndef NOINCS
#include "mconf.h"
#include "prec.h"
#include "ssystem.h"
#endif

#include "int.h"
#include "const.h"

#ifdef DEBUG
#undef DEBUG
#endif
#define DEBUG 0

extern DOUBLE GMs[], KG, C, EMRAT;
extern DOUBLE yn0[], JD0, yn1[], JD1;

/* Define 1 for ephemeris output in double precision even if LDOUBLE=1.
 * The saved integrator state will still be LDOUBLE.
 */
#define OUTDOUBLE 1

DOUBLE vnewt[6*NTOTAL];

DOUBLE *awork; /* Adams integrator work array */
DOUBLE *rwork; /* Runge-Kutta work array */

DOUBLE Rij[NTOTAL][NTOTAL];
DOUBLE Rij3[NTOTAL][NTOTAL];

int fd;
char fname[80];

/* #if MSC */
/* #include <stdio.h> */
/* #include <sys\types.h> */
/* #include <sys\stat.h> */
/* #include <fcntl.h> */
/* #include <io.h> */
/* #include <stdlib.h> */
/* #include <malloc.h> */
/* #endif */

//#if UNIX
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
  //#endif

/* #if VAX */
/* #include <stdio.h> */
/* #include <file.h> */
/* #include <stat.h> */
/* #endif */

/* #if RSX */
/* #include <stdio.h> */
/* #include <fcntl.h> */
/* #define CPERM 1 */
/* #endif */

/* #if SVC */
/* #endif */


//------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <netinet/in.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>
#include <libgen.h>

//------------------------------------------------------------------------------


DOUBLE ssyt, err, ssyh, unused;
int ord, outsteps;
long isamps, nsamps, nsteps;

//------------------------------------------------------------------------------
/* Subroutines to save and restore the integrator state
 */
//------------------------------------------------------------------------------
void savstate()
{

  char *save_fname = "de118.sav";

  /*int f, i;*/
  register int f;

  f = open( save_fname, O_CREAT | O_TRUNC | O_RDWR, S_IWRITE | S_IREAD );

  sinternal(f);

  write( f, &ssyt,     sizeof(DOUBLE) );
  write( f, &err,      sizeof(DOUBLE) );
  write( f, &ssyh,     sizeof(DOUBLE) );
  write( f, &unused,   sizeof(DOUBLE) );
  write( f, yn0,       6*NTOTAL*sizeof(DOUBLE) );
  write( f, &ord,      sizeof(int) );
  write( f, &outsteps, sizeof(int) );
  write( f, awork,     ((NEQ+1)*(MAXORD+2)+(6*NEQ))*sizeof(DOUBLE) );

  /*write( f, awork, i );*/
  close( f);

}
//------------------------------------------------------------------------------
void resstate ( )
{

  char *save_fname = "de118.sav";

  /*int f, i;*/
  register int f;

  f = open( save_fname, O_RDONLY, S_IREAD | S_IWRITE );

  if( f <= 0 )
  {
    printf( "Can't find $s \n", save_fname);
    exit (0);
  }

  rinternal (f, awork);

  read( f, &ssyt,     sizeof(DOUBLE) );
  read( f, &err,      sizeof(DOUBLE) );
  read( f, &ssyh,     sizeof(DOUBLE) );
  read( f, &unused,   sizeof(DOUBLE) );
  read( f, yn0,       6*NTOTAL*sizeof(DOUBLE) );
  read( f, &ord,      sizeof(int) );
  read( f, &outsteps, sizeof(int) );
  read( f, awork,     ((NEQ+1)*(MAXORD+2)+(6*NEQ))*sizeof(DOUBLE) );
  /*read( f, awork, i );*/

  close( f);

}
//------------------------------------------------------------------------------

/* Double precision variables for printf() */
static double dx, dy, dz;
/* DOUBLE adstep();*/

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------

/* Subroutine func() calculates the forces and accelerations.
  * Data in the output vector v[] are in the order
  * d^2x/dt^2, dx/dt, d^2y/dt^2, dy/dt, d^2z/dt^2, dz/dt
  * for each object.  For this problem the velocities dx/dt, ...
  * are simply copied over from the input array yw[].
  */

#if MOON

DOUBLE yw[6*NTOTAL];
DOUBLE ymoon[6];

#endif

// func - вызывается из "adams4/adstep" и "runge.c/rungek" на этапе шагов
// 
// 

//------------------------------------------------------------------------------
func ( t, yin, v )

  DOUBLE t;     /* time */
  DOUBLE yin[]; /* input state: velocity and position */
  DOUBLE v[];   /* output: calculated acceleration, copy of input velocity */
{

  DOUBLE xs, ys, zs, xd, yd, zd, xv, yv, zv, temp;
  int i, j, ii, jj;

#if DEBUG
  printf( "func " );
#endif

#if MOON
  /* Locally replace input variable EMB by barycentric earth
    * and input variable Moon-Earth by barycentric Moon.
    */
  fromemb( yin, yw );

#if DEBUG
  printf( "fromemb " );
#endif
#endif

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /* asteroid positions */
#if AROIDS

  aroids ( t, yw );

#if DEBUG
  printf ("aroids " );
#endif
#endif
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  /* Table of distances between objects i and j */
  distances (yw);

#if DEBUG
  printf( "distances " );
#endif

  fixsun( yw, v );
#if DEBUG
  printf( "fixsun " );
#endif

#if MOON

  //  это для SUN ??
  // 
  for (j = 0; j < 6; j++ ) {

    yin[(6*ISUN) + j] = yw[(6*ISUN) + j];
  }

#endif

  /* Compute Newtonian acceleration. */
  //
  ii = 6 * FMASS; // 6*1

  for( i=FMASS/*1*/; i<NTOTAL; i++ )
  {
    xs = Zero;
    ys = Zero;
    zs = Zero;
    xv = yw[ii+1];
    yv = yw[ii+3];
    zv = yw[ii+5];

    jj = 6*FMASS;

    for( j=FMASS; j<NTOTAL; j++ )
    {
      if( (j != i) && (j != ISUN) )
      {
        xd = yw[jj+1] - xv;
        yd = yw[jj+3] - yv;
        zd = yw[jj+5] - zv;
        temp = GMs[j] * Rij3[i][j];
        xs += xd * temp;
        ys += yd * temp;
        zs += zd * temp;
      }
      jj += 6;
    }

    /* Add the Sun last. */
    // 
    if( i != ISUN )
    {
      xd = yw[6*ISUN+1] - xv;
      yd = yw[6*ISUN+3] - yv;
      zd = yw[6*ISUN+5] - zv;
      temp = GMs[ISUN] * Rij3[i][ISUN];
      xs += xd * temp;
      ys += yd * temp;
      zs += zd * temp;
    }

    vnewt[ii]   = xs;        /* acceleration */
    vnewt[ii+2] = ys;
    vnewt[ii+4] = zs;
    vnewt[ii+1] = yw[ii];    /* velocity */
    vnewt[ii+3] = yw[ii+2];
    vnewt[ii+5] = yw[ii+4];
    ii += 6;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if DOREL

  reltiv (yw, v);

#if DEBUG
  printf( "reltiv " );
#endif

#else

  /* No relativity theory. Return the Newtonian accelerations. */
  //
  ii = 6 * FMASS;

  for( i = FMASS; i < IAROIDS; i++ )
  {
    v[ii]   = vnewt[ii];
    v[ii+2] = vnewt[ii+2];
    v[ii+4] = vnewt[ii+4];
    v[ii+1] = yw[ii];
    v[ii+3] = yw[ii+2];
    v[ii+5] = yw[ii+4];

    ii += 6;
  }
#endif /* DOREL */
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if MOON

  oblate ( t, yw, v, ymoon );

#if DEBUG
  printf( "oblate\n" );
#endif

  /* Add the Newtonian accelerations last. */

  ii = 6*FMASS;
  for( i=FMASS; i<IAROIDS; i++ )
  {
    v[ii] += vnewt[ii];
    v[ii+2] += vnewt[ii+2];
    v[ii+4] += vnewt[ii+4];
    /*
      *	v[ii+1] = yw[ii];
      *	v[ii+3] = yw[ii+2];
      *	v[ii+5] = yw[ii+4];
      */
    ii += 6;
  }

  /* Convert barycentric Earth and Moon to output EMB and M variables. */

  ii = 6 * IEARTH;
  jj = 6 * IMOON;

  for (i = 0; i < 6; i += 2 )
  {
    xd = v[ii+i]; /* Earth */
    yd = v[jj+i]; /* Moon  */
    v[ii+i] = (EMRAT * xd + yd)/(EMRAT+One); /* EMB */
    v[jj+i] = yd - xd; /* M = Moon - Earth */
    v[ii+i+1] = yin[ii+i]; /* copy original velocity */
    v[jj+i+1] = yin[jj+i];
  }

#if 0

  /* Display initial acceleration discrepancies (see findcent.c) */
  chkacc (v);

  exit (0);

#endif

#endif /* MOON */
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}
//------------------------------------------------------------------------------

/* Constrain the barycenter to stay at the origin.
  * This is done by adjusting the Sun.
  */

DOUBLE mustar[NTOTAL];

//------------------------------------------------------------------------------
fixsun ( y, v )

  DOUBLE y[], v[];
{
  DOUBLE xx, yy, zz, vx, vy, vz;
  DOUBLE csqi, s;
  int i, j, k, ii, jj;


  csqi = Half / (C*C);
  for( k=0; k<2; k++ )
  { /* Iterate to find solution of implicit expressions. */


    /* Relativistic GM of each body */

    ii = 6*FMASS;
    for( i=FMASS; i<NTOTAL; i++ )
    {
      vx = y[ii]; /* velocity */
      s = vx * vx;
      vx = y[ii+2];
      s += vx * vx;
      vx = y[ii+4];
      s += vx * vx;  /* square of velocity */
      for( j=FMASS; j<NTOTAL; j++ )
      {
        if( j == i )
          continue;
        s -= GMs[j]*Rij[i][j];
      }
      /*
        *	mustar[i] = GMs[i] * (One + csqi * s);
        */
      mustar[i] = GMs[i] + GMs[i] * csqi * s;
      ii += 6;
    }

    /* Compute center of mass with Sun omitted. */

    xx = Zero;
    yy = Zero;
    zz = Zero;
    vx = Zero;
    vy = Zero;
    vz = Zero;

    ii = 6*FMASS;

    for( i=FMASS; i<NTOTAL; i++ )
    {
      if( i != ISUN )
      {
        s = mustar[i];
        xx += y[ii+1] * s; /* position */
        yy += y[ii+3] * s;
        zz += y[ii+5] * s;
        vx += y[ii] * s; /* velocity */
        vy += y[ii+2] * s;
        vz += y[ii+4] * s;
      }
      ii += 6;
    }

    /* Evaluate the Sun so that the center = 0. */

    s = mustar[ISUN];
    xx = -xx/s;
    yy = -yy/s;
    zz = -zz/s;
    vx = -vx/s;
    vy = -vy/s;
    vz = -vz/s;
    jj = 6*ISUN;
    y[jj+1] = xx;
    y[jj+3] = yy;
    y[jj+5] = zz;
    y[jj] = vx;
    y[jj+2] = vy;
    y[jj+4] = vz;
    v[jj+1] = vx;
    v[jj+3] = vy;
    v[jj+5] = vz;

    /* Adjust the table of distances between bodies i and j.
      * Note, most of this work was already done by func().
      * Only the entries involving the Sun need to be changed.
      */

    ii = 6*FMASS;
    for( j=FMASS; j<NTOTAL; j++ )
    {
      if( j != ISUN )
      {
        vx = xx - y[ii+1];
        s = vx * vx;
        vx = yy - y[ii+3];
        s += vx * vx;
        vx = zz - y[ii+5];
        s += vx * vx;
        vx = SQRT(s);
        vy = One/vx;
        Rij[ISUN][j] = vy;
        Rij[j][ISUN] = vy;
        s = One/(vx*s);
        Rij3[ISUN][j] = s;
        Rij3[j][ISUN] = s;
      }
      ii += 6;
    }

  } /* iteration */

}
//------------------------------------------------------------------------------
#if MOON

/* Convert EMB, M to barycentric Earth and Moon vectors.
  * ymoon[] declared externally.
  */

//------------------------------------------------------------------------------
fromemb ( yinp, y )

  DOUBLE yinp[], y[];
{
  DOUBLE zd, yd;
  int i, ii, jj;

  for( i=0; i<(6*NTOTAL); i++ )
    y[i] = yinp[i];

  for( i=0; i<6; i++ )
    ymoon[i] = yinp[6*IMOON+i]; /* M */

  ii = 6*IEARTH;
  jj = 6*IMOON;

  for( i=0; i<6; i++ )
  {
    zd = yinp[ii+i];		/* EMB */
    yd = yinp[jj+i]/(One+EMRAT);	/* M */
    y[ii+i] = zd -  yd;		/* r_e */
    y[jj+i] = zd + EMRAT * yd;	/* r_m */
  }

}
//------------------------------------------------------------------------------

#endif
//------------------------------------------------------------------------------

/* Make table of distances between ith and jth objects
  */
//------------------------------------------------------------------------------
distances (y)

  DOUBLE y[];
{
  register DOUBLE t, u;
  DOUBLE r, xv, yv, zv;
  int j, i, jj, ii;

  jj = 6*(FMASS+1);

  for(j=(FMASS+1); j<NTOTAL; j++)
  {
    ii = 6*FMASS;
    xv = y[jj+1];
    yv = y[jj+3];
    zv = y[jj+5];
    for(i=FMASS; i<j; i++)
    {
      t = xv - y[ii+1];
      u = t * t;
      t = yv - y[ii+3];
      u += t * t;
      t = zv - y[ii+5];
      u += t * t;
      r = SQRT(u);
      t = One/r;
      Rij[j][i] = t;
      Rij[i][j]  = t;
      t = One/(r*u);
      Rij3[j][i] = t;
      Rij3[i][j] = t;
      ii += 6;
    }
    jj += 6;
  }

#if MOON

  /* Take the input M vector for distance from Earth to Moon.
    * ymoon[] set by previous call to fromemb().
    */

  t = ymoon[1];
  u = t * t;
  t = ymoon[3];
  u += t * t;
  t = ymoon[5];
  u += t * t;
  r = SQRT(u);
  t = One/r;
  Rij[IEARTH][IMOON] = t;
  Rij[IMOON][IEARTH] = t;
  t = One/(r*u);
  Rij3[IEARTH][IMOON] = t;
  Rij3[IMOON][IEARTH] = t;
#endif

}
//------------------------------------------------------------------------------

/* Write date and solution vector to output file.
  */

#if OUTDOUBLE
double yout[6*(ISUN+1)]; /* Output array for disc file. */
#else
DOUBLE yout[6*(ISUN+1)];
#endif

//------------------------------------------------------------------------------
//wfile (t, y)
wfile (int fd, DOUBLE t, DOUBLE y[])

  //int fd;
  //DOUBLE t;
  //DOUBLE y[];
{

  DOUBLE p[3], v[3];
#if OUTDOUBLE
  double dt;
#else
  DOUBLE dt;
#endif
  int i, j, k;

  dt = t;
  /*for( i=0; i<6*(ISUN+1); i++ ) */
  for( i=0; i<6*FMASS; i++ )
    yout[i] = y[i];

  /* Rotate to DE200 coordinate system
    */
  for( i=FMASS; i<=ISUN; i++ )
  {
    k = 6 * i;
    for( j=0; j<3; j++ )
    {
      v[j] = y[k+j+j];
      p[j] = y[k+j+j+1];
    }
#if DE118
    r118_200(v);
    r118_200(p);
#endif
    for( j=0; j<3; j++ )
    {
      yout[k+j+j] = v[j];
      yout[k+j+j+1] = p[j];
    }
  }

#if OUTDOUBLE

  write( fd, &dt, sizeof(double) );
  write( fd, yout, 6*(ISUN+1)*sizeof(double) );
#else
  write( fd, &dt, sizeof(DOUBLE) );
  write( fd, yout, 6*(ISUN+1)*sizeof(DOUBLE) );

#endif

}
//------------------------------------------------------------------------------
// 
//
//
//
//------------------------------------------------------------------------------
void
get_options_print_help () 
{

  printf ("\n");
  printf ("USE: \n");
  printf ("./s_dial -t* -p* -o* \n");
  printf ("\n");      
  
}
//------------------------------------------------------------------------------
void
options_warning_NULL (char c) {

  fprintf (stderr, "\n");
  fprintf (stderr, "WARNING OPTIONS: pointer to -%c == NUL !! \n", c);
  //fprintf (stderr, "\n");

  return;
}
//------------------------------------------------------------------------------
void
get_options (int argc, char **argv, 
             //char     *p_ip_buf, 
             //uint16_t *port, 
             double   *interval_ms,
             int *count_max, 
             int *count_out, int *is_verbose, int *test_num,
             int *user_flag)
{
  char c;

  // optarg - указатель на текущий аргумент, если таковой имеется
  // optind - индекс на след. указатель argv (будет обработан при след. вызове)

  optind = 1; // чтобы бы повторном (вложенном вызове) начать с начала 

  //fprintf (stderr, "\n");
  //fprintf (stderr, "get_options: BEGIN  optind = %d   \n", optind);

  //if (argc == 1) {
  //  get_options_print_help ();
  //  return;
  //}
    

  // перебираем все параметры:
  // (в QNX эти параметры должны идти первыми)

  while ((c = getopt (argc, argv, "uha:p:i:n:o:vt:")) != -1) { 

    //fprintf (stderr, "get_options: optind = %d   c = %c \n", optind, c);

    switch (c) {

      //case 'a': // IP-адрес
      //  if (p_ip_buf) strcpy (p_ip_buf, optarg); 
      //  else              options_warning_NULL (c);
      //  break;

      //case 'p': 
      //  if (port) *port = atoi (optarg); 
      //  else              options_warning_NULL (c);
      //  break;

    case 'i': 
      if (interval_ms) *interval_ms = atoi (optarg); 
      else              options_warning_NULL (c);
      break;

    case 'n': 
      if (count_max) *count_max = atoi (optarg); 
      else            options_warning_NULL (c);
      break;

    case 'o': 
      if (count_out) *count_out = atoi (optarg); 
      else            options_warning_NULL (c);
      break;

    case 'v':     // verbose - многословный
      if (is_verbose) *is_verbose  = 1; // подробная печать 
      else             options_warning_NULL (c);
      break;

    case 't': 
      if (test_num) *test_num = atoi (optarg); 
       else          options_warning_NULL (c);
     break;

    case 'u':    
      if (user_flag) *user_flag  = 1;
      else           options_warning_NULL (c);
      break;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'h': 
      get_options_print_help (); 
      break;


    default:      
      exit( EXIT_FAILURE ); // сюда невозможно попасть!
    }
  }

  //fprintf (stderr, "get_options: END \n");

  return;
}
//------------------------------------------------------------------------------
// 
// 
// 
//------------------------------------------------------------------------------
void 
print_consts () {


  printf( "NEQ %d, NEQC %d, 6*NTOTAL %d \n",  NEQ, NEQC, 6*NTOTAL);

  printf ("\n");

  printf ("DOREL      = %d  Compute relativity corrections, or not \n",  DOREL);
  printf ("MOON       = %d  Include the Moon as a separate body, or not \n",  MOON);
  printf ("LIBRAT     = %d  \n",  LIBRAT);
  printf ("\n");

  printf ("EXTRA_BODY = %d  To add an extra body such as an asteroid, or not \n",  
          EXTRA_BODY);
  printf ("AROIDS     = %d  Include the five asteroids, or not \n",  
          AROIDS);
  printf ("\n");
  
  printf ("IEARTH  = %2d  \n",  IEARTH);
  printf ("IMOON   = %2d  \n",  IMOON);
  printf ("ISUN    = %2d  \n",  ISUN);
  printf ("IAROIDS = %2d  \n",  IAROIDS);
  printf ("\n");

  printf ("NBODY   = %2d  \n",  NBODY);
  printf ("FMASS   = %d  \n",  FMASS);
  printf ("NTOTAL  = %2d  \n",  NTOTAL);
  
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
void 
malloc_alls () {

  DOUBLE xa; 
  int i, ii;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ii = (NEQ+1) * (MAXORD+2) + (7*NEQ);

  awork = (DOUBLE *) malloc (ii * sizeof(DOUBLE));

  if (awork == 0)
  {
    printf ( "malloc(awork) failed\n" );
    exit (0);
  }

  /* Note, it is not really necessary to clear out the work arrays.
    */

  DOUBLE *pdv;

  pdv = awork;
  xa  = Zero;

  for (i = 0; i < ii; i++) {
    *pdv++ = xa;
  }

  printf ("Allocated %d bytes to awork \n", i);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  i = (MAXRUNG + 3) * NEQ * sizeof (DOUBLE);

  rwork = (DOUBLE *) malloc (i);
  if (rwork == 0)
  {
    printf ( "malloc(rwork) failed \n" );
    exit (0);
  }

  printf( "Allocated %d bytes to rwork\n", i );

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  fprintf (stderr, "\n");

  return;
}
//------------------------------------------------------------------------------
void
print_final_positions_and_errors (int iout) {

  DOUBLE xa, ya, za, xb, yb, zb; /* for display */
  int i, ii, restored;


#if LDOUBLE
static union
{
  long double ld;
  unsigned short s[8];
}  funion;
#else
static union
{
  double d;
  unsigned short s[4];
}  funion;
#endif



  printf ("\n");
  printf ("Final x, y, z, positions and errors: \n");
  printf ("\n");

  /*
    ii = 6*FMASS;
    for(i=FMASS; i<IAROIDS; i++)
    */

  ii = 0;
  for (i = 0; i <= ISUN; i++)

  {
    /* Compare Heliocentric positions
      xa = yn0[ii+1] - yn0[6*ISUN+1];
      ya = yn0[ii+3] - yn0[6*ISUN+3];
      za = yn0[ii+5] - yn0[6*ISUN+5];
      xb = yn1[ii+1] - yn1[6*ISUN+1];
      yb = yn1[ii+3] - yn1[6*ISUN+3];
      zb = yn1[ii+5] - yn1[6*ISUN+5];
      */

    if( i == 0 )
    {
      dx = yn0[0];
      dy = yn0[2];
      dz = yn0[4];
      printf ("%22.15lf %22.15lf %22.15lf \n", dx, dy, dz );
    }

    // посчитаем разницу (ошибку) между двумя моментами времени:
    // "расчетный - табуилрованный"

    xa = yn0[ii+1] - yn1[ii+1];
    ya = yn0[ii+3] - yn1[ii+3];
    za = yn0[ii+5] - yn1[ii+5];
    err = SQRT( xa*xa + ya*ya + za*za );

    /* Use ieee.c to print numerical results */

#if LDOUBLE
    funion.ld = yn0[ii+1];
    e64toasc( funion.s, fname, 15 );
    printf( "%22s ", fname );

    funion.ld = yn0[ii+3];
    e64toasc( funion.s, fname, 15 );
    printf( "%22s ", fname );

    funion.ld = yn0[ii+5];
    e64toasc( funion.s, fname, 15 );
    printf( "%22s ", fname );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if( iout )
    {
      funion.ld = err;
      e64toasc( funion.s, fname, 3 );
      printf( "err= %s\n", fname );

    } else {
      printf( "\n");
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#else
    funion.d = yn0[ii+1];
    e53toasc( funion.s, fname, 15 );
    printf( "%s ", fname );

    funion.d = yn0[ii+3];
    e53toasc( funion.s, fname, 15 );
    printf( "%s ", fname );

    funion.d = yn0[ii+5];
    e53toasc( funion.s, fname, 15 );
    printf( "%s ", fname );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if( iout )
    {
      funion.d = err;
      e53toasc( funion.s, fname, 3 );
      printf( "err= %s\n", fname );
    }else {
      printf( "\n");
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

    /*
      dx = (double) yn0[ii+1];
      dy = (double) yn0[ii+3];
      dz = (double) yn0[ii+5];
      printf("%19.15lf%19.15lf%19.15lf", dx, dy, dz );
      dz = (double) err;
      printf(" %.3e\n", dz );
      */

    ii += 6;
  }

  printf ("\n");


}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
main (int argc, char *argv[])
{

  DOUBLE xa, ya, za, xb, yb, zb; /* for display */
  int i, ii, iout, restored;


  int test_num  = /* 0 */1;

  get_options (argc, argv, 
               //NULL, NULL, 
             NULL, NULL, 
             NULL, NULL, &test_num, NULL);


  switch (test_num) {
    
    //case 0:
    //  break;
    
  case 1:
    ord  =  12;       // Adams order
    dx   =  0.125;    // Step size, days
    outsteps = 320;   // Interval between output samples
    nsamps   =  10;   // Number of output samples
    break;
    
  case 2:
    ord  = 10;       
    dx   = 0.125;    
    outsteps = 100;   
    nsamps   =  10;     
    break;
    
  default:
    fprintf (stderr, "\nERROR in test_num !!!! \n\n");
    exit (0);
    break;

  }

  fprintf (stderr, "\n");
  fprintf (stderr, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n");
  fprintf (stderr, "\n");


#if FPESHOW
  fpeins (); // ?? арифметик дебагер
#endif

  malloc_alls ();

  print_consts ();


  /* Read in all physical constants from file named aconst.h.
    * This is optional, since the constants are already compiled
    * in to the program.
    */
  //rdnums ();
  //
  //
  set_const_by_hand (); // ручная подстановка:

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  restored = 0;

  fname[0] = 'n';
  printf( "... Not to restore a previously saved integrator state !! \n" );
  //  
  // Do you want to restore a previously saved integrator state?
  // If the answer is 'z' then the state will be restored and also
  // the program will ask for a new interval between output records.
  //  
  
  ii = fname[0] & 0x7f;
  if( (ii == 'y') || (ii == 'z') )
  {

    resstate ();
    restored = 1;

    goto opnout;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

orlup:

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  printf( "... Adams order = %d !! \n", ord);

  if( (ord > MAXORD) || (ord < 1) )
  {
    printf( "order must be between 1 and %d\n", MAXORD );
    goto orlup;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  printf( "... Step size, days = %lf !! \n", dx);
  
  ssyh = dx;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  printf ( "... initializing  .. " );

  /* Ensure that the relativistic center of the N bodies is at 0.0.
    * This condition is already approximately satisfied by states
    * taken from the JPL ephemeris tapes.
    */
  findcenter( JD0, yn0 );
  findcenter( JD1, yn1 );

  ssyt = JD0;
  err = Zero;
  rkstart( NEQ, rwork );

#if DEBUG
  printf( "rkstart " );
#endif

  //
  // ssyh = dx
  //
  //
  //

  adstart (ssyh /*dx*/, &yn0[0], awork, NEQ, ord, ssyt /*JD0*/ ); // считаем для набора JD0 ?

#if DEBUG
  printf( "adstart\n" );
#endif

  printf( "initialized. \n" );

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /* Each output file record contains the Julian date
    * followed by the current velocity and position states.
    * See ini118.h for the data structure of the state array.
    * The state variables for a test object such as a comet
    * or asteroid should be placed in the array immediately after
    * the Sun state.
    */
opnout:

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /* If a previous state was restored, you had the option of 
    * getting the interval between outputs from the saved state
    * file or of entering a new interval by responding 'z'.
    */

  if( ii != 'y' )
  {

    printf( "... outsteps = %d\n", outsteps ); 
    // Interval between output samples, in steps 

  }

  strcpy (fname, "FILE");

  printf( "... fname = %s \n", fname);

  printf( "... nsamps = %ld \n", nsamps); // Number of output samples 


  fd = open( fname, O_CREAT | O_TRUNC | O_RDWR, S_IWRITE | S_IREAD );

  if( fd <= 0 )
  {
    printf( "Can't create \"%s\"\n", fname );
    exit( 1);
  }



  /* Write initial state vector to output file */
  // 
  if (restored == 0) {

    wfile (fd, ssyt, yn0);
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //   РЕАЛЬНЫЙ РАСЧЕТ
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  printf ("\n");

  nsteps = 0;

  for (isamps=0; isamps < nsamps; isamps++)
  {
    for (iout = 0; iout < outsteps; iout++)
    {
      err += adstep (&ssyt, &yn0[0], NEQC); // шаги вычислений?
      //
      //  вызывается: -> func -> aroids
      // 

      nsteps += 1;
    }

    wfile (fd, ssyt, yn0);
    dx = err;
    dy = ssyt;

    printf ("%5ld %11.2lf %.6e\n", nsteps, dy, dx);
    /* #if IBMPC */


#if MSC
    /* Check for operator abort
      */
    if (kbhit())
    {
      ii = getchar();
      if( (ii & 0x7f) == 'S' )
      {
        printf ("Operator interrupt, closing %s\n", fname);
        break;
      }
    }
#endif

  }

  printf ("\n");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  close    (fd);

  savstate ();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /* Test if the ending time is equal to the time of
    * the compiled-in test state vector. If so, print out
    * the errors.
    */

  xa = ssyt - JD1;

  fprintf (stderr, "JD0= %f ssyt= %f JD1= %f xa= %f", 
           (double) JD0, (double) ssyt, 
           (double) JD1, (double) xa);
  printf ("\n");

  if ( xa < 0 ) {
    xa = -xa;
  }

  /* allow for roundoff in the time variable */

  if( xa > 1e-6 )
    iout = 0; /* error printout not appropriate */
  else
    iout = 1;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  print_final_positions_and_errors (iout);

  printf ("\n");


#if FPESHOW
  fperem ();
#endif

} /* end of main program */
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
