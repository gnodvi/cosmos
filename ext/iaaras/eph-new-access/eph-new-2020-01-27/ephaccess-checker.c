#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#include "libephaccess/ephaccess.h"
#include "ephaccess.h"

const double jd_start = 2415020.5; // 01.01.1900
const double jd_end   = 2488069.5; // 01.01.2100

int testBodyLongInterval (EphAccess *eph,
        const char *body, const char *reference,
        double estimated_distance, double estimated_velocity)
{
  int id_object = ephObjectByName(body);
  int id_reference = ephObjectByName(reference);
  double jd;
  double pos[3], vel[3];

  if (id_object < 0)
  {
    fprintf(stderr, "unknown body %s\n", body);
    return -1;
  }
  if (id_reference < 0)
  {
    fprintf(stderr, "unknown reference %s\n", reference);
    return -1;
  }
  
  // 01.01.1900 -- 01.01.2100
  for (jd = jd_start; jd < jd_end; jd += 1.0)
  {
    double r, v;
    if (ephCalculateRectangular(eph, id_object, id_reference, jd, 0.0,
                                pos, vel))
    {
      fprintf(stderr, "%s\n", ephLastError(eph));
      return -1;
    }
    r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    v = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    
    if (fabs(r - estimated_distance) > estimated_distance * 0.25)
    {
      fprintf(stderr, "trouble with %s: on JD %lf distance = %lf, expected %lf +- 25%%\n",
              body, jd, r, estimated_distance);
      return -1;
    }
    if (fabs(v - estimated_velocity) > estimated_velocity * 0.31)
    {
      fprintf(stderr, "trouble with %s: on JD %lf velocity = %lf, expected %lf +- 31%%\n",
              body, jd, v, estimated_velocity);
      return -1;
    }
  }
  
  printf("%s checked\n", body);
  return 0;
}

int testBodySmoothness (EphAccess *eph, const char *body, const char *reference)
{
  int id_object = ephObjectByName(body);
  int id_reference = ephObjectByName(reference);
  double jd;
  double pos0[3], vel0[3], pos[3], pos_prev[3];
  int j, k;

  if (id_object < 0)
  {
    fprintf(stderr, "unknown body %s\n", body);
    return -1;
  }
  if (id_reference < 0)
  {
    fprintf(stderr, "unknown reference %s\n", reference);
    return -1;
  }
  
  for (jd = jd_start; jd < jd_end; jd += 200.0)
  {
    if (ephCalculateRectangular(eph, id_object, id_reference, jd, 0.0, pos0, vel0))
    {
      fprintf(stderr, "%s\n", ephLastError(eph));
      return -1;
    }
    
    memcpy(pos_prev, pos0, 3 * sizeof(double));
    
    // run small steps for one minute
    // one step = 1/1000 of a minute = 1/(1000*60*24) of a day
    int per_minute = 1000;
    for (k = 1; k < per_minute; k++)
    {
      if (ephCalculateRectangular(eph, id_object, id_reference, jd, k / (per_minute * 60 * 24.0),
                                  pos, NULL))
      {
        fprintf(stderr, "%s\n", ephLastError(eph));
        return -1;
      }
      
      for (j = 0; j < 3; j++)
      {
        if (fabs(((fabs(pos[j] - pos0[j]) / (k * 60.0 / per_minute) + 0.1) / 
                  (fabs(vel0[j]) + 0.1)) - 1.0) > 0.05)
        {
          fprintf(stderr, "trouble with %s: not smooth near JD %lf\n", body, jd);
          return -1;
        }
        if (fabs(((fabs(pos[j] - pos_prev[j]) / (60.0 / per_minute) + 0.1) /
                  (fabs(vel0[j]) + 0.1)) - 1.0) > 0.05)
        {
          printf("pos[%d] = %lf, pos_prev[%d] = %lf, vel = %lf, vel0 = %lf\n",
                  j, pos[j], j, pos_prev[j], (pos[j] - pos_prev[j]) / (60.0 / per_minute),
                  vel0[j]);
          fprintf(stderr, "trouble with %s: scatter near JD %lf\n", body, jd);
          return -1;
        }
      }
      
      memcpy(pos_prev, pos, 3 * sizeof(double));
    }
  }
  printf("%s smoothness OK\n", body);
  return 0;
}

int testLibrationLongInterval (EphAccess *eph)
{
  double jd, angles[3], rates[3];
  
  // 01.01.1900 -- 01.01.2100
  for (jd = jd_start; jd < jd_end; jd += 1.0)
  {
    if (ephCalculateEulerAngles(eph, 0, jd, 0, angles, rates))
    {
      fprintf(stderr, "%s\n", ephLastError(eph));
      return -1;
    }
    
    if (fabs(angles[0]) > 0.1)
    {
      fprintf(stderr, "trouble with libration: on JD %lf phi = %lf, expected |phi| < 0.1\n",
              jd, angles[0]);
      return -1;
    }
    if (fabs(angles[1]) > 0.5)
    {
      fprintf(stderr, "trouble with libration: on JD %lf theta = %lf, expected |theta| < 0.5\n",
              jd, angles[1]);
      return -1;
    }
    if (rates[2] > 2.7e-6 || rates[2] < 2.6e-6)
    {
      fprintf(stderr, "trouble with libration: on JD %lf psi' = %.10lf, expected 2.6e-6 < psi' < 2.7e-6\n",
              jd, rates[2]);
      return -1;
    }
  }
  printf("libration checked\n");
  return 0;
}

int testTDBLongInterval (EphAccess *eph)
{
  double jd, diff;
  
  // 01.01.1900 -- 01.01.2100
  for (jd = jd_start; jd < jd_end; jd += 1.0)
  {
    if (ephCalculateTimeDiff(eph, EPH_TT_MINUS_TDB, jd, 0.0, &diff))
    {
      fprintf(stderr, "%s\n", ephLastError(eph));
      return -1;
    }
    
    if (fabs(diff) > 2e-3)
    {
      fprintf(stderr, "trouble with TT-TDB: on JD %lf it is %.10lf, expected |TT-TDB| < 2e-3\n",
              jd, diff);
      return -1;
    }
  }
  printf("TT-TDB checked\n");
  return 0;
}

int _run (EphAccess *eph, int argc, char **argv)
{
  int i;
  
  if (argc < 2)
  {
    fprintf(stderr, "usage: %s <ephemeris file 1> [<ephemeris file 2> ...]\n",
	    argv[0]);
    return -1;
  }
  
  for (i = 1; i < argc; i++)
    if (ephLoadFile(eph, argv[i]))
    {
      fprintf(stderr, "%s\n", ephLastError(eph));
      return -1;
    }

  
  if (testBodyLongInterval(eph, "earth", "sun", 150e6, 30.0)) return -1;
  if (testBodyLongInterval(eph, "moon", "earth", 384000.0, 1.0)) return -1;
  if (testBodyLongInterval(eph, "mercury", "sun", 57900000, 48.0)) return -1;
  if (testBodyLongInterval(eph, "venus", "sun", 108e6, 45.0)) return -1;
  if (testBodyLongInterval(eph, "mars_bc", "sun", 228e6, 24.0)) return -1;
  if (testBodyLongInterval(eph, "jupiter_bc", "sun", 778e6, 13.0)) return -1;
  if (testBodyLongInterval(eph, "saturn_bc", "sun", 1.4e9, 10.0)) return -1;
  if (testBodyLongInterval(eph, "uranus_bc", "sun", 2.87e9, 6.8)) return -1;
  if (testBodyLongInterval(eph, "neptune_bc", "sun", 4.5e9, 5.4)) return -1;
  if (testBodyLongInterval(eph, "pluto_bc", "sun", 5.9e9, 4.67)) return -1;
  
  if (testLibrationLongInterval(eph)) return -1;
  if (testTDBLongInterval(eph)) return -1;
  
  if (testBodySmoothness(eph, "sun", "ssb")) return -1;
  if (testBodySmoothness(eph, "earth", "sun")) return -1;
  if (testBodySmoothness(eph, "moon", "earth")) return -1;
  if (testBodySmoothness(eph, "mercury", "sun")) return -1;
  if (testBodySmoothness(eph, "mars_bc", "sun")) return -1;
  if (testBodySmoothness(eph, "jupiter_bc", "sun")) return -1;
  if (testBodySmoothness(eph, "saturn_bc", "sun")) return -1;
  if (testBodySmoothness(eph, "uranus_bc", "sun")) return -1;
  if (testBodySmoothness(eph, "neptune_bc", "sun")) return -1;
  if (testBodySmoothness(eph, "pluto_bc", "sun")) return -1;

  return 0;
}

int main(int argc, char** argv)
{
  EphAccess *eph = ephCreate();
  int rc;
  
  if (eph == NULL)
  {
    fprintf(stderr, "initialization error\n");
    return -1;
  }

  rc = _run(eph, argc, argv);
  
  ephDestroy(eph);
  
  return rc;
}
