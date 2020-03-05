//------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#include "epm_bin_access.h"

#ifdef _WIN32
#define snprintf _snprintf
#define strcasecmp _stricmp
#endif

typedef struct
{
  char directory[1024];
  char suffix[100];
  double mjd;
  double interval;
  double step;
  int width;
  int positions_only;
  double dist_scale;
  double time_scale;
  int tt;
} Config;
//------------------------------------------------------------------------------

void init_config (Config *config)
{
  memset(config, 0, sizeof(Config));
  config->width = 20;
  config->step = 1.0;
  config->dist_scale = 1.0;
  config->time_scale = 1.0;
}
//------------------------------------------------------------------------------

void mjd_to_jd_plus_frac (double mjd, double corr, int *jd, double *frac)
{
  *jd = (int)mjd + 2400000;
  *frac = mjd + 2400000 - *jd + 0.5;

  (*frac) += corr;
  
  if (*frac >= 1.0)
  {
    (*frac) -= 1.0;
    (*jd) += 1;
  }
  else if (*frac < 0)
  {
    (*frac) += 1.0;
    (*jd) -= 1;
  }
}
//------------------------------------------------------------------------------

void list_tdb (Config *config, const char *object)
{
  char full_name[1024];
  EpmFile *epmf;
  double mjd;
  double pos[3], vel[3];

  if (snprintf(full_name, sizeof(full_name), "%s/%s.%s",
               config->directory, object, config->suffix) >=
      sizeof(full_name))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }
  
  epmf = epm_file_open(full_name);

  if (epmf == NULL)
    return;
  
  for (mjd = config->mjd;
       mjd <= config->mjd + config->interval;
       mjd += config->step)
  {
    int jd;
    double frac;
    
    mjd_to_jd_plus_frac(mjd, 0.0, &jd, &frac);
    
    if (!epm_calc_raw_values(epmf, jd, frac, pos, vel))
      return;
    
    printf("%17.8lf %+*.*le\n", //,
           mjd,
           config->width, config->width - 7, vel[0] * config->time_scale / 24 / 3600.0);
  }
  epm_file_close(epmf);
}
//------------------------------------------------------------------------------

void list_one (Config *config, const char *object)
{
  char full_name[1024];
  EpmFile *epmf, *epmf_tdb;
  double mjd;
  double pos[3], vel[3];

  if (snprintf(full_name, sizeof(full_name), "%s/%s.%s",
               config->directory, object, config->suffix) >=
      sizeof(full_name))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }
  
  epmf = epm_file_open(full_name);

  if (epmf == NULL)
    return;

  if (config->tt)
  {
    if (snprintf(full_name, sizeof(full_name), "%s/tdb.%s",
                 config->directory, config->suffix) >=
        sizeof(full_name))
    {
      fprintf(stderr, "TDB file name too long\n");
      return;
    }
    
    epmf_tdb = epm_file_open(full_name);

    if (epmf_tdb == NULL)
    {
      epm_file_close(epmf);
      return;
    }
  }
  
  for (mjd = config->mjd;
       mjd <= config->mjd + config->interval;
       mjd += config->step)
  {
    int jd;
    double frac;
    
    mjd_to_jd_plus_frac(mjd, 0.0, &jd, &frac);
    
    if (config->tt)
    {
      if (!epm_calc_raw_values(epmf_tdb, jd, frac, pos, vel))
        return;
      mjd_to_jd_plus_frac(mjd, -vel[0] / 24.0 / 3600.0, &jd, &frac);
    }
    
    if (!epm_calc_raw_values(epmf, jd, frac, pos, vel))
      return;
    
    printf("%17.8lf %+*.*le %+*.*le %+*.*le", //,
           mjd,
           config->width, config->width - 7, pos[0] * config->dist_scale,
           config->width, config->width - 7, pos[1] * config->dist_scale,
           config->width, config->width - 7, pos[2] * config->dist_scale);
    if (config->positions_only)
      printf("\n");
    else
    {
      printf(" %+*.*le %+*.*le %+*.*le\n",
             config->width, config->width - 7,
             vel[0] * config->dist_scale / config->time_scale,
             config->width, config->width - 7,
             vel[1] * config->dist_scale / config->time_scale,
             config->width, config->width - 7,
             vel[2] * config->dist_scale / config->time_scale);
    }
  }
  epm_file_close(epmf);
}
//------------------------------------------------------------------------------

void list_two (Config *config, const char *object1, const char *object2,
               double coef1, double coef2)
{
  char full_name1[1024], full_name2[1024];
  double mjd;
  EpmFile *epmf1, *epmf2, *epmf_tdb;
  double pos1[3], vel1[3], pos2[3], vel2[3];

  if (snprintf(full_name1, sizeof(full_name1), "%s/%s.%s",
               config->directory, object1, config->suffix) >=
      sizeof(full_name1))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }

  if (snprintf(full_name2, sizeof(full_name2), "%s/%s.%s",
               config->directory, object2, config->suffix) >=
      sizeof(full_name2))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }
  
  epmf1 = epm_file_open(full_name1);
  if (epmf1 == NULL)
    return;

  epmf2 = epm_file_open(full_name2);
  if (epmf2 == NULL)
  {
    epm_file_close(epmf1);
    return;
  }

  if (config->tt)
  {
    if (snprintf(full_name1, sizeof(full_name1), "%s/tdb.%s",
                 config->directory, config->suffix) >=
        sizeof(full_name1))
    {
      fprintf(stderr, "TDB file name too long\n");
      epm_file_close(epmf1);
      epm_file_close(epmf2);
      return;
    }
    
    epmf_tdb = epm_file_open(full_name1);

    if (epmf_tdb == NULL)
    {
      epm_file_close(epmf1);
      epm_file_close(epmf2);
      return;
    }
  }
  
  
  for (mjd = config->mjd;
       mjd <= config->mjd + config->interval;
       mjd += config->step)
  {
    int jd;
    double frac;
    
    mjd_to_jd_plus_frac(mjd, 0.0, &jd, &frac);
    
    if (config->tt)
    {
      if (!epm_calc_raw_values(epmf_tdb, jd, frac, pos1, vel1))
        return;
      mjd_to_jd_plus_frac(mjd, -vel1[0] / 24.0 / 3600.0, &jd, &frac);
    }
    
    if (!epm_calc_raw_values(epmf1, jd, frac, pos1, vel1))
      return;
    if (!epm_calc_raw_values(epmf2, jd, frac, pos2, vel2))
      return;

    printf("%17.8lf %+*.*le %+*.*le %+*.*le", //,
           mjd,
           config->width, config->width - 7,
           (pos1[0] * coef1 + pos2[0] * coef2) * config->dist_scale,
           config->width, config->width - 7,
           (pos1[1] * coef1 + pos2[1] * coef2) * config->dist_scale,
           config->width, config->width - 7,
           (pos1[2] * coef1 + pos2[2] * coef2) * config->dist_scale);
    if (config->positions_only)
      printf("\n");
    else
    {
      printf(" %+*.*le %+*.*le %+*.*le\n",
             config->width, config->width - 7,
             (vel1[0] * coef1 + vel2[0] * coef2) * config->dist_scale / config->time_scale,
             config->width, config->width - 7,
             (vel1[1] * coef1 + vel2[1] * coef2) * config->dist_scale / config->time_scale,
             config->width, config->width - 7,
             (vel1[2] * coef1 + vel2[2] * coef2) * config->dist_scale / config->time_scale);
    }
  }
  epm_file_close(epmf1);
  epm_file_close(epmf2);
}
//------------------------------------------------------------------------------

void list_three (Config *config,
                 const char *object1, const char *object2, const char *object3,
                 double coef1, double coef2, double coef3)
{
  char full_name1[1024], full_name2[1024], full_name3[1024];
  double mjd;
  EpmFile *epmf1, *epmf2, *epmf3, *epmf_tdb;
  double pos1[3], vel1[3], pos2[3], vel2[3], pos3[3], vel3[3];

  if (snprintf(full_name1, sizeof(full_name1), "%s/%s.%s",
               config->directory, object1, config->suffix) >=
      sizeof(full_name1))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }

  if (snprintf(full_name2, sizeof(full_name2), "%s/%s.%s",
               config->directory, object2, config->suffix) >=
      sizeof(full_name2))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }

  if (snprintf(full_name3, sizeof(full_name3), "%s/%s.%s",
               config->directory, object3, config->suffix) >=
      sizeof(full_name3))
  {
    fprintf(stderr, "file name too long\n");
    return;
  }
  
  epmf1 = epm_file_open(full_name1);
  if (epmf1 == NULL)
    return;

  epmf2 = epm_file_open(full_name2);
  if (epmf2 == NULL)
  {
    epm_file_close(epmf1);
    return;
  }

  epmf3 = epm_file_open(full_name3);
  if (epmf3 == NULL)
  {
    epm_file_close(epmf1);
    epm_file_close(epmf2);
    return;
  }

  if (config->tt)
  {
    if (snprintf(full_name1, sizeof(full_name1), "%s/tdb.%s",
                 config->directory, config->suffix) >=
        sizeof(full_name1))
    {
      fprintf(stderr, "TDB file name too long\n");
      return;
    }
    
    epmf_tdb = epm_file_open(full_name1);

    if (epmf_tdb == NULL)
    {
      epm_file_close(epmf1);
      epm_file_close(epmf2);
      epm_file_close(epmf3);
      return;
    }
  }
  
  for (mjd = config->mjd;
       mjd <= config->mjd + config->interval;
       mjd += config->step)
  {
    int jd;
    double frac;
    
    mjd_to_jd_plus_frac(mjd, 0.0, &jd, &frac);
    
    if (config->tt)
    {
      if (!epm_calc_raw_values(epmf_tdb, jd, frac, pos1, vel1))
        return;
      mjd_to_jd_plus_frac(mjd, -vel1[0] / 24.0 / 3600.0, &jd, &frac);
    }
    
    if (!epm_calc_raw_values(epmf1, jd, frac, pos1, vel1))
      return;
    if (!epm_calc_raw_values(epmf2, jd, frac, pos2, vel2))
      return;
    if (!epm_calc_raw_values(epmf3, jd, frac, pos3, vel3))
      return;

    printf("%17.8lf %+*.*le %+*.*le %+*.*le",
           mjd,
           config->width, config->width - 7,
           (pos1[0] * coef1 + pos2[0] * coef2 + pos3[0] * coef3) * config->dist_scale,
           config->width, config->width - 7,
           (pos1[1] * coef1 + pos2[1] * coef2 + pos3[1] * coef3) * config->dist_scale,
           config->width, config->width - 7,
           (pos1[2] * coef1 + pos2[2] * coef2 + pos3[2] * coef3) * config->dist_scale);
    if (config->positions_only)
      printf("\n");
    else
    {
      printf(" %+*.*le %+*.*le %+*.*le\n",
             config->width, config->width - 7,
             (vel1[0] * coef1 + vel2[0] * coef2 + vel3[0] * coef3) *
             config->dist_scale / config->time_scale,
             config->width, config->width - 7,
             (vel1[1] * coef1 + vel2[1] * coef2 + vel3[1] * coef3) *
             config->dist_scale / config->time_scale,
             config->width, config->width - 7,
             (vel1[2] * coef1 + vel2[2] * coef2 + vel3[2] * coef3) *
             config->dist_scale / config->time_scale);
    }
  }
  epm_file_close(epmf1);
  epm_file_close(epmf2);
  epm_file_close(epmf3);
}
//------------------------------------------------------------------------------

void display_help ()
{
  fprintf(stderr,
          "Usage: epm-calculator [parameters]\n"
          "Obligatory parameters:\n"
          "   -from-mjd <real>\n"
          "      Modified Julian date to start from\n"
          "   -body [tdb | earth | moon | moonlibr | EMB | sun | venus | ... ]\n"
          "      Celestial body or another object:\n"
          "        tdb = TT-TDB difference\n"
          "        moonlibr = lunar libration angles (in radians) and rates\n"
          "        EMB = Earth-Moon barycenter\n"
          "   -interval-days <real>\n"
          "      Interval length in days\n"
          "Optional parameters:\n"
          "   -step-days <real>\n"
          "      Step size in days (default 1 day)\n"
          "   -theory [epm2008|epm2011m|epm2015]\n"
          "      Theory to use for the calculations (default epm2015)\n"
          "   -frame [geocentric | barycentric | heliocentric]\n"
          "      Reference frame (default barycentric)\n"
          "   -timescale [tt | tdb]\n"
          "      Time scale (default TDB)\n"
          "   -width <integer>\n"
          "      Width of output numbers, including sign, decimal point, and exponent (default 20)\n"
          "   -distance-units [au|km]\n"
          "      Output positions in AU or km (default AU)\n"
          "   -time-units [sec|day]\n"
          "      Output velocities as a unit of distance per second or per day (default day)\n"
          "   -only-pos\n"
          "      Output only positions, skip velocities\n"
          "   -directory <path>\n"
          "      Directory containing the binary ephemeris files; default EPM2008/BIN,\n"
          "      EPM2011m/BIN, or EPM2015/BIN depending of chosen theory\n"
          );
}
//------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
  int i;
  char *body = NULL;
  const char *frame = "barycentric";
  int theory = 2015;
  Config config;
  double rho;
  int km = 0;

  init_config(&config);
  config.mjd = 1e15;
  config.interval = 0.0;
  
  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-body") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting body name after -body\n");
        return -1;
      }
      body = argv[i];
    }
    else if (strcmp(argv[i], "-from-mjd") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting MJD after -from-mjd\n");
        return -1;
      }
      if (sscanf(argv[i], "%lf", &config.mjd) != 1)
      {
        fprintf(stderr, "%s is not a valid MJD\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-interval-days") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting number of days after -interval-days\n");
        return -1;
      }
      if (sscanf(argv[i], "%lf", &config.interval) != 1)
      {
        fprintf(stderr, "%s is not a valid number\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-directory") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting directory after -directory\n");
        return -1;
      }
      strncpy(config.directory, argv[i], sizeof(config.directory));
    }
    else if (strcmp(argv[i], "-step-days") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting step size in days after -step-days\n");
        return -1;
      }
      if (sscanf(argv[i], "%lf", &config.step) != 1)
      {
        fprintf(stderr, "%s is not a valid number\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-frame") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting reference frame after -frame\n");
        return -1;
      }
      frame = argv[i];
    }
    else if (strcmp(argv[i], "-timescale") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting time scale after -timescale\n");
        return -1;
      }
      if (strcmp(argv[i], "tt") == 0)
        config.tt = 1;
      else if (strcmp(argv[i], "tdb") == 0)
        config.tt = 0;
      else
      {
        fprintf(stderr, "unknown timescale: %s\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-width") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting integer after -width\n");
        return -1;
      }
      if (sscanf(argv[i], "%d", &config.width) != 1)
      {
        fprintf(stderr, "%s is not a valid integer\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-distance-units") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting au or km after -distance-units\n");
        return -1;
      }
      if (strcmp(argv[i], "au") == 0)
        ;
      else if (strcmp(argv[i], "km") == 0)
        km = 1;
      else
      {
        fprintf(stderr, "unknown distance units: %s\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-time-units") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting sec or day after -time-units\n");
        return -1;
      }
      if (strcmp(argv[i], "sec") == 0)
        config.time_scale = 24 * 3600;
      else if (strcmp(argv[i], "day") == 0)
        ;
      else
      {
        fprintf(stderr, "unknown time units: %s\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-theory") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting theory after -theory\n");
        return -1;
      }
      if (strcasecmp(argv[i], "epm2008") == 0)
        theory = 2008;
      else if (strcasecmp(argv[i], "epm2011m") == 0)
        theory = 2011;
      else if (strcasecmp(argv[i], "epm2015") == 0)
        theory = 2015;
      else
      {
        fprintf(stderr, "unknown theory: %s\n", argv[i]);
        return -1;
      }
    }
    else if (strcmp(argv[i], "-only-pos") == 0)
    {
      config.positions_only = 1;
    }
    else
    {
      fprintf(stderr, "unknown argument: %s\n", argv[i]);
      display_help();
      return -1;
    }
  }
 
  if (body == NULL)
  {
    fprintf(stderr, "body not specified\n");
    display_help();
    return -1;
  }

  if (config.mjd > 1e14)
  {
    fprintf(stderr, "MJD not specified\n");
    display_help();
    return -1;
  }

  if (config.interval < 1e-15)
  {
    fprintf(stderr, "interval not specified\n");
    display_help();
    return -1;
  }

  if (theory == 2008)
  {
    rho = 82.3005676536174207;
    if (km)
      config.dist_scale = 149597870.696547;
    if (config.directory[0] == 0)
      strncpy(config.directory, "EPM2008/BIN", sizeof(config.directory));
    strncpy(config.suffix, "08b", sizeof(config.suffix));
  }
  else if (theory == 2011)
  {
    rho = 82.3005676344;
    if (km)
      config.dist_scale = 149597870.69588;
    if (config.directory[0] == 0)
      strncpy(config.directory, "EPM2011m/BIN", sizeof(config.directory));
    strncpy(config.suffix, "bin", sizeof(config.suffix));
  }
  else if (theory == 2015)
  {
    rho = 82.30056888602054;
    if (km)
      config.dist_scale = 149597870.7;
    if (config.directory[0] == 0)
      strncpy(config.directory, "EPM2015/BIN", sizeof(config.directory));
    strncpy(config.suffix, "bin", sizeof(config.suffix));
  }
  else
  {
    fprintf(stderr, "internal error\n");
    return -1;
  }

  for (i = 0; body[i] != 0; i++)
    body[i] = tolower(body[i]);

  if (strcasecmp(body, "tdb") == 0)
  {
    list_tdb(&config, "tdb");
  }
  else if (strcasecmp(body, "moonlibr") == 0)
  {
    config.dist_scale = 1.0;
    list_one(&config, "moonlibr");
  }
  else if (strcasecmp(body, "moon") == 0)
  {
    if (strcasecmp(frame, "geocentric") == 0)
      list_one(&config, "moon");
    else if (strcasecmp(frame, "barycentric") == 0)
      // barycentric Moon = EMB + (geocentric Moon) *
      //                    (mass of Earth) / (mass of Earth + mass of Moon)
      //                  = EMB + (geocentric Moon) * (rho - 1) / rho
      list_two(&config, "earth_m", "moon", 1.0, (rho - 1.0) / rho);
    else if (strcasecmp(frame, "heliocentric") == 0)
      // heliocentric Moon = barycentric Moon - barycentric Sun
      list_three(&config, "earth_m", "moon", "sun", 1.0, (rho - 1.0) / rho, -1.0);
    else
    {
      fprintf(stderr, "unsupported frame: %s\n", frame);
      return -1;
    }
  }
  else if (strcasecmp(body, "earth") == 0)
  {
    if (strcasecmp(frame, "geocentric") == 0)
    {
      fprintf(stderr, "geocentric Earth is an interesting choice\n");
      return -1;
    }
    else if (strcasecmp(frame, "barycentric") == 0)
    {
      // barycentric Earth = EMB - (geocentric Moon) *
      //                    (mass of Moon) / (mass of Earth + mass of Moon)
      //                   = EMB - (geocentric Moon) / rho
      list_two(&config, "earth_m", "moon", 1.0, -1.0 / rho);
    }
    else if (strcasecmp(frame, "heliocentric") == 0)
      // heliocentric Earth = barycentric Earth - barycentric Sun
      list_three(&config, "earth_m", "moon", "sun", 1.0, -1.0 / rho, -1.0);
    else
    {
      fprintf(stderr, "unsupported frame: %s\n", frame);
      return -1;
    }
  }
  else
  {
    if (strcasecmp(body, "EMB") == 0)
      body = "earth_m";
    
    if (strcasecmp(frame, "geocentric") == 0)
      // geocentric body = barycentric body - barycentric Earth =
      //                   barycentric body - EMB + (geocentric Moon) / rho
      list_three(&config, body, "earth_m", "moon", 1.0, -1.0, 1.0 / rho);
    else if (strcasecmp(frame, "heliocentric") == 0)
    {
      if (strcasecmp(body, "sun") == 0)
      {
        fprintf(stderr, "heliocentric Sun is an interesting choice\n");
        return -1;
      }
      // heliocentric body = barycentric body - barycentric Sun
      list_two(&config, body, "sun", 1.0, -1.0);
    }
    else if (strcasecmp(frame, "barycentric") == 0)
      list_one(&config, body);
    else
    {
      fprintf(stderr, "unsupported frame: %s\n", frame);
      return -1;
    }
  }

  return 0;
}
//------------------------------------------------------------------------------
