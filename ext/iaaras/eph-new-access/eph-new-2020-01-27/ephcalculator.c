/* 
 * Файл         : ephcalculator.c
 * Исполнитель  : ИПА РАН
 * Заказчик     : ОАО НПК СПП
 * Проект       : КР ФЭЛП, ОКР "Эфемериды"
 * Год          : 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ephaccess.h"
//#include "libephaccess/ephaccess.h"

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define snprintf _snprintf
#endif

typedef struct
{
  double jd0, jd1;
  double interval;
  double step;
  int width;
  int positions_only;
} Config;

void init_config (Config *config)
{
  memset(config, 0, sizeof(Config));
  config->width = 20;
  config->step = 1.0;
}

void display_help ()
{
  fprintf(stderr,
          "Usage: ephcalculator [parameters]\n"
          "Obligatory parameters:\n"
          "   -from-jd <real> <real>\n"
          "      Julian date to start from (two given numbers are added up)\n"
          "   -object [ moonlibr | tt-tdb | earth | moon | EMB | sun | venus |\n"
          "             mars_bc | jupiter_bc | ...]\n"
          "      Object of interest\n"
          "   -reference [ EMB | earth | sun | SSB | ...]\n"
          "      Reference object (does not apply to moonlibr or tt-tdb)\n"
          "   -interval-days <real>\n"
          "      Interval length in days\n"
          "   -file <string>\n"
          "      File name. Multiple files are allowed, each with its own -file.\n"
          "      Allowed formats are:\n"
          "      SPK (.bsp extension), binary PCK (.bpc)\n"
          "Optional parameters:\n"
          "   -step-days <real>\n"
          "      Step size in days (default 1 day)\n"
          "   -width <integer>\n"
          "      Width of output numbers, including sign, decimal point, and exponent (default 20)\n"
          "   -distance-units [au|km]\n"
          "      Output positions in AU or km (default km)\n"
          "   -time-units [sec|day]\n"
          "      Output velocities as a unit of distance per second or per day (default day)\n"
          "   -only-pos\n"
          "      Output only positions, skip velocities\n"
          );
}

void _date_to_text (double jd0, double jd1, char jd_text[64])
{
  int jd0_int = (int)jd0;
  double jd0_frac = jd0 - jd0_int;
  int jd1_int = (int)jd1;
  double jd1_frac = jd1 - jd1_int;
  
  int sum_int = jd0_int + jd1_int;
  double sum_frac = jd0_frac + jd1_frac;
  
  if (sum_frac > 1.0)
  {
    sum_int += 1;
    sum_frac -= 1.0;
  }
  if (sum_frac < 0.0)
  {
    sum_int -= 1;
    sum_frac += 1.0;
  }
  
  snprintf(jd_text, 64, "%d.%07d%07d", sum_int,
          (int)(sum_frac * 1e7),
          (int)((sum_frac - ((int)(sum_frac * 1e7) * 1e-7)) * 1e14));
}

#define CHECK(call) do { \
   if ((rc = call) != 0) { \
     fprintf(stderr, "%s\n", ephLastError(eph)); \
     return rc; \
   } } while (0)
  

int _printLibrations (EphAccess *eph, Config *config)
{
  double ang[6];
  double jd;
  char jd_text[64];
  int rc;

  for (jd = config->jd1;
       jd <= config->jd1 + config->interval;
       jd += config->step)
  {
    _date_to_text(config->jd0, jd, jd_text);

    CHECK(ephCalculateEulerAngles(eph, 0, config->jd0, jd, ang,
            config->positions_only ? NULL : ang + 3));

    printf("%s %+*.*le %+*.*le %+*.*le", //,
           jd_text,
           config->width, config->width - 7, ang[0],
           config->width, config->width - 7, ang[1],
           config->width, config->width - 7, ang[2]);
    if (config->positions_only)
      printf("\n");
    else
    {
      printf(" %+*.*le %+*.*le %+*.*le\n",
             config->width, config->width - 7, ang[3],
             config->width, config->width - 7, ang[4],
             config->width, config->width - 7, ang[5]);
    }
  }
  ephDestroy(eph);
  return 0;
}

int _printTTTDB (EphAccess *eph, Config *config)
{
  double jd;
  char jd_text[64];
  int rc;

  for (jd = config->jd1;
       jd <= config->jd1 + config->interval;
       jd += config->step)
  {
    double diff;
    _date_to_text(config->jd0, jd, jd_text);

    CHECK(ephCalculateTimeDiff(eph, EPH_TT_MINUS_TDB,
            config->jd0, jd, &diff));

    printf("%s %+*.*le\n", jd_text, config->width, config->width - 7, diff);
  }
  ephDestroy(eph);
  return 0;
}

int _printRectangular (EphAccess *eph, Config *config, int object, int reference)
{
  double posvel[6];
  double jd;
  char jd_text[64];
  int rc;

  for (jd = config->jd1;
       jd <= config->jd1 + config->interval;
       jd += config->step)
  {
    _date_to_text(config->jd0, jd, jd_text);

    CHECK(ephCalculateRectangular(eph, object, reference, config->jd0, jd,
            posvel, config->positions_only ? NULL : posvel + 3));

    printf("%s %+*.*le %+*.*le %+*.*le", //,
           jd_text,
           config->width, config->width - 7, posvel[0],
           config->width, config->width - 7, posvel[1],
           config->width, config->width - 7, posvel[2]);
    if (config->positions_only)
      printf("\n");
    else
    {
      printf(" %+*.*le %+*.*le %+*.*le\n",
             config->width, config->width - 7, posvel[3],
             config->width, config->width - 7, posvel[4],
             config->width, config->width - 7, posvel[5]);
    }
  }
  ephDestroy(eph);
  return 0;
}

int main(int argc, char *argv[])
{
  EphAccess *eph = ephCreate();
  Config config;
  int rc, i;
  char *object = NULL;
  char *reference = NULL;
  
  if (eph == NULL)
  {
    // no memory? should not happen
    fprintf(stderr, "initialization error\n");
    return -1;
  }

  init_config(&config);
  config.jd0 = 1e15;
  config.interval = 0.0;
  
  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-object") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting object name after -object\n");
        return -1;
      }
      object = argv[i];
    }
    else if (strcmp(argv[i], "-reference") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting reference name after -reference\n");
        return -1;
      }
      reference = argv[i];
    } 
    else if (strcmp(argv[i], "-file") == 0)
    {
      if (++i == argc)
      {
        fprintf(stderr, "expecting file name after -file\n");
        return -1;
      }
      CHECK(ephLoadFile(eph, argv[i]));
    }
    else if (strcmp(argv[i], "-from-jd") == 0)
    {
      if (i + 2 >= argc)
      {
        fprintf(stderr, "expecting two numbers after -from-jd\n");
        return -1;
      }
      i++;
      if (sscanf(argv[i], "%lf", &config.jd0) != 1)
      {
        fprintf(stderr, "%s is not a valid JD\n", argv[i]);
        return -1;
      }
      i++;
      if (sscanf(argv[i], "%lf", &config.jd1) != 1)
      {
        fprintf(stderr, "%s is not a valid JD\n", argv[i]);
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
      if (strcasecmp(argv[i], "km") == 0)
        ephSetDistanceUnits(eph, EPH_KM);
      else if (strcasecmp(argv[i], "au") == 0)
        ephSetDistanceUnits(eph, EPH_AU);
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
      if (strcasecmp(argv[i], "sec") == 0)
        ephSetTimeUnits(eph, EPH_SEC);
      else if (strcasecmp(argv[i], "day") == 0)
        ephSetTimeUnits(eph, EPH_DAY);
      else
      {
        fprintf(stderr, "unknown time units: %s\n", argv[i]);
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
 
  if (object == NULL)
  {
    fprintf(stderr, "object not specified\n");
    display_help();
    return -1;
  }
  
  if (config.jd0 > 1e14)
  {
    fprintf(stderr, "JD not specified\n");
    display_help();
    return -1;
  }

  if (config.interval < 1e-15)
  {
    fprintf(stderr, "interval not specified\n");
    display_help();
    return -1;
  }
  
  if (strcasecmp(object, "moonlibr") == 0)
    return _printLibrations(eph, &config);
  
  if (strcasecmp(object, "tt-tdb") == 0)
    return _printTTTDB(eph, &config);
  
  if (reference == NULL)
  {
    fprintf(stderr, "reference not specified\n");
    display_help();
    return -1;
  }
  
  {
    int id_object = ephObjectByName(object);
    int id_reference = ephObjectByName(reference);
    
    if (id_object < 0)
    {
      fprintf(stderr, "unknown object %s\n", object);
      return -1;
    }
    if (id_reference < 0)
    {
      fprintf(stderr, "unknown reference %s\n", reference);
      return -1;
    }
    return _printRectangular(eph, &config, id_object, id_reference);
  }
}
