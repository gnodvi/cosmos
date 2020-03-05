/*
 * readfile.c
 * readstar routine
 * StarChart Version 3.2 copyright (c) March 1990 by Craig Counterman 
 * original StarChart Software Suite copyright (c) 1987 by Alan Paeth
 *
 * All rights reserved. Redistribution granted for non-commercial
 * non-profit use only. Disclaimer: users of this work understand that
 * (a) the authors' cannot undertake to support this software (b) users
 * agree to acknowledge the use of the software in any published work
 * arising from its application and (c) any subsequent redistribution of
 * this work retains this warranty placard. 
 *
 * No representation is made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty, to the extent permitted by applicable law.
 *
 */

static char rcsid[]="$Header: /cvsroot/starchart/starchart/readfile.c,v 1.1.1.1 2000/07/27 19:00:12 gfk Exp $";

#include <stdio.h>
#include <math.h>

char *malloc();


#ifndef SYSV
#include <strings.h>
#else
#include <string.h>
#endif
#include <ctype.h>

#include "star3.h"

#ifndef READMODE
#define READMODE "r"
#endif
#define OPENFAIL 0
#define LINELEN 82


/* PI / 180 = .0174532925199 */
#define DCOS(x) (cos((x)*.0174532925199))
#define DSIN(x) (sin((x)*.0174532925199))
#define DTAN(x) (tan((x)*.0174532925199))
#define DASIN(x) (asin(x)/.0174532925199)
#define DATAN2(x,y) (atan2(x,y)/.0174532925199)
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))


/* read file function indirection to allow for use of alternate file readers */
extern int (*readfile)();
extern char *cur_file_name;

/* Readstar globals */
extern double obj_lat, obj_lon, obj_mag;
extern char obj_type[3], obj_color[3], obj_label[3];
extern char obj_constell[4], obj_name[LINELEN];
extern char *obj_commnt, fileline[LINELEN];


/* constellation abbreviations */
char *con_table[] = {
  "   ",
  "AND",
  "ANT",
  "APS",
  "AQL",
  "AQR",
  "ARA",
  "ARI",
  "AUR",
  "BOO",
  "CAE",
  "CAM",
  "CAP",
  "CAR",
  "CAS",
  "CEN",
  "CEP",
  "CET",
  "CHA",
  "CIR",
  "CMA",
  "CMI",
  "CNC",
  "COL",
  "COM",
  "CRA",
  "CRB",
  "CRT",
  "CRU",
  "CRV",
  "CVN",
  "CYG",
  "DEL",
  "DOR",
  "DRA",
  "EQU",
  "ERI",
  "FOR",
  "GEM",
  "GRU",
  "HER",
  "HOR",
  "HYA",
  "HYI",
  "IND",
  "LAC",
  "LEO",
  "LEP",
  "LIB",
  "LMI",
  "LUP",
  "LYN",
  "LYR",
  "MEN",
  "MIC",
  "MON",
  "MUS",
  "NOR",
  "OCT",
  "OPH",
  "ORI",
  "PAV",
  "PEG",
  "PER",
  "PHE",
  "PIC",
  "PSA",
  "PSC",
  "PUP",
  "PYX",
  "RET",
  "SCL",
  "SCO",
  "SCT",
  "SER",
  "SEX",
  "SGE",
  "SGR",
  "TAU",
  "TEL",
  "TRA",
  "TRI",
  "TUC",
  "UMA",
  "UMI",
  "VEL",
  "VIR",
  "VOL",
  "VUL"
  };


/* typedefs for exact sizes of int */
typedef char int_8;
typedef short int int_16;
typedef long int int_32;

/* BINFULL structure */
struct bfull_struct {
  int_32 lat;			/* RA in seconds * 1000 */
  int_32 lon;			/* Dec in seconds * 1000 */
  int_16 mag;			/* Mag * 1000 */
  char tycolb[6];		/* Type, color, label fields */
  int_8 consindx;		/* Index number of constellation */
  int_16 strlen;		/* length of name and comment field */
} binfull_in;

char name_comment[LINELEN];

/* BINOBJ structure */
struct bobj_struct {
  int_32 lat;			/* RA in seconds * 1000 */
  int_32 lon;			/* Dec in seconds * 1000 */
  int_16 mag;			/* Mag * 1000 */
  char type[2];			/* e.g. 'SD', 'CO' */
} binobj_in;

/* BINSTAR structure */
struct bstar_struct {
  int_32 lat;			/* RA in seconds * 1000 */
  int_32 lon;			/* Dec in seconds * 1000 */
  int_16 mag;			/* Mag * 1000 */
} binstar_in;

/* readstar reads from the file the information for one object, and
   loads the following variables:
double obj_lat, obj_lon, obj_mag;
char obj_type[] ="SS", obj_color[3], obj_label[3];
char obj_constell[3], obj_name[LINELEN];
char *obj_commnt, fileline[LINELEN];

only lat, lon, and mag are required.  type should default to 'SS',
color, label, constell default to "  ", and the rest default to ""
*/
int readstar(file, ftype)
     FILE *file;
     int ftype;
{
  char *ptr;
  double rah, ram, ras, dld, dlm, dl, inten;
  int i, j;
  int nchars;
  char m1;

  if ((ftype != LINEREAD) && (ftype != BINFULL)
      && (ftype != BINOBJ) && (ftype != BINSTAR) && (ftype != GSCTYPE))
    return (TRUE);
  /* only LINEREAD, BINFULL, BINOBJ, BINSTAR and GSC supported at this time */

  if (ftype == BINSTAR) {
    if (fread((char *) &binstar_in, sizeof(binstar_in), 1, file) != 1) {
      if (feof(file)) return TRUE;
      perror("Error reading input file");
      exit(2);
    }

    obj_lat = ((double) binstar_in.lat) / 3600000L;
    obj_lon = ((double) binstar_in.lon) / 3600000L;
    obj_mag = ((double) binstar_in.mag) / 1000L;
    obj_type[0] = 'S';
    obj_type[1] = 'S';
    obj_color[0] = ' ';
    obj_color[1] = ' ';
    obj_label[0] = ' ';
    obj_label[1] = ' ';
    obj_constell[0] = ' ';
    obj_constell[1] = ' ';
    obj_constell[2] = ' ';
    obj_name[0] = '\0';
    obj_commnt = "";

    strcpy(fileline, "");
  } else if (ftype == BINOBJ) {
    if (fread((char *) &binobj_in, sizeof(binobj_in), 1, file) != 1) {
      if (feof(file)) return TRUE;
      perror("Error reading input file");
      exit(2);
    }

    obj_lat = ((double) binobj_in.lat) / 3600000L;
    obj_lon = ((double) binobj_in.lon) / 3600000L;
    obj_mag = ((double) binobj_in.mag) / 1000L;
    obj_type[0] = binobj_in.type[0];
    obj_type[1] = binobj_in.type[1];
    obj_color[0] = ' ';
    obj_color[1] = ' ';
    obj_label[0] = ' ';
    obj_label[1] = ' ';
    obj_constell[0] = ' ';
    obj_constell[1] = ' ';
    obj_constell[2] = ' ';
    obj_name[0] = '\0';
    obj_commnt = "";

    strcpy(fileline, "");
  } else if (ftype == BINFULL) {
    if (fread((char *) &binfull_in, sizeof(binfull_in), 1, file) != 1) {
      if (feof(file)) return TRUE;
      perror("Error reading input file");
      exit(2);
    }

    if (binfull_in.strlen == 0)
      strcpy(name_comment, "");
    else {
      if (fread((char *) name_comment, binfull_in.strlen, 1, file) != 1) {
	perror("Error reading input file");
	exit(2);
      }
      name_comment[binfull_in.strlen] = '\0';
    }

    obj_lat = ((double) binfull_in.lat) / 3600000L;
    obj_lon = ((double) binfull_in.lon) / 3600000L;
    obj_mag = ((double) binfull_in.mag) / 1000L;
    obj_type[0] = binfull_in.tycolb[0];
    obj_type[1] = binfull_in.tycolb[1];
    obj_color[0] = binfull_in.tycolb[2];
    obj_color[1] = binfull_in.tycolb[3];
    obj_label[0] = binfull_in.tycolb[4];
    obj_label[1] = binfull_in.tycolb[5];
    strcpy(obj_constell,con_table[binfull_in.consindx]);

    ptr = name_comment;
    i = 0;
    while (*ptr == ' ') ptr++;
    while (*ptr != ',' && *ptr != '\n' && *ptr)
      obj_name[i++] = *ptr++;
    obj_name[i] = '\0';
    if ((*ptr == ',') && (*++ptr) && name_comment[0]) obj_commnt = ptr;
    else obj_commnt = "";

    strcpy(fileline, "");
  } else { /* LINEREAD */

/*
 * file formats:
 * new
064509-1643-14SDA1a CMASirius
051432-0812015SDB8b ORIRigel
 * old
064509-1643-146SSSirius
051432-08120015SSRigel
 */

    fgets(fileline, LINELEN, file);
    if (feof(file)) return(TRUE);	/* IS AN ERROR or eof */
    nchars = 0;
    while (fileline[nchars++]);
    nchars--;
    nchars--;

/*
 * sscanf of floats is TOOO slow:
 *     sscanf(fileline, "%2f%2f%2f%c%2f%2f ... );
 * use alternate:
 */
#define F2(i) (((fileline[i]-'0')*10.0+fileline[i+1]-'0'))
#define F3(i) (((fileline[i]-'0')*100.0+(fileline[i+1]-'0')*10+fileline[i+2]-'0'))
#define F4(i) (((fileline[i]-'0')*1000.0+(fileline[i+1]-'0')*100+(fileline[i+2])-'0')*10+fileline[i+3]-'0')
#define F3M(i) (((fileline[i]-'A'+10.0)*100+(fileline[i+1]-'0')*10+fileline[i+2]-'0'))
    rah = F2(0);
    ram = F2(2);
    ras = F2(4);
    dld = F2(7);
    dlm = F2(9);
/*
 * common code
 */
#define DLDEGSEC 3600.0
#define DLMINSEC 60.0
#define RAHRSSEC 54000.0
#define RAMINSEC 900.0
#define RASECSEC 15.0


    obj_lon = (RAHRSSEC*rah + RAMINSEC*ram + RASECSEC*ras)/DLDEGSEC;
    dl = (DLDEGSEC*dld + DLMINSEC*dlm)/DLDEGSEC;
    obj_lat = (fileline[6]  == '-') ? -dl : dl;
    
    /* set unknowns to blanks */
    obj_color[0] = ' ';
    obj_color[1] = ' ';
    obj_color[2] = '\0';
    obj_label[0] = ' ';
    obj_label[1] = ' ';
    obj_label[2] = '\0';
    obj_constell[0] = ' ';
    obj_constell[1] = ' ';
    obj_constell[2] = ' ';
    obj_constell[3] = '\0';
    
    if  (isdigit(fileline[14])) {
    /*
     * old reduced Yale catalog
     */
      inten = F3(12);
      if (fileline[11] == '0' || fileline[11] == '+') obj_mag = inten/100.0;
      else if (fileline[11] == '-') obj_mag = -inten/100.0;
      else obj_mag = F4(11)/1000.0;	/* new feature for stars >= 10.0 mag */
      
      if (nchars > 15) {
	obj_type[0] = fileline[15];
	obj_type[1] = fileline[16];
	ptr = &fileline[MIN(17,nchars)];
	i = 0;
	while (*ptr == ' ') ptr++;
#ifdef ATARI_ST
	while (*ptr != ',' && *ptr != '\n' && *ptr != '\r' && *ptr)
	  obj_name[i++] = *ptr++;
#else
	while (*ptr != ',' && *ptr != '\n' && *ptr)
	  obj_name[i++] = *ptr++;
#endif
	obj_name[i] = '\0';
	if (*++ptr) obj_commnt = ptr;
	else obj_commnt = "";
      } else {
	obj_type[0] = obj_type[1] = 'S'; /* Default SS single star */
	obj_name[0] = '\0';
	obj_commnt = "";
      }
    } else {
      /*
       * new reduced Yale catalog
       */
      m1 = fileline[11];
      obj_mag = ((m1 == '-') ? -F2(12)/10.0 :
		 (m1 <= '9') ? F3(11)/100.0 : F3M(11)/100.0);
      /* let's get Sirius */
    
      /*
       * extract color, label, constellation, name, and comment
       * Would be faster to just guarentee that the data file is correct
       */
      if (nchars > 22) {
	obj_constell[0] = fileline[20];
	obj_constell[1] = fileline[21];
	obj_constell[2] = fileline[22];
	obj_constell[3] = '\0';
      }
      if (nchars > 19) {
	obj_label[0] = fileline[18];
	obj_label[1] = fileline[19];
	obj_label[2] = '\0';
      }
      if (nchars > 17) {
	obj_color[0] = fileline[16]; 
	obj_color[1] = fileline[17];
	obj_color[2] = '\0';
      }
      if (nchars > 15) {
	obj_type[0] = fileline[14];
	obj_type[1] = fileline[15];
      }

      ptr = &fileline[MIN(23,nchars)];
      i = 0;
      while (*ptr == ' ') ptr++;
#ifdef ATARI_ST
      while (*ptr != ',' && *ptr != '\n' && *ptr != '\r' && *ptr)
	obj_name[i++] = *ptr++;
#else
      while (*ptr != ',' && *ptr != '\n' && *ptr)
	obj_name[i++] = *ptr++;
#endif
      obj_name[i] = '\0';
      if (*++ptr) obj_commnt = ptr;
      else obj_commnt = "";
    }
  }

  return(FALSE); /* NO error */
}

/* Macintosh, under MPW, currently won't do buf_readstar,
 atari it doesn't help much. */
#ifndef NO_BUF_READSTAR
#ifdef ATARI_ST
#include<types.h>
#include<stat.h>
#else
#include<sys/types.h>
#include<sys/stat.h>
#endif

/* > MAXMAPFILES because more files may be read
   if chart is changed interactively */
struct fbuf {
  char *name;
  char *data;
  char *cur_point;
  char *last_point;
} buf_files[10*MAXMAPFILES];

int nbuf_files = 0;

char *last_read = "";
char *cur_fpt = NULL;
char *eodata = NULL;
int cur_buf = 0;

int buf_readstar(file, ftype)
     FILE *file;
     int ftype;
{
  int i;
  struct stat fstats;
  unsigned fsze;

  char *ptr;
  double rah, ram, ras, dld, dlm, dl, inten;
  int nchars;
  char m1;

  if (ftype == GSCTYPE) return readstar(file, ftype);
  /* We don't do buffered reads of the GSC data. */

  /* if we've read this file before, access from memory
     usually will be last read
     Otherwise, try to allocate memory for file */

  if (strcmp(cur_file_name, last_read)) {	/* New file, most likely */
    for (i = 0; i < nbuf_files; i++)
      if (!strcmp(cur_file_name, buf_files[i].name)) break;

    if (i == nbuf_files) {	/* Need to read the file */
      if (buf_files[nbuf_files].name =
	  (char *) malloc((unsigned) strlen(cur_file_name)+1))
	strcpy(buf_files[nbuf_files].name, cur_file_name);
      else {
	/* out of memory */
	readfile = readstar;	/* Give up on buf_readstar */
	D_comment("buffered file read: Out of memory\n");

	/* call readstar now */
	return readstar(file, ftype);
      };

#ifndef ATARI_ST
      fstat(fileno(file), &fstats);
#else
      /* MWC 3.06 fstat is not implemented correctly. st_size is
	 always 0 (according to the manual !!) despite it's known. */
      stat(cur_file_name, &fstats);
#endif
      fsze = fstats.st_size;
      if (fsze != 0) {		/* e.g. /dev/null */
	buf_files[nbuf_files].data = (char *) malloc(fsze);
	if ((buf_files[nbuf_files].data == NULL)
	    || (((fsze = fread(buf_files[nbuf_files].data, 1,
			       (int) fsze, file)) == 0)  /* no bytes read */
		&& (!feof(file)))) { /* And it's not the end of the file */
	  readfile = readstar;	/* Give up on buf_readstar */
	  D_comment("buffered file read: Out of memory\n");

				/* call readstar now */
	  return readstar(file, ftype);
	}
      } else
	buf_files[nbuf_files].data = "";

      buf_files[nbuf_files].cur_point = buf_files[nbuf_files].data;
      cur_fpt = buf_files[nbuf_files].cur_point;
      buf_files[nbuf_files].last_point =
	buf_files[nbuf_files].cur_point + fsze;
      eodata = buf_files[nbuf_files].last_point;
      last_read = buf_files[nbuf_files].name;
      cur_buf = nbuf_files;
      nbuf_files++;
    } else { 			/* Old file */
/*      D_comment("Reading from memory\n");*/
      last_read = buf_files[i].name;
      cur_fpt = buf_files[i].cur_point = buf_files[i].data;
      eodata = buf_files[i].last_point;
      cur_buf = i;
    }
  }

  /* Now, cur_fpt points to file data.  Read the data as in readstar above */
  if ((ftype != LINEREAD) && (ftype != BINFULL)
      && (ftype != BINOBJ) && (ftype != BINSTAR))
    return (TRUE);
  /* only LINEREAD, BINFULL, BINOBJ and BINSTAR supported at this time */

  if (cur_fpt == eodata) return TRUE; /* End of data, end of file */


  if (ftype == BINSTAR) {
    for (i = 0, ptr = (char *) &binstar_in; i < sizeof(binstar_in); i++) 
      *ptr++ = *cur_fpt++;
    buf_files[cur_buf].cur_point += sizeof(binstar_in);

    obj_lat = ((double) binstar_in.lat) / 3600000L;
    obj_lon = ((double) binstar_in.lon) / 3600000L;
    obj_mag = ((double) binstar_in.mag) / 1000L;
    obj_type[0] = 'S';
    obj_type[1] = 'S';
    obj_color[0] = ' ';
    obj_color[1] = ' ';
    obj_label[0] = ' ';
    obj_label[1] = ' ';
    obj_constell[0] = ' ';
    obj_constell[1] = ' ';
    obj_constell[2] = ' ';
    obj_name[0] = '\0';
    obj_commnt = "";

    strcpy(fileline, "");
  } else if (ftype == BINOBJ) {
    for (i = 0, ptr = (char *) &binobj_in; i < sizeof(binobj_in); i++) 
      *ptr++ = *cur_fpt++;
    buf_files[cur_buf].cur_point += sizeof(binobj_in);

    obj_lat = ((double) binobj_in.lat) / 3600000L;
    obj_lon = ((double) binobj_in.lon) / 3600000L;
    obj_mag = ((double) binobj_in.mag) / 1000L;
    obj_type[0] = binobj_in.type[0];
    obj_type[1] = binobj_in.type[1];
    obj_color[0] = ' ';
    obj_color[1] = ' ';
    obj_label[0] = ' ';
    obj_label[1] = ' ';
    obj_constell[0] = ' ';
    obj_constell[1] = ' ';
    obj_constell[2] = ' ';
    obj_name[0] = '\0';
    obj_commnt = "";

    strcpy(fileline, "");
  } else if (ftype == BINFULL) {
    for (i = 0, ptr = (char *) &binfull_in; i < sizeof(binfull_in); i++) 
      *ptr++ = *cur_fpt++;
    buf_files[cur_buf].cur_point += sizeof(binfull_in);

    if (binfull_in.strlen == 0)
      strcpy(name_comment, "");
    else {
      for (i = 0; i < binfull_in.strlen; i++) 
	name_comment[i] = *cur_fpt++;
      buf_files[cur_buf].cur_point += binfull_in.strlen;

      name_comment[binfull_in.strlen] = '\0';
    }

    obj_lat = ((double) binfull_in.lat) / 3600000L;
    obj_lon = ((double) binfull_in.lon) / 3600000L;
    obj_mag = ((double) binfull_in.mag) / 1000L;
    obj_type[0] = binfull_in.tycolb[0];
    obj_type[1] = binfull_in.tycolb[1];
    obj_color[0] = binfull_in.tycolb[2];
    obj_color[1] = binfull_in.tycolb[3];
    obj_label[0] = binfull_in.tycolb[4];
    obj_label[1] = binfull_in.tycolb[5];
    strcpy(obj_constell,con_table[binfull_in.consindx]);

    ptr = name_comment;
    i = 0;
    while (*ptr == ' ') ptr++;
    while (*ptr != ',' && *ptr != '\n' && *ptr)
      obj_name[i++] = *ptr++;
    obj_name[i] = '\0';
    if ((*ptr == ',') && (*++ptr) && name_comment[0]) obj_commnt = ptr;
    else obj_commnt = "";

    strcpy(fileline, "");
  } else { /* LINEREAD */

/*
 * file formats:
 * new
064509-1643-14SDA1a CMASirius
051432-0812015SDB8b ORIRigel
 * old
064509-1643-146SSSirius
051432-08120015SSRigel
 */

    ptr = fileline;
    while ((*ptr++ = *cur_fpt++) != '\n');
    *ptr = '\0';
    buf_files[cur_buf].cur_point = cur_fpt;

    nchars = 0;
    while (fileline[nchars++]);
    nchars--;
    nchars--;

/*
 * sscanf of floats is TOOO slow:
 *     sscanf(fileline, "%2f%2f%2f%c%2f%2f ... );
 * use alternate:
 */
#define F2(i) (((fileline[i]-'0')*10.0+fileline[i+1]-'0'))
#define F3(i) (((fileline[i]-'0')*100.0+(fileline[i+1]-'0')*10+fileline[i+2]-'0'))
#define F4(i) (((fileline[i]-'0')*1000.0+(fileline[i+1]-'0')*100+(fileline[i+2])-'0')*10+fileline[i+3]-'0')
#define F3M(i) (((fileline[i]-'A'+10.0)*100+(fileline[i+1]-'0')*10+fileline[i+2]-'0'))
    rah = F2(0);
    ram = F2(2);
    ras = F2(4);
    dld = F2(7);
    dlm = F2(9);
/*
 * common code
 */
#define DLDEGSEC 3600.0
#define DLMINSEC 60.0
#define RAHRSSEC 54000.0
#define RAMINSEC 900.0
#define RASECSEC 15.0


    obj_lon = (RAHRSSEC*rah + RAMINSEC*ram + RASECSEC*ras)/DLDEGSEC;
    dl = (DLDEGSEC*dld + DLMINSEC*dlm)/DLDEGSEC;
    obj_lat = (fileline[6]  == '-') ? -dl : dl;
    
    /* set unknowns to blanks */
    obj_color[0] = ' ';
    obj_color[1] = ' ';
    obj_color[2] = '\0';
    obj_label[0] = ' ';
    obj_label[1] = ' ';
    obj_label[2] = '\0';
    obj_constell[0] = ' ';
    obj_constell[1] = ' ';
    obj_constell[2] = ' ';
    obj_constell[3] = '\0';
    
    if  (isdigit(fileline[14])) {
    /*
     * old reduced Yale catalog
     */
      inten = F3(12);
      if (fileline[11] == '0' || fileline[11] == '+') obj_mag = inten/100.0;
      else if (fileline[11] == '-') obj_mag = -inten/100.0;
      else obj_mag = F4(11)/1000.0;	/* new feature for stars >= 10.0 mag */
      
      if (nchars > 15) {
	obj_type[0] = fileline[15];
	obj_type[1] = fileline[16];
	ptr = &fileline[MIN(17,nchars)];
	i = 0;
	while (*ptr == ' ') ptr++;
	while (*ptr != ',' && *ptr != '\n' && *ptr)
	  obj_name[i++] = *ptr++;
	obj_name[i] = '\0';
	if (*++ptr) obj_commnt = ptr;
	else obj_commnt = "";
      } else {
	obj_type[0] = obj_type[1] = 'S'; /* Default SS single star */
	obj_name[0] = '\0';
	obj_commnt = "";
      }
    } else {
      /*
       * new reduced Yale catalog
       */
      m1 = fileline[11];
      obj_mag = ((m1 == '-') ? -F2(12)/10.0 :
		 (m1 <= '9') ? F3(11)/100.0 : F3M(11)/100.0);
      /* let's get Sirius */
    
      /*
       * extract color, label, constellation, name, and comment
       * Would be faster to just guarentee that the data file is correct
       */
      if (nchars > 22) {
	obj_constell[0] = fileline[20];
	obj_constell[1] = fileline[21];
	obj_constell[2] = fileline[22];
	obj_constell[3] = '\0';
      }
      if (nchars > 19) {
	obj_label[0] = fileline[18];
	obj_label[1] = fileline[19];
	obj_label[2] = '\0';
      }
      if (nchars > 17) {
	obj_color[0] = fileline[16]; 
	obj_color[1] = fileline[17];
	obj_color[2] = '\0';
      }
      if (nchars > 15) {
	obj_type[0] = fileline[14];
	obj_type[1] = fileline[15];
      }

      ptr = &fileline[MIN(23,nchars)];
      i = 0;
      while (*ptr == ' ') ptr++;
      while (*ptr != ',' && *ptr != '\n' && *ptr)
	obj_name[i++] = *ptr++;
      obj_name[i] = '\0';
      if (*++ptr) obj_commnt = ptr;
      else obj_commnt = "";
    }
  }

  return(FALSE); /* NO error */
}
#endif /* NO_BUF_READSTAR not defined */
