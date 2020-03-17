/* -*- linux-c -*- */
/* glstarview.h

   Copyright (C) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define VERSION "0.6"
#define NICK "Ernie"
#define DATE "Wed Apr 12 14:23:12 CDT 2006"

/* ASCII keycodes */
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

/* other stuff */
#define PI 3.14159265358979323846

#define STARSIZE 0.03f
#define NTRI 8

#define DELTAXROT 15.0f
#define DELTAYROT 15.0f
#define DELTAZ 3.0f

#define SCALE 15.0

#define ROTSCALE 1.0
#define ZOOMSCALE 0.2

/* the scene structure */
typedef struct{
	int nstar;
	GLfloat z;
	GLfloat xrot;
	GLfloat yrot;
	GLfloat starsize;
	int fullscreen;
	int drawaxes;
	int printtime;
	int printlog;
	int ntri;
	int paused;
	int step;
	int mousemode;
	int changed;
	double t;
	char tstring[128];
	char log[3500];
	double dt;
} scene_t;

/* prototypes */
void print_version(FILE *stream);
void print_usage(FILE *stream);
void print_bindings(FILE *stream);
void InitGL(int Width, int Height);
void ReSizeGLScene(int Width, int Height);
void DrawGLScene(void);
void keyPressed(unsigned char key, int x, int y);
void specialKeyPressed(int key, int x, int y);
void mousefunc(int button, int state, int x, int y);
void motionfunc(int x, int y);

/* macros */
#define MAX(a, b) ((a)>=(b)?(a):(b))
#define MIN(a, b) ((a)<=(b)?(a):(b))
#define dprintf(args...) if (debug) fprintf(stderr, args)
