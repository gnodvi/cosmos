/* -*- linux-c -*- */
/* glstarview.c
   
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
/*
  This code was created by Jeff Molofee '99 (ported to Linux/GLUT by Richard Campbell '99)

  If you've found this code useful, please let me know.

  Visit me at www.demonews.com/hosted/nehe 
  (email Richard Campbell at ulmont@bellsouth.net)
*/

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "glstarview.h"

/* So many global variables... */
int window, xlast=0, ylast=0, debug=0, width=640, height=480;
float r[1024][3], color[1024][3];
GLUquadricObj *quadratic;
scene_t scene;
GLuint base;

/* print the version */
void print_version(FILE *stream)
{
	fprintf(stream, "** GLStarView %s (%s) [%s] **\n", VERSION, NICK, DATE);
}

/* print the usage */
void print_usage(FILE *stream)
{
	fprintf(stream, "USAGE:\n");
	fprintf(stream, "  glstarview [options...]\n");
	fprintf(stream, "\n");
	fprintf(stream, "OPTIONS:\n");
	fprintf(stream, "  -f --fullscreen : full screen mode\n");
	fprintf(stream, "  -p --paused     : start paused\n");
	fprintf(stream, "  -a --noaxes     : start with axes off\n");
	fprintf(stream, "  -t --time       : start with time display on\n");
	fprintf(stream, "  -l --log        : start with log display on\n");
	fprintf(stream, "  -d --debug      : turn on debugging\n");
	fprintf(stream, "  -V --version    : print version info\n");
	fprintf(stream, "  -h --help       : display this help text\n");
	fprintf(stream, "\n");
	print_bindings(stream);
}

/* print the bindings */
void print_bindings(FILE *stream)
{
	fprintf(stream, "KEYBINDINGS:\n");
	fprintf(stream, "  PAGE (UP/DOWN) : zoom (out/in)\n");
	fprintf(stream, "  ARROW KEYS     : rotate\n");
	fprintf(stream, "  h              : print this help text\n");
	fprintf(stream, "  q, ESCAPE      : quit\n");
	fprintf(stream, "  p              : pause/unpause\n");
	fprintf(stream, "  SPACE          : step through frames (while paused)\n");
	fprintf(stream, "  b              : make stars bigger\n");
	fprintf(stream, "  s              : make stars smaller\n");
	fprintf(stream, "  m              : use more polygons to render stars\n");
	fprintf(stream, "  f              : use fewer polygons to render stars\n");
	fprintf(stream, "  a              : toggle axes\n");
	fprintf(stream, "  t              : toggle time printout\n");
	fprintf(stream, "  l              : toggle log printout\n");
	fprintf(stream, "  <              : slow down animation\n");
	fprintf(stream, "  >              : speed up animation\n");
	fprintf(stream, "\n");
	fprintf(stream, "MOUSEBINDINGS:\n");
	fprintf(stream, "  DRAG                   : rotate\n");
	fprintf(stream, "  [SHIFT] DRAG (up/down) : zoom (in/out)\n");
}

/* A general OpenGL initialization function.  Sets all of the initial parameters. */
void InitGL(int Width, int Height)
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	
	glClearDepth(1.0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	quadratic=gluNewQuadric();
	gluQuadricNormals(quadratic, GLU_SMOOTH);
	gluQuadricTexture(quadratic, GL_TRUE);
	
	gluPerspective(45.0f, (GLfloat)Width/(GLfloat)Height, 0.1f, 1000.0f);
	
	glMatrixMode(GL_MODELVIEW);
}

/* The function called when our window is resized */
void ReSizeGLScene(int Width, int Height)
{
	/* Prevent A Divide By Zero If The Window Is Too Small */
	if (Height==0) {
		Height=1;
	}
	
	width = Width;
	height = Height;

	glViewport(0, 0, Width, Height);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	gluPerspective(45.0f, (GLfloat)Width/(GLfloat)Height, 0.1f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);

	scene.changed = 1;
}

/* The main drawing function */
void DrawGLScene(void)
{
	int i=0;
	static char *str, *ret="blah", line[16384]="";
	double usleepdt=100.0;
	static double usleepcount=0.0;
	
	usleepcount += usleepdt;

	/* read data if not paused */
	if (usleepcount >= scene.dt && (scene.paused == 0 || scene.step == 1)) {
		usleepcount = 0.0;
		scene.step = 0;
		
		/* parse Starlab format */
		while ((ret = fgets(line, 16384, stdin)) != NULL && !strstr(line, "(P")) {
		}
		
		if (ret != NULL) {
			scene.nstar = 0;
			scene.changed = 1;
			
			while ((ret = fgets(line, 16384, stdin)) != NULL && strcmp(line, ")Particle\n")) {
				str = strtok(line, " =");
				
				/* get log */
				if (strcmp(line, "(Log\n") == 0) {
					i = 0;
					while ((ret = fgets(line, 16384, stdin)) != NULL && strcmp(line, ")Log\n")) {
						i += snprintf(&(scene.log[i]), (3499-i)<0?0:(3499-i), "%s", line);
					}
					if (i != 0) {
						dprintf("%s", scene.log);
					}
				}

				/* get time */
				if (strcmp(line, "(Dynamics\n") == 0) {
					while ((ret = fgets(line, 16384, stdin)) != NULL && strcmp(line, ")Dynamics\n")) {
						str = strtok(line, " =");
						if (strcmp(str, "system_time") == 0) {
							str = strtok(NULL, " =");
							scene.t = atof(str);
							snprintf(scene.tstring, 127, "t=%.4g", scene.t);
							dprintf("%s\n", scene.tstring);
						}
					}
				}
				
				/* get particle positions */
				if (strcmp(line, "(Particle\n") == 0) {
					scene.nstar++;
					while ((ret = fgets(line, 16384, stdin)) != NULL && strcmp(line, ")Particle\n")) {
						str = strtok(line, " =");
						if (strcmp(str, "(Dynamics\n")==0) {
							while ((ret = fgets(line, 16384, stdin)) != NULL && strcmp(line, ")Dynamics\n")) {
								str = strtok(line, " =");
								if (strcmp(str, "r") == 0) {
									str = strtok(NULL, " =");
									r[i-1][0] = atof(str);
									str = strtok(NULL, " =");
									r[i-1][1] = atof(str);
									str = strtok(NULL, " =");
									r[i-1][2] = atof(str);
								}
							}
						}
						
						if (strcmp(str, "i") == 0) {
							str = strtok(NULL, " =");
							i = atoi(str);
						}
					}
				}
			}
		}
	}
	
	usleep(usleepdt);

	/* draw scene if necessary */
	if (scene.changed) {
		scene.changed = 0;
		
		/* erase scene */
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		glLoadIdentity();
		
		/* print text */
		if (scene.printtime) {
			fprintf(stdout, "glstarview: t=%s\n", scene.tstring);
		}
		if (scene.printlog) {
			fprintf(stdout, "glstarview: %s", scene.log);
		}

		/* get ready to render 3D objects */
		glTranslatef(0.0f, 0.0f, scene.z+0.1);
		
		/* rotate view--I know these two rotation matrices don't commute */
		glRotatef(scene.xrot, 1.0f, 0.0f, 0.0f);
		glRotatef(scene.yrot, 0.0f, 1.0f, 0.0f);
		
		/* draw axes */
		if (scene.drawaxes) {
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(1.0, 0.0, 0.0);
			glEnd();
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 1.0, 0.0);
			glEnd();
			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_LINES);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 1.0);
			glEnd();
		}
		
		/* draw stars */
		for (i=0; i<scene.nstar; i++) {
			glLoadIdentity();
			glTranslatef(0.0f, 0.0f, scene.z);
			glRotatef(scene.xrot, 1.0f, 0.0f, 0.0f);
			glRotatef(scene.yrot, 0.0f, 1.0f, 0.0f);
			glTranslatef(r[i][0]/SCALE, r[i][1]/SCALE, r[i][2]/SCALE);
			glColor3f(color[i][0], color[i][1], color[i][2]);
			/* glutSolidSphere(scene.starsize, scene.ntri, scene.ntri); */
			gluSphere(quadratic, scene.starsize, scene.ntri, scene.ntri);
			/* glTranslatef(-r[i][0]/SCALE, -r[i][1]/SCALE, -r[i][2]/SCALE); */
		}
		
		glutSwapBuffers();
	}
}

/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y) 
{
	/* avoid thrashing this call */
	usleep(100);
	
	if (key == ESCAPE || key == 'q') {
		glutDestroyWindow(window);
		exit(0);
	} else if (key == 'h') {
		print_bindings(stderr);
	} else if (key == 'p') {
		if (scene.paused == 0) {
			scene.paused = 1;
		} else {
			scene.paused = 0;
		}
	} else if (key == ' ') {
		scene.step = 1;
	} else if (key == 'b') {
		scene.starsize *= 1.5;
		if (scene.starsize > 100.0) {
			scene.starsize = 100.0;
		}
		scene.changed = 1;
	} else if (key == 's') {
		scene.starsize /= 1.5;
		if (scene.starsize < 1.0e-3) {
			scene.starsize = 1.0e-3;
		}
		scene.changed = 1;
	} else if (key == 'm') {
		scene.ntri *= 2;
		if (scene.ntri > 24) {
			scene.ntri = 24;
		}
		scene.changed = 1;
	} else if (key == 'f') {
		scene.ntri /= 2;
		if (scene.ntri < 4) {
			scene.ntri = 4;
		}
		scene.changed = 1;
	} else if (key == 'a') {
		if (scene.drawaxes == 0) {
			scene.drawaxes = 1;
		} else {
			scene.drawaxes = 0;
		}
		scene.changed = 1;
	} else if (key == 't') {
		if (scene.printtime == 0) {
			scene.printtime = 1;
		} else {
			scene.printtime = 0;
		}
	} else if (key == 'l') {
		if (scene.printlog == 0) {
			scene.printlog = 1;
		} else {
			scene.printlog = 0;
		}
	} else if (key == '>') {
		scene.dt /= 2.0;
		if (scene.dt < 10.0) {
			scene.dt = 10.0;
		}
	} else if (key == '<') {
		scene.dt *= 2.0;
		if (scene.dt > 1.0e6) {
			scene.dt = 1.0e6;
		}
	}
}

/* The function called whenever a normal key is pressed. */
void specialKeyPressed(int key, int x, int y) 
{
	/* avoid thrashing this procedure */
	usleep(100);
	
	switch (key) {    
	case GLUT_KEY_PAGE_UP:
		scene.z-=DELTAZ;
		scene.changed = 1;
		break;
	case GLUT_KEY_PAGE_DOWN:
		scene.z+=DELTAZ;
		scene.changed = 1;
		break;
	case GLUT_KEY_UP:
		scene.xrot -= DELTAXROT;
		scene.changed = 1;
		break;
	case GLUT_KEY_DOWN:
		scene.xrot += DELTAXROT;
		scene.changed = 1;
		break;
	case GLUT_KEY_LEFT:
		scene.yrot -= DELTAYROT;
		scene.changed = 1;
		break;
	case GLUT_KEY_RIGHT:
		scene.yrot += DELTAYROT;
		scene.changed = 1;
		break;
	default:
		break;
	}	
}

void mousefunc(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN) {
		xlast = x;
		ylast = y;
		if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
			scene.mousemode = 2;
		} else {
			scene.mousemode = 1;
		}
	}
}

void motionfunc(int x, int y)
{
	double dx, dy;
	
	dx = ((double) x) - ((double) xlast);
	dy = ((double) y) - ((double) ylast);
	
	xlast = x;
	ylast = y;
	
	if (scene.mousemode == 1) {
		scene.xrot += ROTSCALE * dy;
		scene.yrot += ROTSCALE * dx;
		scene.changed = 1;
	} else if (scene.mousemode == 2) {
		scene.z -= ZOOMSCALE * dy;
		scene.changed = 1;
	}
}

int main(int argc, char **argv) 
{
	int i;
	double R, G, B, scale;
	const char *short_opts = "fpatldVh";
	const struct option long_opts[] = {
		{"fullscreen", no_argument, NULL, 'f'},
		{"paused", no_argument, NULL, 'p'},
		{"noaxes", no_argument, NULL, 'a'},
		{"time", no_argument, NULL, 't'},
		{"log", no_argument, NULL, 'l'},
		{"debug", no_argument, NULL, 'd'},
		{"version", no_argument, NULL, 'V'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}
	};
	
	/* set a few initial variables */
	scene.nstar = 0;
	scene.xrot = 0.0f;
	scene.yrot = 0.0f;
	scene.z = -6.0f;
	scene.starsize = STARSIZE;
	scene.fullscreen = 0;
	scene.drawaxes = 1;
	scene.printtime = 0;
	scene.printlog = 0;
	scene.ntri=NTRI;
	scene.paused = 0;
	scene.step = 0;
	scene.mousemode = 1;
	scene.changed = 0;
	scene.dt = 300;

	while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (i) {
		case 'f':
			scene.fullscreen = 1;
			break;
		case 'p':
			scene.paused = 1;
			break;
		case 'a':
			scene.drawaxes = 0;
			break;
		case 't':
			scene.printtime = 1;
			break;
		case 'l':
			scene.printlog = 1;
			break;
		case 'd':
			debug = 1;
			break;
		case 'V':
			print_version(stdout);
			return(0);
		case 'h':
			print_version(stdout);
			fprintf(stdout, "\n");
			print_usage(stdout);
			return(0);
		default:
			break;
		}
	}
	
	/* check to make sure there was nothing crazy on the command line */
	if (optind < argc) {
		print_usage(stdout);
		return(1);
	}

	/* assign particle colors */
	srand(823);

	color[0][0] = 1.0;
	color[0][1] = 0.0;
	color[0][2] = 0.0;

	color[1][0] = 0.0;
	color[1][1] = 1.0;
	color[1][2] = 0.0;

	color[2][0] = 0.0;
	color[2][1] = 0.0;
	color[2][2] = 1.0;

	color[3][0] = 1.0;
	color[3][1] = 1.0;
	color[3][2] = 0.0;

	color[4][0] = 1.0;
	color[4][1] = 0.0;
	color[4][2] = 1.0;

	color[5][0] = 0.0;
	color[5][1] = 1.0;
	color[5][2] = 1.0;
	
	for (i=6; i<1024; i++) {
		R = ((double) rand())/((double) RAND_MAX);
		G = ((double) rand())/((double) RAND_MAX) * (1.0 - R);
		B = 1.0 - R - G;
		if (R >= G && R >= B) {
			scale = 1.0 + ((double) rand())/((double) RAND_MAX) * (MIN(2.0,1.0/R) - 1.0);
		} else if (G >= R && G >= B) {
			scale = 1.0 + ((double) rand())/((double) RAND_MAX) * (MIN(2.0,1.0/G) - 1.0);
		} else {
			scale = 1.0 + ((double) rand())/((double) RAND_MAX) * (MIN(2.0,1.0/B) - 1.0);
		}
		color[i][0] = R * scale;
		color[i][1] = G * scale;
		color[i][2] = B * scale;
	}

	/* prepare rendering */
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(0, 0);
	window = glutCreateWindow("glstarview");
	glutDisplayFunc(&DrawGLScene);
	
	if (scene.fullscreen) {
		glutFullScreen();
	}

	glutIdleFunc(&DrawGLScene);
	glutReshapeFunc(&ReSizeGLScene);
	glutKeyboardFunc(&keyPressed);
	glutSpecialFunc(&specialKeyPressed);
	glutMouseFunc(&mousefunc);
	glutMotionFunc(&motionfunc);
	
	InitGL(width, height);
	
	glutMainLoop();
	
	/* should never get here */
	return(1);
}
