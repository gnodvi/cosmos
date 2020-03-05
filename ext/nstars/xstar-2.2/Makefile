#-------------------------------------------------------------------------------
#
#  Possible options:
#
#  -DNUM_COLORS=<num>
#		number of colors to use.  Default is 48, max is 64.  It
#           probably should be a multiple of 3 to give pure RGB values.
#
#  -DHASH_TABLE_BITS=<num>
#		number of bits used to store indexes into the
#		hash table.  The default is 15, the maximum (without
#		some recoding) is 16.  The larger this number, the
#		larger the hash table will be and the more point
#		locations will saved.
#
#		The default of 15 uses 64k for the hash table and
#		128k for points.
#
#  -DNEEDS_TYPEDEFS
#		define if the typedefs uint_t etc are not automatically
#		defined for you.  (HPUX 8.0 and 9.0 need this)
#
#  -DUSE_USLEEP
#		use usleep() instead of nap().  (needed for HPUX 8.0)
#
#  -DRAND48_PROTO
#		give prototypes for the *rand48() functions.  (needed
#		for HPUX 8.0)
#
#  -DUSE_RANDOM
#           Use the BSD random() functions instead of the sysv rand48()
#           functions.
#

CFLAGS = -O


# DEPLIBS = $(DEPXLIB)
LOCAL_LIBRARIES = $(XLIB)
LDLIBS = -lm -lX11

#-------------------------------------------------------------------------------

PROGS = xstar

all: $(PROGS)

CC = LANG=C gcc -g  -DNEEDS_TYPEDEFS #-DDEBUG

#-------------------------------------------------------------------------------

MOVE_NAMES	= move_stars_euler1 move_stars_taylor3 move_stars_rk4 \
	move_stars_gpemce8 move_stars_ab4 move_stars_ab7 \
	move_stars_am7 move_stars_rk4b

MOVE_SRCS = $(addsuffix   .c, ${MOVE_NAMES})
MOVE_OBJS = $(addprefix ./T/, $(addsuffix  .o, ${MOVE_NAMES}))

#-------------------------------------------------------------------------------

SYS_NAMES	= dump_sys init_sys_4 init_sys_8 init_sys_8b \
	init_sys_1a init_sys_1b init_sys_1c

SYS_SRCS = $(addsuffix   .c, ${SYS_NAMES})
SYS_OBJS = $(addprefix ./T/, $(addsuffix  .o, ${SYS_NAMES}))

#-------------------------------------------------------------------------------

MAIN_NAMES    = xstar init_stars init_colors global_vars \
	  check_bounce animate collide sys_const default_init \
	  set_sys_param set_xmva update_screen plot_collapsars \
	  reset_hstep set_star_disp

MAIN_SRCS = $(addsuffix   .c, ${MAIN_NAMES})
MAIN_OBJS = $(addprefix ./T/, $(addsuffix  .o, ${MAIN_NAMES}))

#-------------------------------------------------------------------------------

#MOVE_OBJS = $(MOVE_SRCS:.c=.o)
#SYS_OBJS  = $(SYS_SRCS:.c=.o)
#MAIN_OBJS = $(MAIN_SRCS:.c=.o)

OBJS    =  $(MAIN_OBJS) $(MOVE_OBJS) $(SYS_OBJS) 

xstar:	$(OBJS) $(DEPLIBS)
	@echo  $(MOVE_OBJS)
	$(CC) -o $@ $(OBJS) $(LDOPTIONS) $(LOCAL_LIBRARIES) $(LDLIBS) $(EXTRA_LOAD_FLAGS)


./T/%.o : %.c
	$(CC) -c -o $@ $<


INCS1 = xstar.h xstar_ext.h 

./T/animate.o::           $(INCS1)
./T/check_bounce.o::      $(INCS1)
./T/collide.o::           $(INCS1)
./T/default_init.o::      $(INCS1)
./T/dump_sys.o::          $(INCS1)
./T/global_vars.o::       $(INCS1)
./T/init_colors.o::       $(INCS1)
./T/init_stars.o::        $(INCS1)
./T/init_sys_1a.o::       $(INCS1)
./T/init_sys_1b.o::       $(INCS1)
./T/init_sys_1c.o::       $(INCS1)
./T/init_sys_4.o::        $(INCS1)
./T/init_sys_8.o::        $(INCS1)
./T/init_sys_8b.o::       $(INCS1)
./T/move_stars_ab4.o::    $(INCS1) plot_star.h
./T/move_stars_ab7.o::    $(INCS1) plot_star.h
./T/move_stars_am7.o::    $(INCS1) plot_star.h
./T/move_star~euler1.o::  $(INCS1) plot_star.h
./T/move_star~pemce8.o::  $(INCS1) plot_star.h
./T/move_stars_rk4.o::    $(INCS1) plot_star.h
./T/move_stars_rk4b.o::   $(INCS1) plot_star.h
./T/move_star~aylor3.o::  $(INCS1) plot_star.h
./T/plot_collapsars.o::   $(INCS1)
./T/reset_hstep.o::       $(INCS1)
./T/set_star_disp.o::     $(INCS1)
./T/set_sys_param.o::     $(INCS1)
./T/set_xmva.o::          $(INCS1)
./T/sys_const.o::         $(INCS1)
./T/update_screen.o::     $(INCS1)
./T/xstar.o::             $(INCS1) vroot.h patchlevel.h

#-------------------------------------------------------------------------------

#HDRS    = xstar.h patchlevel.h vroot.h
#HDREXTS = xstar_ext.h plot_star.h


#AUX	= README xstar.1 Imakefile \
#	  README.xgrav COPYING theory_of_op.ltr n-body.ps Makefile.simple

#shar:
#	shar -l 50 -o xstar.part $(AUX) $(MAIN_SRCS) $(MOVE_SRCS) $(SYS_SRCS) $(HDRS) $(HDREXTS)

#tar:
#	tar cvzf xstar.tar.gz $(AUX) $(MAIN_SRCS) $(MOVE_SRCS) $(SYS_SRCS) $(HDRS) $(HDREXTS)

#$(OBJS)::	$(HDRS)

#$(MOVE_SRCS:.c=.o)::	plot_star.h

#-------------------------------------------------------------------------------

PROGS_exe = $(addsuffix .exe, ${PROGS})

clean:
	rm -f *.o T/*.o $(PROGS)  $(PROGS_exe)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
