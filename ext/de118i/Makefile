#-------------------------------------------------------------------------------
# Unix make file for ssystem
# Be sure to define machine in mconf.h
# and arithmetic precision in prec.h.
#
# If your system supports 80-bit long double arithmetic
# then add the following subroutines: sinl.c, atanl.c, tanl.c,
# asinl.c, polevll.c, floorl.c. Otherwise, leave them out.
#
#-------------------------------------------------------------------------------

PROGS = ssystem # ssystem_new

all: $(PROGS)

#-------------------------------------------------------------------------------

# Name your compiler.
#CC=gcc
CC = LANG=C gcc 

# Optimizations:
CFLAGS= -O2

# Subroutine libraries:
LIBS= -lc -lm
#-lc881u -lm881u

#-------------------------------------------------------------------------------

OBJS= T/adams4.o T/runge.o \
T/oblate.o T/precess.o T/epsiln.o T/nut1t.o T/zatan2.o \
T/rdnums.o T/ieee.o T/econst.o T/mtherr.o T/findcent.o \
T/reltiv.o T/jplmp.o \
T/asinl.o T/atanl.o T/sinl.o T/tanl.o T/polevll.o  T/aconst.o   # T/sqrtl.o T/floorl.o 

#-------------------------------------------------------------------------------
# For real Unix:

ssystem: T/ssystem.o $(OBJS)
	$(CC) $(CFLAGS) -o $@  $^ $(LIBS)

#ssystem_new: T/ssystem_new.o $(OBJS)
#	$(CC) $(CFLAGS) -o $@  $^ $(LIBS)


# For djgcc on MSDOS:
#ssystem: $(OBJS)
#	>ssystem.rf $(CFLAGS) $(OBJS) $(LIBS) -o ssystem
#	$(CC) @ssystem.rf
#	aout2exe ssystem

#-------------------------------------------------------------------------------

./T/%.o : %.c
	$(CC) -c -o $@ $<


T/ssystem.o: ssystem.c  ssystem.h ini118d.h mconf.h prec.h
	$(CC) -c -o $@ $<

T/ssystem_new.o: ssystem.c  ssystem.h ini118d.h mconf.h prec.h
	$(CC) -c -o $@ $< -D_NEW

#-------------------------------------------------------------------------------

T/adams4.o: adams4.c int.h prec.h mconf.h

T/runge.o: runge.c int.h prec.h mconf.h

T/aconst.o: aconst.c prec.h mconf.h ini118d.h

T/oblate.o: oblate.c ssystem.h mconf.h prec.h

T/precess.o: precess.c mconf.h prec.h

T/epsiln.o: epsiln.c mconf.h prec.h

T/nut1t.o: nut1t.c mconf.h prec.h

T/zatan2.o: zatan2.c

T/rdnums.o: rdnums.c

T/ieee.o: ieee.c ehead.h

T/econst.o: econst.c ehead.h

T/mtherr.o: mtherr.c mconf.h

T/reltiv.o: reltiv.c

T/findcent.o: findcent.c

T/jplmp.o: jplmp.c

T/asinl.o: asinl.c mconf.h

T/atanl.o: atanl.c mconf.h

T/sinl.o: sinl.c mconf.h

T/tanl.o: tanl.c mconf.h

T/polevll.o: polevll.c mconf.h

#-------------------------------------------------------------------------------
# IBM PC version
#sqrtl.o: sqrtl.387
#	as -o sqrtl.o sqrtl.387
#floorl.o: floorl.387
#	as -o floorl.o floorl.387

# Use this for x86 linux
#floorl.o: floorlelf.387
#	as -o floorl.o floorlelf.387
#sqrtl.o: sqrtlelf.387
#	as -o sqrtl.o sqrtlelf.387

# 68k version
#sqrtl.o: sqrtl.sun
#	as -o sqrtl.o sqrtl.sun
#floorl.o: floorl.sun
#	as -o floorl.o floorl.sun
#
# The C version of sqrtl does not give a strictly rounded result.
#sqrtl.o: sqrtl.c mconf.h
#
#floorl.o: floorl.c mconf.h
#
#-------------------------------------------------------------------------------

clean:
	rm -f *.o T/*.o $(PROGS)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
