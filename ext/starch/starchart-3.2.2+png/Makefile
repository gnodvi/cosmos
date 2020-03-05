#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

#PROGS = libabmd.so main 

PROGS = planet  starpost starpng stardsp

all: $(PROGS)

#-------------------------------------------------------------------------------

stardir = .

STAR   = \"${stardir}/D/yale.star\"
NEB    = \"${stardir}/D/messier.star\"
BOUND  = \"${stardir}/D/boundaries.star\"
PATT   = \"${stardir}/D/pattern.star\"
CNAME  = \"${stardir}/D/cnames.star\"
PLANET = \"./planet.star\"
CONS   = \"${stardir}/D/con.locs\"


FILEFLAGS = -DSTARFILE=$(STAR) -DNEBFILE=$(NEB) -DBOUNDFILE=$(BOUND) -DPATTERNFILE=$(PATT) -DCNAMEFILE=$(CNAME) -DCONSTFILE=$(CONS) -DPLANETFILE=$(PLANET) 


CC = LANG=C gcc $(FILEFLAGS)  

#-------------------------------------------------------------------------------


starpng_OBJECTS =  starmain.o parse_input.o readfile.o starm2.o \
starsupp.o starcust.o starpng.o

stardsp_OBJECTS =  starmain.o parse_input.o readfile.o starm2.o \
starsupp.o starcust.o stardsp.o

planet_OBJECTS =  planet.o

#-------------------------------------------------------------------------------

starpost:  starmain.o parse_input.o readfile.o starm2.o starsupp.o starcust.o starpost.o 
	$(CC) -o $@ $^  -lm -lgd -lz 

starpng:  $(starpng_OBJECTS) 
	$(CC) -o $@ $^  -lm -lgd -lz 

stardsp:  $(stardsp_OBJECTS) 
	$(CC) -o $@ $^  -lm -lgd -lz 


planet: planet.o
	$(CC) -o $@ $^  -lm -lgd -lz 

#-------------------------------------------------------------------------------

starmain.o: starmain.c  
	$(CC) -c -o $@ $<  -I. 

parse_input.o: parse_input.c  
	$(CC) -c -o $@ $<  -I. 

readfile.o: readfile.c  
	$(CC) -c -o $@ $<  -I. 

starm2.o: starm2.c  
	$(CC) -c -o $@ $<  -I. 

.starsuppo: starsupp.c  
	$(CC) -c -o $@ $<  -I. 

starcust.o: starcust.c  
	$(CC) -c -o $@ $<  -I. 

starpost.o: starpost.c  
	$(CC) -c -o $@ $<  -I. 

#-------------------------------------------------------------------------------

planet.o: planet.c  
	$(CC) -c -o $@ $<  -I. 

#-------------------------------------------------------------------------------

clean:
	rm -f *.o T/*.o $(PROGS) 

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
