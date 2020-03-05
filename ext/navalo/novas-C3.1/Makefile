#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


ALLOBJS = T/novas.o T/novascon.o T/nutation.o T/readeph0.o T/solsys1.o \
	T/solsys2.o T/solsys3.o T/example.o T/checkout-mp.o \
	T/checkout-stars-full.o T/checkout-stars.o T/cio_file.o T/eph_manager.o

PROGS = $(ALLOBJS) example 


all: $(PROGS)

CC = LANG=C gcc -g  

#-------------------------------------------------------------------------------

LIB = T/novas.o T/eph_manager.o T/novascon.o T/nutation.o T/solsys1.o T/readeph0.o

example: example.c $(LIB)
	$(CC) -o $@ $< $(LIB)  -lm -I.


#-------------------------------------------------------------------------------

T/novas.o: novas.c novas.h
	$(CC) -c -o $@ $<

T/novascon.o: novascon.c novascon.h
	$(CC) -c -o $@ $<

T/nutation.o: nutation.c nutation.h
	$(CC) -c -o $@ $<

T/readeph0.o: readeph0.c 
	$(CC) -c -o $@ $<

T/solsys1.o: solsys1.c 
	$(CC) -c -o $@ $<

T/solsys2.o: solsys2.c 
	$(CC) -c -o $@ $<

T/solsys3.o: solsys3.c 
	$(CC) -c -o $@ $<

T/example.o: example.c 
	$(CC) -c -o $@ $<

T/checkout-mp.o: checkout-mp.c 
	$(CC) -c -o $@ $<

T/checkout-stars-full.o: checkout-stars-full.c 
	$(CC) -c -o $@ $<

T/checkout-stars.o: checkout-stars.c 
	$(CC) -c -o $@ $<

T/cio_file.o: cio_file.c 
	$(CC) -c -o $@ $<

T/eph_manager.o: eph_manager.c eph_manager.h
	$(CC) -c -o $@ $<


#-------------------------------------------------------------------------------

# libabmd.so: T/l_abmd.o T/l_api.o T/l_queue.o T/l_rk.o 
# 	gcc  -shared -o $@ $^ 


#-------------------------------------------------------------------------------

clean:
	rm -f *.o T/*.o $(PROGS) 

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
