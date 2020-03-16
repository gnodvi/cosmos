#-------------------------------------------------------------------------------
#


LIBEPHACCESS = T/daf_reader.o  T/ephaccess.o 

PROGS = $(LIBEPHACCESS) ephcalculator ephaccess-checker libephaccess.so 

all: $(PROGS)

CC = LANG=C gcc -g  

#-------------------------------------------------------------------------------

ephcalculator: ephcalculator.c $(LIBEPHACCESS)
	$(CC) -o $@ $^  -lm


ephaccess-checker: ephaccess-checker.c $(LIBEPHACCESS)
	$(CC) -o $@ $^  -lm


T/daf_reader.o: daf_reader.c daf_reader.h
	$(CC) -c -o $@ $<

T/ephaccess.o: ephaccess.c
	$(CC) -c -o $@ $<



libephaccess.so:  $(LIBEPHACCESS)
	gcc  -shared -o $@ $^ 

#-------------------------------------------------------------------------------

clean:
	rm -f *.o T/*.o $(PROGS) *.pyc

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
