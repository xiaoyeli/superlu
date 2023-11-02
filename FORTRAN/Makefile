.SUFFIXES:
.SUFFIXES: .f90 .F90 .c .o .f .F
include ../make.inc

#######################################################################
#  This makefile creates the Fortran example interface to use the
#  C routines in SuperLU.
#######################################################################

HEADER   = ../SRC
LIBS	+= -lm

# double real
DFEXMOBJ = f_main.o hbcode1.o c_fortran_dgssv.o
TESTOMPOBJ= test_omp.o hbcode1.o c_fortran_dgssv.o

# double complex
ZFEXMOBJ = z_f_main.o zhbcode1.o c_fortran_zgssv.o

all:	dfexm zfexm test_omp

superlu_config.fh: ../SRC/superlu_config.h
	sed 's/\/.*$$//' < $< > temp.fh
	sed '/typedef/ d' < temp.fh > $@
	rm temp.fh

dfexm: $(DFEXMOBJ) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(DFEXMOBJ) $(LIBS) -o $@

zfexm: $(ZFEXMOBJ) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(ZFEXMOBJ) $(LIBS) -o $@

test_omp: $(TESTOMPOBJ) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(TESTOMPOBJ)  $(LIBS) -o $@

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f90.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

.F90.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

.F.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o *fexm test_omp

