#
# Makefile para montar el programa DFT4He3d 
#  (Density Functional Theory 4He 3dimensional program)
#
# (Version de compilacion para Pentium IV)
#
COMP = ifort
CFLAGS = -c  -static -module ./modules
CFLAGS = -c -O3 -tpp6 -axK -xK -tune pn3 -arch pn3 -static -reentrancy threaded -threads -module ./modules
CFLAGS = -c -O0 -tpp7 -static -reentrancy threaded -threads -module ./modules
CFLAGS = -c -g -check -static -module ./modules 
CFLAGS = -c -O3 -tpp7 -unroll2048 -reentrancy threaded -threads -module ./modules
CFLAGS = -c -O3 -tpp7 -unroll8192 -align -axN -xN -arch pn4 -static -module ./modules
CFLAGS = -c -O3 -tpp7  -arch pn4  -module ./modules
CFLAGS = -c -fast -fpp -fopenmp -parallel -qopt-matmul -unroll -module ./modules -xAVX -axAVX -arch AVX -mkl
#CFLAGS = -c -fast -fpp -fopenmp -parallel -opt-matmul -unroll -module ./modules -xAVX -axAVX -arch AVX
#CFLAGS = -c -g  -check -static -module ./modules 

#LD_FLAGS = -lguide_stats -lguide -lcprts -lcxa -lcxaguard \
#-lifcoremt -lunwind -threads
#LD_FLAGS = -threads -parallel -opt-matmul
LD_FLAGS = -threads -parallel -qopt-matmul -mkl

#   Libraries for the FFT
LIB1=fftw3
LIB4=m
LIB2=fftw3_threads
LIB3=pthread
#LIB5=svml


#   Name of the program
PROGNAME=vectorimp_absor


#   Fortran objects
objs=init_deriv_parallel.o	modules.o	V_impur.o	DFT4He3d.o	derden.o	dimen.o\
		 energy.o	fforma.o	fft.o		initcg.o	morse.o\
		mates.o		poten.o		printoutc.o	r_cm.o		instates.o	diag.o\
		readenc.o	respar.o	term_alfa.o	timer.o		titols.o\
		vareps.o	varmu.o		tstgrid.o	s13adf.o	pdergc.o\
		steprk.o	steppc.o	potenimp.o	potenimpini.o	newder.o	redef.o
#
.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(objs)
	$(COMP)	-o $(PROGNAME) $(objs)  $(LD_FLAGS)
#	$(COMP)	-o $(PROGNAME) $(objs) -l$(LIB4) -l$(LIB1) -l$(LIB2) -l$(LIB3)  $(LD_FLAGS)
#	$(COMP)	-o $(PROGNAME) $(objs) -l$(LIB5) -l$(LIB4) -l$(LIB1) -l$(LIB2) -l$(LIB3)  $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
modules:	
	$(COMP) $(CFLAGS) $(mods)

clear:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
clean:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
