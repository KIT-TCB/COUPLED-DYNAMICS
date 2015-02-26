#CC= icc
#CC= gcc
CCOPTS=  -I../include -I../src/kernel -I. 
#LDOPTS= -lm
#LDOPTS= -lblas -llapack
#LDOPTS= -lmk_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack -lpthread
#LDOPTS= -llapack -lblas -lg2c
OBJ= broyden.o fermi.o eglcao.o eglcao2.o gammamat.o mulliken.o orthogonalize.o output.o skpar.o slkode.o slktrafo.o gradient.o

libdftb1: $(OBJ)
	#$(CC) $(LDOPTS) $(OBJ)
	ar rvs libdftb.a $(OBJ)
.c.o:
	$(CC) $(CCOPTS) -Wall -c -g $<
clean:
	rm -f libdftb.a $(OBJ)
