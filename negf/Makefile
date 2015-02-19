FC= gfortran 
FCOPTS= -ffree-line-length-999 
#FCOPTS=  -I../gromacs-4.0.3/include -I../gromacs-4.0.3/src/kernel
OBJ= negf_lorentz.o rksuite.o

libnegf: $(OBJ)
	ar rvs libnegf.a $(OBJ)
negf_lorentz.o: negf_lorentz.f90
	$(FC) $(FCOPTS) -Wall -c -g negf_lorentz.f90
rksuite.o: rksuite.f
	$(FC) $(FCOPTS) -Wall -c -g rksuite.f
clean:
	rm -f libnegf.a $(OBJ)

