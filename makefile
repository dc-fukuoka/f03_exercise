all:
	ifort -c $(MKLROOT)/include/mkl_vsl.f90
	ifort -g -O3 -fopenmp test.F90 -mkl

clean:
	rm -f a.out *.o *.mod *~
