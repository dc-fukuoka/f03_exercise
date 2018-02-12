all:
	ifort -c $(MKLROOT)/include/mkl_vsl.f90
	ifort -warn all -g -O3 -mavx -fopenmp test.F90 -mkl
#	gfortran -g -O3 -mavx -fno-range-check  -c $(MKLROOT)/include/mkl_vsl.f90
#	gfortran -g -O3 -mavx -fopenmp test.F90 -lmkl_gnu_thread -lmkl_rt -lmkl_core -lmkl_intel_lp64

clean:
	rm -f a.out *.o *.mod *~
