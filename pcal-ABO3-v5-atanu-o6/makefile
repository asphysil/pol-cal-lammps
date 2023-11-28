#!/usr/bin/bash

FC=gfortran
FCFLAG= -Wall -Wextra -pedantic -fimplicit-none -fbacktrace 
LIBPATH=
LIB= ${LIBPATH} #-llapack -lblas
objects=dtype_constant.o dtype_struc.o\
        neighbour_list.o\
        relative2O_disp.o\
	readfiles_list.o\
	writefiles_list.o\
	polarizationCal_md_v2.o  

exe.x:$(objects) 
	$(FC) ${FCFLAG} -o exeO6O12.x $(objects) $(LIB) 

dtype_constant.o:dtype_constant.f90 
	$(FC) ${FCFLAG} -c dtype_constant.f90 $(LIB)
dtype_struc.o:dtype_struc.f90 
	$(FC) ${FCFLAG} -c dtype_struc.f90 $(LIB)
relative2O_disp.o: relative2O_disp.f90
	$(FC) ${FCFLAG} -c relative2O_disp.f90 $(LIB)	
neighbour_list.o:neighbour_list.f90 
	$(FC) ${FCFLAG} -c neighbour_list.f90 $(LIB)
readfiles_list.o:readfiles_list.f90 
	$(FC) ${FCFLAG} -c readfiles_list.f90 $(LIB)
writefiles_list.o:writefiles_list.f90 
	$(FC) ${FCFLAG} -c writefiles_list.f90 $(LIB)
polarizationCal_md_v2.o:polarizationCal_md_v2.f90
	$(FC) ${FCFLAG} -c polarizationCal_md_v2.f90 $(LIB) 
clean:
	rm $(objects) *.mod  exeO6O12.x
