#
#
FC = /usr/bin/f90
RM = /usr/bin/rm
NCINC = /usr/local/netcdf/include
NCLIB = /usr/local/netcdf/lib

energy:   	energy.o netcdfio.o
		$(FC) -o energy energy.o netcdfio.o -L$(NCLIB) -lnetcdf
netcdfio.o:	netcdfio.f90
		$(FC) -c -I$(NCINC) netcdfio.f90
energy.o:	energy.f90 netcdfio.o 
		$(FC) -c -I$(NCINC) energy.f90 netcdfio.o
clean:		*.o *.mod
		$(RM) *.o *.mod
