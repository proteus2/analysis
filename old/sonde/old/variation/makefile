#
#
FC = /usr/bin/f90
RM = /usr/bin/rm
NCINC = /usr/local/netcdf/include
NCLIB = /usr/local/netcdf/lib

variation:	variation.o netcdfio.o
		$(FC) -o variation variation.o netcdfio.o -L$(NCLIB) -lnetcdf
netcdfio.o:	netcdfio.f90
		$(FC) -c -I$(NCINC) netcdfio.f90
variation.o:	variation.f90 netcdfio.o 
		$(FC) -c -I$(NCINC) variation.f90 netcdfio.o
clean:		*.o *.mod
		$(RM) *.o *.mod
