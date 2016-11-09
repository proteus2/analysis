#
#
FC = /usr/bin/f90
RM = /usr/bin/rm
NCINC = /usr/local/netcdf/include
NCLIB = /usr/local/netcdf/lib

od_addstbas:   	od_addstbas.o netcdfio.o
		$(FC) -fpe3 -O2 -o od_addstbas od_addstbas.o netcdfio.o -L$(NCLIB) -lnetcdf
netcdfio.o:	netcdfio.f90
		$(FC) -c -I$(NCINC) netcdfio.f90
od_addstbas.o:	od_addstbas.f90 netcdfio.o 
		$(FC) -fpe3 -O2 -c -I$(NCINC) od_addstbas.f90 netcdfio.o
clean:		*.o *.mod
		$(RM) *.o *.mod
