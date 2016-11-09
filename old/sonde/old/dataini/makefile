#
#
FC = /usr/bin/f90
RM = /usr/bin/rm
NCINC = /usr/local/netcdf/include
NCLIB = /usr/local/netcdf/lib

analyod:    	analyod.o cubicspline.o netcdfio.o regress.o
		$(FC) -o analyod analyod.o cubicspline.o regress.o  netcdfio.o -L$(NCLIB) -lnetcdf
netcdfio.o:	netcdfio.f90
		$(FC) -c -I$(NCINC) netcdfio.f90
regress.o:	regress.f90
		$(FC) -c regress.f90 
cubicspline.o:	cubicspline.f90
		$(FC) -c cubicspline.f90 cubicspline.o
analyod.o:	  	analyod.f90 netcdfio.o cubicspline.o regress.o
		$(FC) -c -I$(NCINC) analyod.f90 netcdfio.o cubicspline.o regress.o
clean:		*.o *.mod
		$(RM) *.o *.mod
