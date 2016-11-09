#
#
FC = /usr/bin/f90
RM = /usr/bin/rm
NCINC = /usr/local/netcdf/include
NCLIB = /usr/local/netcdf/lib

runnn:	anal.o netcdfio.o pwrspd.o fft.o regress.o
		$(FC) -o runnn anal.o netcdfio.o -L$(NCLIB) -lnetcdf pwrspd.o fft.o regress.o
netcdfio.o:	netcdfio.f90
		$(FC) -c -I$(NCINC) netcdfio.f90
regress.o:	regress.f90
		$(FC) -c regress.f90
fft.o:	fft.f90
		$(FC) -c fft.f90
pwrspd.o:	pwrspd.f90 fft.o regress.o
		$(FC) -c pwrspd.f90  fft.o regress.o
anal.o:	anal.f90 netcdfio.o pwrspd.o fft.o regress.o
		$(FC) -c -I$(NCINC) anal.f90  netcdfio.o pwrspd.o fft.o regress.o
clean:	*.o *.mod
		$(RM) *.o *.mod
