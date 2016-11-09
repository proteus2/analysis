#
#
FC = pgf90
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib


# exec

runn:		anal.o  netio.o um_anal.o um_axis.o pwrspd.o fft.o regress.o
		$(FC) -o runn anal.o netio.o -L$(NCLIB) -lnetcdf um_anal.o um_axis.o pwrspd.o fft.o regress.o

# main

anal.o:	anal.f90  netio.o um_anal.o um_axis.o pwrspd.o fft.o pwrspd.o
		$(FC) -c -I$(NCINC) anal.f90  netio.o um_anal.o um_axis.o pwrspd.o fft.o regress.o

# modules

netio.o:	netio.f90
		$(FC) -c -I$(NCINC) netio.f90

um_anal.o:	um_anal.f90
		$(FC) -c um_anal.f90

um_axis.o:	um_axis.f90 netio.o
		$(FC) -c -I$(NCINC) um_axis.f90 netio.o

regress.o:	regress.f90
		$(FC) -c regress.f90

fft.o:	fft.f90
		$(FC) -c fft.f90

pwrspd.o:	pwrspd.f90 fft.o regress.o
		$(FC) -c pwrspd.f90 fft.o regress.o

# clean

clean:	*.o *.mod
		$(RM) *.o *.mod
