#
FC = pgf90
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

w_interp:	w_interp.o
		$(FC) -o w_interp w_interp.o -L$(NCLIB) -lnetcdf 
w_interp.o:	w_interp.f90
		$(FC) -c w_interp.f90 -I$(NCINC) -L$(NCLIB) -lnetcdf 
clean:		*.o
		$(RM) *.o

