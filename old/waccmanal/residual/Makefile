FC = pgf90
RM = /bin/rm
NCLIB = $(NETCDF)/lib
NCINC = $(NETCDF)/include
OPT = #-omp

residual:	residual.o interpolation.o mud2crf90.o mudcom.o smth9.o
		$(FC) $(OPT) -o residual residual.o interpolation.o \
		mud2crf90.o mudcom.o smth9.o -L$(NCLIB) -lnetcdf 
smth9.o:	smth9.F90
		$(FC) $(OPT) -c smth9.F90
interpolation.o:	interpolation.F90
		$(FC) $(OPT) -c interpolation.F90
mud2crf90.o:	mud2crf90.F mudcom.o
		$(FC) $(OPT) -c mud2crf90.F
mudcom.o:	mudcom.F
		$(FC) $(OPT) -c mudcom.F
residual.o:	residual.F90 interpolation.o mud2crf90.o smth9.o
		$(FC) $(OPT) -c residual.F90 -I$(NCINC)
clean:		*.mod *.o
		$(RM) *.mod *.o
