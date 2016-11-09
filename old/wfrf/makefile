#
#
FC = pgf90
RM = rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

cal_wfrf:	cal_wfrf.o cal_xsq.o params.o quadpack.o 
		$(FC) -o cal_wfrf cal_wfrf.o cal_xsq.o params.o quadpack.o -L$(NCLIB) -lnetcdf
cal_wfrf.o:	cal_wfrf.f90 cal_xsq.o params.o quadpack.o
		$(FC) -c -I$(NCINC) cal_wfrf.f90 cal_xsq.o params.o quadpack.o
cal_xsq.o:	cal_xsq.f90 params.o quadpack.o
		$(FC) -c cal_xsq.f90 params.o quadpack.o 
quadpack.o:	quadpack.f
		$(FC) -c quadpack.f
params.o:	params.f90 
		$(FC) -c params.f90  
clean:		*.o *.mod
		$(RM) *.o *.mod
