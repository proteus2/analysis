#
#
FC = pgf90
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

MOD_F = m_deriv.f90 m_integ.f90 m_avg.f90 m_utiletc.f90 m_regrid.f90
MOD_O = m_deriv.o   m_integ.o   m_avg.o   m_utiletc.o   m_regrid.o

SUB_F = 
SUB_O =

ALL_SUB_F = $(MOD_F) $(SUB_F)

F_INC = c_math.inc c_phys.inc


# exec

runn:		anal.o  netio.o um_anal.o um_axis.o pwrspd.o fft.o regress.o mud2crf90.o mudcom.o SUBS_COMP
		$(FC) -o runn anal.o netio.o -L$(NCLIB) -lnetcdf um_anal.o um_axis.o pwrspd.o fft.o regress.o mud2crf90.o mudcom.o $(MOD_O) $(SUB_O)

# main

anal.o:	anal.f90  netio.o um_anal.o um_axis.o pwrspd.o fft.o mud2crf90.o SUBS_COMP
		$(FC) -c -I$(NCINC) anal.f90

# modules

netio.o:	netio.f90
		$(FC) -c -I$(NCINC) netio.f90

um_anal.o:	um_anal.f90
		$(FC) -c um_anal.f90

um_axis.o:	um_axis.f90 netio.o
		$(FC) -c -I$(NCINC) um_axis.f90

regress.o:	regress.f90
		$(FC) -c regress.f90

fft.o:	fft.f90
		$(FC) -c fft.f90

pwrspd.o:	pwrspd.f90 fft.o regress.o
		$(FC) -c pwrspd.f90

mudcom.o:	mudcom.F
		$(FC) -c mudcom.F
mud2crf90.o:	mud2crf90.F mudcom.o
		$(FC) -c mud2crf90.F

SUBS_COMP:	$(F_INC)
		$(FC) -c $(ALL_SUB_F)

# clean

clean:	
		$(RM) runn *.o *.mod
