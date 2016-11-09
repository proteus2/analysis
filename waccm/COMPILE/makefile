#
#FC = $(FC_DBG)
RM = rm -f
NCINC ?= $(NETCDF)/include
NCLIB ?= $(NETCDF)/lib
NCLIBNAME ?= -lnetcdf -lnetcdff
LIB_FFTPACK ?= $(UTIL)/num_analysis/FFT/fftpack5.1/$(VER_COMP)/lib/libfftpack.a

EXEC_ALL := tem_mm_cam
MOD_N := nrtype nr nrutil m_cam netio
MOD_N := $(MOD_N) fft tem
SUB_N := spline splint tridag_ser locate fdcoef

MOD_F = $(addsuffix .f90,$(MOD_N))
MOD_O = $(MOD_F:.f90=.o)
SUB_F = $(addsuffix .f90,$(SUB_N))
SUB_O = $(SUB_F:.f90=.o)

ALL_SUB_F = $(MOD_F) $(SUB_F)
ALL_SUB_O = $(MOD_O) $(SUB_O)

LIB_SUBS = libsubs.a

F_INC = c_math.inc c_phys.inc


# main

tem_mm_cam:	tem_mm_cam.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

# modules

$(LIB_SUBS):	$(ALL_SUB_F)
		$(FC) -c -I$(NCINC) $(ALL_SUB_F)
		ar -rv $@ $(ALL_SUB_O)
		rm -f $(ALL_SUB_O)

# clean

clean:	
		$(RM) $(EXEC_ALL) *.a *.o *.mod

