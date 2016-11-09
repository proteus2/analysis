#
#FC = $(FC_DBG)
RM = rm -f
NCINC ?= $(NETCDF)/include
NCLIB ?= $(NETCDF)/lib
NCLIBNAME ?= -lnetcdf $(LDFLAG_NC2)
LIB_FFTPACK ?= $(UTIL)/num_analysis/FFT/fftpack5.1/$(VER_COMP)/lib/libfftpack.a

EXEC_ALL := fft_x_ra fft_t2_ra wavedcp_var reconstr_varsa test1
MOD_N := m_reanal netio
MOD_N := $(MOD_N) fft
SUB_N := 

MOD_F = $(addsuffix .f90,$(MOD_N))
MOD_O = $(MOD_F:.f90=.o)
SUB_F = $(addsuffix .f90,$(SUB_N))
SUB_O = $(SUB_F:.f90=.o)

ALL_SUB_F = $(MOD_F) $(SUB_F)
ALL_SUB_O = $(MOD_O) $(SUB_O)

LIB_SUBS = libsubs.a


# main

fft_x_ra:	fft_x_ra.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

fft_t2_ra:	fft_t2_ra.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

wavedcp_var:	wavedcp_var.f90 $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

reconstr_varsa:	reconstr_varsa.f90 $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

test1:	test1.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

# modules

$(LIB_SUBS):	$(ALL_SUB_F)
		$(FC) -c -I$(NCINC) $(ALL_SUB_F)
		ar -rv $@ $(ALL_SUB_O)
		rm -f $(ALL_SUB_O)

# clean

clean:	
		$(RM) $(EXEC_ALL) *.a *.o *.mod

