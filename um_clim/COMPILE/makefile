#
#FC = $(FC_DBG)
RM = rm -f
NCINC ?= $(NETCDF)/include
NCLIB ?= $(NETCDF)/lib
NCLIBNAME ?= -lnetcdf -lnetcdff
#NCLIBNAME ?= -lnetcdf
LIB_FFTPACK ?= $(UTIL)/num_analysis/FFT/fftpack5.1/$(VER_COMP)/lib/libfftpack.a

EXEC_ALL := v_intp w-cont_um tem_z_um tem_z_cp_um epf_z_um fft_x_um fft_x-wgt_um fft_t2_um fft_t2i_um reconstr_var epf_z_fc2_um epf_z_cp_fc2_um epfsa0_z_fc2_um epfsa1_z_fc2_um epfsa0_z_cp_fc2_um epfsa1_z_cp_fc2_um reconstr_epf reconstr_varsa wavedcp_var reconstr_epfsa reconstr_epfsa1
MOD_N := nrtype nr nrutil m_hadgem netio
MOD_N := $(MOD_N) fft tem epf
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

w-cont_um:	w-cont_um.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

tem_z_um:	tem_z_um.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

tem_z_cp_um:	tem_z_cp_um.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

epf_z_um:	epf_z_um.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

fft_x_um:	fft_x_um.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

fft_x-wgt_um:	fft_x-wgt_um.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

fft_t2_um:	fft_t2_um.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

fft_t2i_um:	fft_t2i_um.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_math.inc c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

reconstr_var:	reconstr_var.f90  $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

reconstr_varsa:	reconstr_varsa.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

wavedcp_var:	wavedcp_var.f90 $(LIB_SUBS) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

epf_z_fc2_um:	epf_z_fc2_um.f90 $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

epf_z_cp_fc2_um:	epf_z_cp_fc2_um.f90 $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

epfsa0_z_fc2_um:	epfsa0_z_fc2_um.f90 $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

epfsa1_z_fc2_um:	epfsa1_z_fc2_um.f90 $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

epfsa0_z_cp_fc2_um:	epfsa0_z_cp_fc2_um.f90 $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

epfsa1_z_cp_fc2_um:	epfsa1_z_cp_fc2_um.f90 $(LIB_SUBS) $(LIB_FFTPACK) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) -L$(NCLIB) $(NCLIBNAME)

reconstr_epf:	reconstr_epf.f90 $(LIB_SUBS) c_phys.inc
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

reconstr_epfsa:	reconstr_epfsa.f90 $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

reconstr_epfsa1:	reconstr_epfsa1.f90 $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

v_intp:		v_intp.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

# modules

$(LIB_SUBS):	$(ALL_SUB_F)
		$(FC) -c -I$(NCINC) $(ALL_SUB_F)
		ar -rv $@ $(ALL_SUB_O)
		rm -f $(ALL_SUB_O)

# clean

clean:	
		$(RM) $(EXEC_ALL) *.a *.o *.mod

