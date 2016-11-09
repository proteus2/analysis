#
#FC = $(FC_DBG)
RM = rm -f
NCINC ?= $(NETCDF)/include
NCLIB ?= $(NETCDF)/lib
NCLIBNAME ?= -lnetcdf -lnetcdff
LIB_FFTPACK ?= $(UTIL)/num_analysis/FFT/fftpack5.1/$(VER_COMP)/lib/libfftpack.a
EZNCINC ?= $(EasyNC)/include
EZNCLIB ?= $(EasyNC)/lib
FLG_EZNC = -I$(EZNCINC) -L$(EZNCLIB) -L$(NCLIB) -leasync $(NCLIBNAME)

EXEC_ALL := vintp2z prtextr prtextr_z fft_x fft_y2 fft_t3 recon_xy recon_xy-w1 recon_xy-w2 recon_xy-w3 recon_xy-w4 recon_xy-w5 psd_wn psd_dir mfsa_wn mfs_dir mean_mean fft_x-wgt intp_rawdata nbe_wrf_z fgf_wrf_z
MOD_N := nrtype nr nrutil
MOD_N := $(MOD_N) fft nbe fgf kl_filter
SUB_N := spline splint tridag_ser locate fdcoef

MOD_F = $(addsuffix .f90,$(MOD_N))
MOD_O = $(MOD_F:.f90=.o)
SUB_F = $(addsuffix .f90,$(SUB_N))
SUB_O = $(SUB_F:.f90=.o)

ALL_SUB_F = $(MOD_F) $(SUB_F)
ALL_SUB_O = $(MOD_O) $(SUB_O)

LIB_SUBS = libsubs.a


# main

vintp2z:	vintp2z.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

prtextr:	prtextr.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

prtextr_z:	prtextr_z.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

fft_x:	fft_x.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

fft_x-wgt:	fft_x-wgt.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

fft_y:	fft_y.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

fft_y2:	fft_y2.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

fft_t3:	fft_t3.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

recon_xy:	recon_xy.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

recon_xy-w1:	recon_xy-w1.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

recon_xy-w2:	recon_xy-w2.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

recon_xy-w3:	recon_xy-w3.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

recon_xy-w4:	recon_xy-w4.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

recon_xy-w5:	recon_xy-w5.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

psd_wn:	psd_wn.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

psd_dir:	psd_dir.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

mfsa_wn:	mfsa_wn.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

mfs_dir:	mfs_dir.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

mean_mean:	mean_mean.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

intp_rawdata:	intp_rawdata.f90  $(LIB_SUBS) $(LIB_FFTPACK)
		$(FC) $< -o $@ $(LIB_SUBS) $(LIB_FFTPACK) $(FLG_EZNC)

nbe_wrf_z:	nbe_wrf_z.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

fgf_wrf_z:	fgf_wrf_z.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) $(FLG_EZNC)

# modules

$(LIB_SUBS):	$(ALL_SUB_F)
		$(FC) -c -I$(NCINC) $(ALL_SUB_F)
		ar -rv $@ $(ALL_SUB_O)
		rm -f $(ALL_SUB_O)

# clean

clean:	
		$(RM) $(EXEC_ALL) *.a *.o *.mod

