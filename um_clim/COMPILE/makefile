#
#FC = $(FC_DBG)

EXEC_ALL := v_intp w-cont_um tem_z_um tem_z_cp_um epf_z_um fft_x_um fft_x-wgt_um fft_t2_um fft_t2i_um reconstr_var epf_z_fc2_um epf_z_cp_fc2_um epfsa0_z_fc2_um epfsa1_z_fc2_um epfsa0_z_cp_fc2_um epfsa1_z_cp_fc2_um reconstr_epf reconstr_varsa wavedcp_var reconstr_epfsa reconstr_epfsa1 dp_um vintp_z2p_um qboi_dyn_um
MOD_DEP_N := nrtype nr nrutil m_hadgem netio
MOD_N := fft tem epf
SUB_N := spline splint tridag_ser locate fdcoef
LIB_FFTPACK ?= $(UTIL)/num_analysis/FFT/fftpack5.1/$(VER_COMP)/lib/libfftpack.a
INCL_c_math := fft_t2i_um qboi_dyn_um
INCL_c_phys := tem_z_um tem_z_cp_um  \
               fft_x_um fft_x-wgt_um fft_t2_um fft_t2i_um reconstr_var  \
               epf_z_fc2_um epf_z_cp_fc2_um epfsa0_z_fc2_um epfsa1_z_fc2_um  \
               epfsa0_z_cp_fc2_um epfsa1_z_cp_fc2_um reconstr_epf wavedcp_var  \
               vintp_z2p_um qboi_dyn_um
USE_FFTPACK := tem_z_um tem_z_cp_um epf_z_um  \
               fft_x_um fft_x-wgt_um fft_t2_um fft_t2i_um reconstr_var  \
               epf_z_fc2_um epf_z_cp_fc2_um epfsa0_z_fc2_um epfsa1_z_fc2_um  \
               epfsa0_z_cp_fc2_um epfsa1_z_cp_fc2_um reconstr_varsa

NCINC := $(filter-out -g,$(shell nc-config --fflags))
NCLIB := $(shell nc-config --flibs)

COMPBASE_D := $(VER_COMP)_nc$(filter-out netCDF,$(shell nc-config --version))
EXEC_D := $(COMPBASE_D)/bin
LIB_D := $(COMPBASE_D)/lib
LIB_SUBS := libsubs.a
SUBSLIB := -L$(LIB_D) -lsubs -module $(LIB_D)

ALL_SUB_N := $(MOD_DEP_N) $(MOD_N) $(SUB_N)
ALL_SUB_F := $(addsuffix .f90,$(ALL_SUB_N))
ALL_SUB_O := $(ALL_SUB_F:.f90=.o)
EXEC_FULL := $(patsubst %,$(EXEC_D)/%,$(EXEC_ALL))

ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$(INCL_c_math)))
  INCL_F += c_math.inc
endif
ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$(INCL_c_phys)))
  INCL_F += c_phys.inc
endif
ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$(USE_FFTPACK)))
  LIB_EXTERNAL += $(LIB_FFTPACK)
endif

# main

$(EXEC_ALL):	%:	$(EXEC_D)/%

$(EXEC_FULL):	$(EXEC_D)/%:	%.f90 $(LIB_D)/$(LIB_SUBS) $(INCL_F)  \
                                $(LIB_EXTERNAL)
		mkdir -p $(EXEC_D)
		$(FC) $< -o $@ $(SUBSLIB) $(NCLIB) $(LIB_EXTERNAL)

# library

$(LIB_D)/$(LIB_SUBS):	$(ALL_SUB_F)
		mkdir -p $(LIB_D)
		cp $? $(LIB_D)/
		( cd $(LIB_D) ; $(FC) -c $? $(NCINC) )
		( cd $(LIB_D) ; ar -rv $(notdir $@ ) $(?:.f90=.o) )
		( cd $(LIB_D) ; rm -f $? $(?:.f90=.o) )

# clean

clean:	
		rm -rf $(LIB_D) $(EXEC_D)
		rmdir $(COMPBASE_D)

