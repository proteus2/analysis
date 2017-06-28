#
#FC = $(FC_DBG)
EXEC_ALL := fft_x_ra fft_t2_ra wavedcp_var reconstr_varsa
LIB_FFTPACK ?= $(UTIL)/num_analysis/FFT/fftpack5.1/$(VER_COMP)/lib/libfftpack.a
MOD_DEP_N := m_reanal netio
MOD_N := fft
SUB_N := 
LNK_FFTPACK := fft_x_ra fft_t2_ra

NCINC := $(filter-out -g,$(shell nc-config --fflags))
NCLIB := $(shell nc-config --flibs)

COMPBASE_D := $(VER_COMP)_nc$(filter-out netCDF,$(shell nc-config --version))
EXEC_D := bin/$(COMPBASE_D)
LIB_D := lib/$(COMPBASE_D)
LIB_SUBS := libsubs.a
ifeq ($(findstring gfortran, $(FC) ),gfortran)
  SUBSLIB := -L$(LIB_D) -lsubs -J$(LIB_D)
else
  SUBSLIB := -L$(LIB_D) -lsubs -module $(LIB_D)
endif

ALL_SUB_N := $(MOD_DEP_N) $(MOD_N) $(SUB_N)
ALL_SUB_F = $(addsuffix .f90,$(ALL_SUB_N))
ALL_SUB_O = $(notdir $(addsuffix .o,$(ALL_SUB_N)) )
EXEC_FULL := $(patsubst %,$(EXEC_D)/%,$(EXEC_ALL))

ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$(LNK_FFTPACK)))
  LNK_O += $(LIB_FFTPACK)
endif

# main

$(EXEC_ALL):	%:	$(EXEC_D)/%

$(EXEC_FULL):	$(EXEC_D)/%:	%.f90 $(LIB_D)/$(LIB_SUBS) $(LNK_O)
		mkdir -p $(EXEC_D)
		$(FC) $< -o $@ $(SUBSLIB) $(LNK_O) $(NCLIB)

# modules

$(LIB_D)/$(LIB_SUBS):	$(ALL_SUB_F)
		mkdir -p $(LIB_D)
		( cd $(LIB_D) && $(FC) -c $(patsubst %,$(CURDIR)/%,$(ALL_SUB_F) ) -I$(CURDIR) $(NCINC) )
		( cd $(LIB_D) && ar -rv $(notdir $@ ) $(ALL_SUB_O) )
		( cd $(LIB_D) && rm -f $(ALL_SUB_O) )

# clean

clean:	
		rm -rf $(LIB_D) $(EXEC_D)

