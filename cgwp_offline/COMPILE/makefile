#
#FC = $(FC_DBG) # -mcmodel=medium -shared-intel

EXEC_ALL := offline_hg2 offline_hg2_ra test test1
NR_N := nrtype nr nrutil spline splint tridag_ser locate svdcmp-sq svbksb pythag
MOD_DEP_N := param_gwp switch_dump netio m_common m_hadgem
MOD_N := mflx-ctop_sc05 prop_diss m_zonalavg
SUB_N := read_erai_lon_dat
INCL_c_math := prop_diss offline_hg2 offline_hg2_ra test test1

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

ALL_SUB_N := $(addprefix NR/,$(NR_N)) $(MOD_DEP_N) $(MOD_N) $(SUB_N)
ALL_SUB_F := $(addsuffix .f90,$(ALL_SUB_N))
ALL_SUB_O := $(notdir $(addsuffix .o,$(ALL_SUB_N)) )
EXEC_FULL := $(patsubst %,$(EXEC_D)/%,$(EXEC_ALL))

ifneq ($(ALL_SUB_N),$(filter-out $(INCL_c_math),$(ALL_SUB_N)))
  INCL_F_SUB += c_math.inc
endif
ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$(INCL_c_math)))
  INCL_F += c_math.inc
endif

# main

$(EXEC_ALL):	%:	$(EXEC_D)/%

$(EXEC_FULL):	$(EXEC_D)/%:	%.f90 $(LIB_D)/$(LIB_SUBS) $(INCL_F)
		mkdir -p $(EXEC_D)
		$(FC) $< -o $@ $(SUBSLIB) $(NCLIB)

# library

$(LIB_D)/$(LIB_SUBS):	$(ALL_SUB_F) $(INCL_F_SUB)
		mkdir -p $(LIB_D)
		cp $? $(LIB_D)/
		( cd $(LIB_D) ; $(FC) -c $(notdir $(ALL_SUB_F) ) $(NCINC) )
		( cd $(LIB_D) ; ar -rv $(notdir $@ ) $(ALL_SUB_O) )
		( cd $(LIB_D) ; rm -f $(notdir $(ALL_SUB_F) ) $(ALL_SUB_O) )

# clean

clean:	
		rm -rf $(LIB_D) $(EXEC_D)

