#
FC = $(FC_DBG) # -mcmodel=medium -shared-intel

EXEC_ALL := test
MOD_DEP_N := param_gwp switch_dump netio m_common
MOD_N := mflx-ctop_sc05 prop_diss
SUB_N := 
INCL_c_math := prop_diss test

NCINC := $(filter-out -g,$(shell nc-config --fflags))
NCLIB := $(shell nc-config --flibs)

COMPBASE_D := $(VER_COMP)_nc$(filter-out netCDF,$(shell nc-config --version))
EXEC_D := $(COMPBASE_D)/bin
LIB_D := $(COMPBASE_D)/lib
LIB_SUBS := libsubs.a
ifeq ($(findstring gfortran, $(FC) ),gfortran)
  SUBSLIB := -L$(LIB_D) -lsubs -J$(LIB_D)
else
  SUBSLIB := -L$(LIB_D) -lsubs -module $(LIB_D)
endif

ALL_SUB_N := $(MOD_DEP_N) $(MOD_N) $(SUB_N)
ALL_SUB_F := $(addsuffix .f90,$(ALL_SUB_N))
ALL_SUB_O := $(ALL_SUB_F:.f90=.o)
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
		( cd $(LIB_D) ; $(FC) -c $(ALL_SUB_F) $(NCINC) )
		( cd $(LIB_D) ; ar -rv $(notdir $@ ) $(ALL_SUB_O) )
		( cd $(LIB_D) ; rm -f $(ALL_SUB_F) $(ALL_SUB_O) )

# clean

clean:	
		rm -rf $(LIB_D) $(EXEC_D)
		rmdir $(COMPBASE_D)

