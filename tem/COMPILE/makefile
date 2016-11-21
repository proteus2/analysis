#
FC = $(FC_DBG) # -mcmodel=medium -shared-intel
EXEC_ALL := tem_ra tem_ra-lastmiss tem_qg tem_qg_betap_gp  \
            epf_fc2_ra reconstr_epf epfsa0_fc2_ra epfsa1_fc2_ra  \
            reconstr_epfsa reconstr_epfsa-a reconstr_epfsa1 reconstr_epfsa2  \
            warc3d_qg_ra warc3d_qg_gp_ra
MOD_DEP_N := m_reanal netio
MOD_N := tem epf warc3d
SUB_N := fdcoef
INCL_c_phys := tem_ra tem_ra-lastmiss tem_qg tem_qg_betap_gp  \
               epf_fc2_ra epfsa0_fc2_ra epfsa1_fc2_ra  \
               warc3d_qg_ra warc3d_qg_gp_ra

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

ifneq ($(ALL_SUB_N),$(filter-out $(INCL_c_phys),$(ALL_SUB_N)))
  INCL_F_SUB += c_phys.inc
endif
ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$(INCL_c_phys)))
  INCL_F += c_phys.inc
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
		( cd $(LIB_D) ; $(FC) -c $? $(NCINC) )
		( cd $(LIB_D) ; ar -rv $(notdir $@ ) $(?:.f90=.o) )
		( cd $(LIB_D) ; rm -f $? $(?:.f90=.o) )

# clean

clean:	
		rm -rf $(LIB_D) $(EXEC_D)
		rmdir $(COMPBASE_D)

