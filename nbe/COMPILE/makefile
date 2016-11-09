#
#FC = $(FC_DBG)
RM = rm -f
NCINC ?= $(NETCDF)/include
NCLIB ?= $(NETCDF)/lib
NCLIBNAME ?= -lnetcdf -lnetcdff

EXEC_ALL := nbe_ra
MOD_N := m_reanal netio
MOD_N := $(MOD_N) nbe
SUB_N := 

MOD_F = $(addsuffix .f90,$(MOD_N))
MOD_O = $(MOD_F:.f90=.o)
SUB_F = $(addsuffix .f90,$(SUB_N))
SUB_O = $(SUB_F:.f90=.o)

ALL_SUB_F = $(MOD_F) $(SUB_F)
ALL_SUB_O = $(MOD_O) $(SUB_O)

LIB_SUBS = libsubs.a


# main

nbe_ra:	nbe_ra.f90  $(LIB_SUBS)
		$(FC) $< -o $@ $(LIB_SUBS) -L$(NCLIB) $(NCLIBNAME)

# modules

$(LIB_SUBS):	$(ALL_SUB_F)
		$(FC) -c -I$(NCINC) $(ALL_SUB_F)
		ar -rv $@ $(ALL_SUB_O)
		rm -f $(ALL_SUB_O)

# clean

clean:	
		$(RM) $(EXEC_ALL) *.a *.o *.mod

