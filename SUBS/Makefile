#
#
MODNAME = integ


FC = pgf90
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

MOD_F = m_$(MODNAME).f90
MOD_O = m_$(MODNAME).o  


# exec

test:		TEST_C MOD_C NET_C
		( cd TEST  ; \
		  $(FC) -o test test.o $(MOD_O) netio.o -L$(NCLIB) -lnetcdf  ; \
		  mv test .. )

# main

TEST_C:	test.f90 MOD_C NET_C
		( cp test.f90 TEST/  ; \
		  cd TEST  ; \
		  $(FC) -c test.f90 $(MOD_O) netio.o -I$(NCINC) )

# modules

MOD_C:	$(MOD_F)
		( cp $(MOD_F) TEST/  ; \
		  cd TEST  ; \
		  $(FC) -c $(MOD_F) )

NET_C:	netio.f90
		( cp netio.f90 TEST/  ; \
		  cd TEST  ; \
		  $(FC) -c netio.f90 -I$(NCINC) )

# clean

clean:
		$(RM) test TEST/*.o TEST/*.mod

