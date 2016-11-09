#
#
FC = pgf90
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

DIR_HERE = $(PWD)
SUBD = $(ANALSUBS)

MOD_F = m_deriv.f90 m_integ.f90 m_avg.f90 m_utiletc.f90 m_temeq.f90
MOD_O = m_deriv.o   m_integ.o   m_avg.o   m_utiletc.o   m_temeq.o

SUB_F =
SUB_O = 

ALL_SUB_F = $(MOD_F) $(SUB_F)

F_INC = c_math.inc c_phys.inc


# exec

runn:		MAIN_COMP SUBS_COMP
		( cd COMPILE  ;  \
		  $(FC) -o runn main.o $(MOD_O) $(SUB_O) -L$(NCLIB) -lnetcdf )
		cp COMPILE/runn ./

# main

MAIN_COMP:	wr.F90 SUBS_COMP COMPILE
		cp wr.F90 COMPILE/main.F90
		( cd COMPILE  ;  \
		  $(FC) -c -I$(NCINC) main.F90 $(MOD_O) )

COMPILE:	
		mkdir COMPILE

# modules

SUBS_COMP:	SUBSRC
		( cd COMPILE  ;  \
		  $(FC) -c $(ALL_SUB_F) )

SUBSRC:	COMPILE
		( cd $(ANALSUBS)  ;  \
		  cp $(ALL_SUB_F) $(F_INC) $(DIR_HERE)/COMPILE/ )

# clean

clean:	COMPILE
		$(RM) runn COMPILE/*.o COMPILE/*.mod COMPILE/runn

