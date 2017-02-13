FC = gfortran
RM = /bin/rm
#
MAINEXE := Example.x
MAINSRC := $(shell echo $(MAINEXE:.x=.F90))
MAINOBJ := $(shell echo $(MAINEXE:.x=.o))
#
F90SRCS := Base.F90 InvLaplace.F90 InvLaplaceExtern.F90
F90OBJS := $(shell echo $(F90SRCS:.F90=.o))
F90MODS := $(shell echo $(F90SRCS:.F90=.mod) | tr A-Z a-z)
#
all: $(MAINOBJ) $(F90OBJS)
	$(FC) -o $(MAINEXE) $(MAINOBJ) $(F90OBJS)
clean:
	$(RM) -f $(MAINEXE) $(MAINOBJ) $(F90OBJS) $(F90MODS)
#
$(MAINOBJ): $(MAINSRC) $(F90OBJS)
	$(FC) -I./ -c $<
Base.o: Base.F90
	$(FC) -I./ -c $<
InvLaplace.o: InvLaplace.F90 Base.o InvLaplaceExtern.o
	$(FC) -I./ -c $<
InvLaplaceExtern.o: InvLaplaceExtern.F90 Base.o
	$(FC) -I./ -c $<
