FC = ifort
RM = /bin/rm
#
MAINEXE := Example1.x
MAINSRC := $(shell echo $(MAINEXE:.x=.F90))
MAINOBJ := $(shell echo $(MAINEXE:.x=.o))
#
F90SRCS := Base.F90 DataStorage.F90 LMNonFit.F90 LMExtern.F90
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
DataStorage.o: DataStorage.F90 Base.o
	$(FC) -I./ -c $<
LMExtern.o: LMExtern.F90 DataStorage.o Base.o
	$(FC) -I./ -c $<
LMNonFit.o: LMNonFit.F90 LMExtern.o Base.o
	$(FC) -I./ -c $<
