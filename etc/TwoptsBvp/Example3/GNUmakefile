FC = ifort
RM = /bin/rm
#
MAINEXE := DTWPBVPBC.x
MAINSRC := $(shell echo $(MAINEXE:.x=.F90))
MAINOBJ := $(shell echo $(MAINEXE:.x=.o))
#
F90SRCS := Base.F90 BVPShared.F90 BVPExtern.F90 BVPC.F90
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
BVPShared.o: BVPShared.F90 Base.o
	$(FC) -I./ -c $<
BVPExtern.o: BVPExtern.F90 BVPShared.o Base.o
	$(FC) -I./ -c $<
BVPC.o: BVPC.F90 BVPExtern.o BVPShared.o Base.o
	$(FC) -I./ -c $<
