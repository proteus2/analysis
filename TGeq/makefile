FC = pgf90 -mcmodel=medium
RM = /bin/rm -f
IMSLDIR = $(ANAL)/IMSL/routine
NRDIR = $(ANAL)/NR
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

MOD_F = t_g_eqn.f90 netio.f90
MOD_O = t_g_eqn.o netio.o

IMSLSUB_F = $(IMSLDIR)/leq2c.for $(IMSLDIR)/leqt1c.for $(IMSLDIR)/uertst.for $(IMSLDIR)/uspkd.for $(IMSLDIR)/ugetio.for
IMSLSUB_O = comp/leq2c.o comp/leqt1c.o comp/uertst.o comp/uspkd.o comp/ugetio.o
NRSUB_F = $(NRDIR)/nrtype.f90 $(NRDIR)/nr.f90 $(NRDIR)/nrutil.f90 $(NRDIR)/odeint.f90 $(NRDIR)/bsstep.f90 $(NRDIR)/stifbs.f90 $(NRDIR)/stoerm.f90 $(NRDIR)/mmid.f90 $(NRDIR)/pzextr.f90 $(NRDIR)/rzextr.f90 $(NRDIR)/simpr.f90 $(NRDIR)/ludcmp.f90 $(NRDIR)/lubksb.f90 $(NRDIR)/locate.f90 $(NRDIR)/spline.f90 $(NRDIR)/splint.f90 $(NRDIR)/tridag_ser.f90 $(NRDIR)/lfit.f90 $(NRDIR)/fpoly.f90 $(NRDIR)/gaussj.f90 $(NRDIR)/covsrt.f90
NRSUB_O = comp/nrtype.o comp/nr.o comp/nrutil.o comp/odeint.o comp/bsstep.o comp/stifbs.o comp/stoerm.o comp/mmid.o comp/pzextr.o comp/rzextr.o comp/simpr.o comp/ludcmp.o comp/lubksb.o comp/locate.o comp/spline.o comp/splint.o comp/tridag_ser.o comp/lfit.o comp/fpoly.o comp/gaussj.o comp/covsrt.o

#

run:		tgeq.o $(IMSLSUB_O) $(NRSUB_O) MOD_C
		$(FC) -o run tgeq.o $(IMSLSUB_O) $(NRSUB_O) $(MOD_O) -L$(NCLIB) -lnetcdf

PRECOMP:	$(IMSLSUB_F) $(NRSUB_F)
		$(RM) run *.o *.mod
		$(FC) -c $(IMSLSUB_F) $(NRSUB_F)
		mv *.o comp/

tgeq.o:	tgeq.f90 MOD_C
		$(FC) -c tgeq.f90 $(MOD_O) -I$(NCINC)

MOD_C:	$(MOD_F)
		$(FC) -c $(MOD_F) -I$(NCINC)

clean:
		$(RM) run *.o *.mod

