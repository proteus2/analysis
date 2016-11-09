FC = pgf90 -mcmodel=medium
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

MOD_F = netcdfio.f90 fft.f90 regress.f90 parafit2d.f90
MOD_O = netcdfio.o   fft.o   regress.o   parafit2d.f90


t2mf:				cp_t2mf run
m2:				cp_m2 run
mf_c_dir:			cp_mf_c_dir run
mf_c_dir_plot-data:	cp_mf_c_dir_plot-data run
mfa-klf:			cp_mfa-klf run
mfa-c_phi:			cp_mfa-c_phi run
test:				cp_test run

#

cp_t2mf:			t2mf.f90
		cp t2mf.f90 run.f90
cp_m2:			m2.f90
		cp m2.f90 run.f90
cp_mf_c_dir:		mf_c_dir.f90
		cp mf_c_dir.f90 run.f90
cp_mf_c_dir_plot-data:	mf_c_dir_plot-data.f90
		cp mf_c_dir_plot-data.f90 run.f90
cp_mfa-klf:			mfa-klf.f90
		cp mfa-klf.f90 run.f90
cp_mfa-c_phi:		mfa-c_phi.f90
		cp mfa-c_phi.f90 run.f90
cp_test:			test.f90
		cp test.f90 run.f90

#

run:		run.o MOD_C
		$(FC) -o run run.o $(MOD_O) -L$(NCLIB) -lnetcdf

run.o:	run.f90 MOD_C
		$(FC) -c run.f90 $(MOD_O) -I$(NCINC)

MOD_C:	$(MOD_F)
		$(FC) -c $(MOD_F) -I$(NCINC)

clean:
		$(RM) run *.o *.mod

