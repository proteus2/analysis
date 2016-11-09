#
#
FC = pgf90
RM = /bin/rm
NCINC = $(NETCDF)/include
NCLIB = $(NETCDF)/lib

stat:		stat.o ssapack.o eofpack.o eof_subr.o specanal.o imsl_subr.o svdpack.o svdcmp.o waveletpack.o wavelet.o chisqr.o regress.o random.o netcdfio.o fft.o svdcmp.o
		$(FC) -o stat stat.o ssapack.o eofpack.o eof_subr.o specanal.o imsl_subr.o svdpack.o svdcmp.o waveletpack.o wavelet.o chisqr.o regress.o random.o netcdfio.o -L$(NCLIB) -lnetcdf fft.o
chisqr.o:	chisqr.f
		$(FC) -c chisqr.f
regress.o:	regress.f90
		$(FC) -c regress.f90
random.o:   random.f90
		$(FC) -c random.f90
wavelet.o:	wavelet.f
		$(FC) -c wavelet.f
waveletpack.o:	waveletpack.f90 regress.o netcdfio.o
		$(FC) -c waveletpack.f90 regress.o netcdfio.o
ssapack.o:	ssapack.f90 netcdfio.o
		$(FC) -c ssapack.f90 netcdfio.o
eofpack.o:	eofpack.f90 netcdfio.o fft.o
		$(FC) -c eofpack.f90 netcdfio.o fft.o
eof_subr.o:	eof_subr.f
		$(FC) -c eof_subr.f
svdpack.o:	svdpack.f90 netcdfio.o
		$(FC) -c svdpack.f90 netcdfio.o
svdcmp.o:	svdcmp.f
		$(FC) -c svdcmp.f
fft.o:	fft.f90
		$(FC) -c fft.f90
specanal.o:	specanal.f imsl_subr.o
		$(FC) -c specanal.f imsl_subr.o
imsl_subr.o:	imsl_subr.f
		$(FC) -c imsl_subr.f
netcdfio.o: netcdfio.f90
		$(FC) -c -I$(NCINC) netcdfio.f90
stat.o:	stat.f90 ssapack.o eofpack.o svdpack.o specanal.o waveletpack.o netcdfio.o random.o fft.o
		$(FC) -c -I$(NCINC) stat.f90 netcdfio.o ssapack.o eofpack.o svdpack.o specanal.o waveletpack.o random.o fft.o
clean:	*.o *.mod
		$(RM) *.o *.mod
