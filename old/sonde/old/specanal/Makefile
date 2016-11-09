FC = /usr/bin/f90
RM = /usr/bin/rm
NCLIB = /usr/local/netcdf/lib
NCINC = /usr/local/netcdf/include
NGLIB = /usr/local/ncarg/lib

specanal:	specanal.o varspec.o hilbert.o meanwavnum.o propagate.o window.o cfftpack.o momentum.o \
		netcdfio.o graphicspec.o
		$(FC) -o specanal specanal.o varspec.o hilbert.o meanwavnum.o propagate.o window.o \
		cfftpack.o momentum.o netcdfio.o graphicspec.o -L$(NCLIB) -lnetcdf  \
		-L$(NGLIB) -lncarg -lncarg_gks -lncarg_c -lX11
varspec.o:	varspec.f90
		$(FC) -c varspec.f90
hilbert.o:	hilbert.f90
		$(FC) -c hilbert.f90
window.o:	window.f90
		$(FC) -c window.f90
cfftpack.o:	cfftpack.f90
		$(FC) -c cfftpack.f90
meanwavnum.o:	meanwavnum.f90 window.o cfftpack.o
		$(FC) -c meanwavnum.f90
propagate.o:	propagate.f90 window.o cfftpack.o
		$(FC) -c propagate.f90
momentum.o:	momentum.f90 
		$(FC) -c momentum.f90
netcdfio.o:	netcdfio.f90
		$(FC) -I$(NCINC) -c netcdfio.f90
graphicspec.o:	graphicspec.f90 varspec.o
		$(FC) -c graphicspec.f90
specanal.o:	specanal.f90 varspec.o hilbert.o meanwavnum.o propagate.o window.o cfftpack.o momentum.o \
		netcdfio.o graphicspec.o
		$(FC) -c specanal.f90
clean:		*.mod *.o
		$(RM) *.mod *.o
