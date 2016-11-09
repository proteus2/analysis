FC = /usr/bin/f90
RM = /usr/bin/rm
NCLIB = /usr/local/netcdf/lib
NCINC = /usr/local/netcdf/include
NGLIB = /usr/local/ncarg/lib

psdanal:	psdanal.o varpsd.o varnfit.o pwrspd.o cfftpack.o nonfit.o graphicpsd.o netcdfio.o window.o
		$(FC) -o psdanal psdanal.o varpsd.o varnfit.o pwrspd.o cfftpack.o nonfit.o graphicpsd.o \
		netcdfio.o window.o -L$(NGLIB) -lncarg -lncarg_gks -lncarg_c -lX11 -L$(NCLIB) -lnetcdf
varpsd.o:	varpsd.f90
		$(FC) -c varpsd.f90
varnfit.o:	varnfit.f90
		$(FC) -c varnfit.f90
pwrspd.o:	pwrspd.f90 cfftpack.o window.o 
		$(FC) -c pwrspd.f90 
cfftpack.o:	cfftpack.f90
		$(FC) -c cfftpack.f90
window.o:	window.f90
		$(FC) -c window.f90
nonfit.o:	nonfit.f90
		$(FC) -c nonfit.f90
graphicpsd.o:	graphicpsd.f90 varpsd.o 
		$(FC) -c graphicpsd.f90 
netcdfio.o:	netcdfio.f90
		$(FC) -I$(NCINC) -c netcdfio.f90
psdanal.o:	psdanal.f90 varpsd.o varnfit.o netcdfio.o pwrspd.o cfftpack.o nonfit.o graphicpsd.o window.o 
		$(FC) -c psdanal.f90
clean:		*.o *.mod
		$(RM) *.o *.mod
