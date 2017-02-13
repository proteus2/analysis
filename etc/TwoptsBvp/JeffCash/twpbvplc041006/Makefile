.SUFFIXES: .o .f

FORTRAN = f77
# suggested for DEC Alpha FFLAGS = -O3 -fpe1
FFLAGS =
LFLAGS =

CFILES = dtwpbvpbl.o dtwpbvpbc.o twpbvpl.o twpbvpc.o \
	dtwpbvpbl dtwpbvpbc

all:		dtwpbvpbl dtwpbvpbc 

dtwpbvpbl:	twpbvpl.o dtwpbvpbl.o
		$(FORTRAN) $(LFLAGS) -o dtwpbvpbl twpbvpl.o dtwpbvpbl.o

dtwpbvpbc:	twpbvpc.o dtwpbvpbc.o
		$(FORTRAN) $(LFLAGS) -o dtwpbvpbc twpbvpc.o dtwpbvpbc.o

clean:
		touch $(CFILES) 
		rm $(CFILES) 

.f.o:
	$(FORTRAN) $(FFLAGS) -c $*.f
