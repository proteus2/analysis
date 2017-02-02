CC = gcc
RM = /bin/rm
#
MAINEXE := sidicfp.x
MAINSRC := $(shell echo $(MAINEXE:.x=.c))
MAINOBJ := $(shell echo $(MAINEXE:.x=.o))
#
CSRCS := ncexpr.c
COBJS := $(shell echo $(CSRCS:.c=.o))
#
all: $(MAINOBJ) $(COBJS)
	$(CC) -o $(MAINEXE) $(MAINOBJ) $(COBJS) -lm
clean:
	$(RM) -f $(MAINEXE) $(MAINOBJ) $(COBJS)
#
$(MAINOBJ): $(MAINSRC) $(COBJS)
	$(CC) -I./ -c $<
ncexpr.o: ncexpr.c
	$(CC) -I./ -c $<
sidicfp.o: sidicfp.c ncexpr.o
	$(CC) -I./ -c $<
