############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP3 : Communications collectives
##
## Auteur          : Dimitri Lecas (CNRS/IDRIS) <dimitri.lecas@idris.fr>
## 
###############################################################################

# Compilateurs, options de compilation et d'édition de liens
include ../arch/make_inc

OBJS  = unanimite.o

# Règle implicite de compilation
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CCLAGS_TP3) $<

all: unanimite

unanimite: $(OBJS)
	$(CC) -o $@ $(LDFLAGS_TP3) $(OBJS)
	$(MPIEXEC_TP3) ./unanimite

clean:
	rm -f $(OBJS) unanimite core
