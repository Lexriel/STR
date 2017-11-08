############################ -*- Mode: Makefile -*- ###########################
## Makefile --- T.P. 7 du cours MPI
##
## Auteur          : Dimitri Lecas (CNRS/IDRIS - France)
##
###############################################################################

# Compilateur et options de compilation et de link
include ../../arch/make_inc

OBJS    = lire_fichier.o

# Règle implicite de compilation
.SUFFIXES:
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CCFLAGS_TP7) $<

lire_fichier:creation_fichier $(OBJS)
	$(CC) $(LDFLAGS_TP7) -o $@ $(OBJS)
	$(MPIEXEC_TP7) ./lire_fichier
	rm -f $(OBJS) core

creation_fichier:
	cd creation;make clean;make

exec: lire_fichier
	$(MPIEXEC_TP7) ./lire_fichier
	rm -f $(OBJS) core

verification:
	$(CC) verification.c $(LDCFLAGS) -o verification
	$(EXEC_TP7_VERIF) ./verification
	rm -f verification

clean:
	rm -f *.o lire_fichier fichier*.dat verification core*

cleanall:
	make clean
	rm -f donnees.dat
