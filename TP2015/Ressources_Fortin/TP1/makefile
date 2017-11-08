CC       = mpicc
CFLAGS   =
LDFLAGS  = -lX11
HOSTFILE = hosts
MPIRUN   = mpirun -n 2
OPTS     = --display-map -hostfile $(HOSTFILE)

# exécutable qui sera lancé par la commande make run
EXE      = ./exercice1


# règle pour construire des fichiers objets
%.o: %.c
	$(CC) -c $(CFLAGS) $(LDFLAGS) $<


# exécution conventionnelle de la commande make
all: exercice1 exercice2

# exécutables
exercice1: exercice1.c
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $<
exercice2: exercice2.c
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $<

# faire tourner un des exécutables
run:
	$(MPIRUN) $(OPTS) $(EXE)
# mpirun -n 2 --display-map -hostfile hosts

# règles de nettoyage
.PHONY: clean realclean

clean:
	rm -rf *~ *.o

realclean:
	rm -rf *~ *.o
	rm -rf exercice1 exercice2
