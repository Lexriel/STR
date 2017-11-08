CC       = mpicc
CFLAGS   =
LDFLAGS  = -lX11
HOSTFILE = hosts
MPIRUN   = mpirun -n 16
OPTS     = --display-map # -hostfile $(HOSTFILE)

OBJS     = outils.o produit_Ax_lignes.o
OBJS2    = outils.o vecteur_propre_lignes.o
OBJS3    = outils.o suite_produits_lignes.o
SRCS     = outils.c vecteur_propre_lignes.c produit_Ax_lignes.c suite_produits_lignes.c
SRCS     = ${SRCS:.c=.o}

# exécutable qui sera lancé par la commande make run
EXE      = ./produit_Ax_lignes


# règle pour construire des fichiers objets
%.o : %.c
	$(CC) -c $(CFLAGS) $(LDFLAGS) $<


# exécution conventionnelle de la commande make
all : produit_Ax_lignes vecteur_propre_lignes suite_produits_lignes


# exécutables
produit_Ax_lignes : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(OBJS) -lm

vecteur_propre_lignes : $(OBJS2)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(OBJS2) -lm

suite_produits_lignes : $(OBJS3)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(OBJS3) -lm


# faire tourner un des exécutables
run :
	$(MPIRUN) $(OPTS) $(EXE)


# règles de nettoyage
.PHONY : clean realclean

clean :
	rm -rf *~ *.o

realclean :
	rm -rf *~ *.o
	rm -rf produit_Ax_lignes vecteur_propre_lignes suite_produits_lignes
