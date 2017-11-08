############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP8 : Equation de poisson
##
## Auteur          : Dimitri Lecas (CNRS/IDRIS)
## 
###############################################################################

# Compilateur et options de compilation et de link
include ../../arch/make_inc

OBJS = parallel.o poisson.o calcul.o
OBJS1 = read.o
OBJS2 = calcul_exact.o
OBJS3 = module_params.o module_parallel_mpi.o module_calcul.o poisson_non_bloquant.o

# Règle implicite de compilation
.SUFFIXES: .o .c
.c.o:
	$(CC) -c  $(CCFLAGS_TP8) $<

all: poisson

poisson: $(OBJS)
	$(CC) $(LDFLAGS_TP8) -o $@ $(OBJS)

poisson_non_bloquant:$(OBJS3)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS3)

read :  $(OBJS1)
	$(CC) $(LDFLAGS_TP8) -o $@ $(OBJS1)

calcul_exact : $(OBJS2)
	$(CC) $(LDFLAGS_TP8) -o $@ $(OBJS2)

verification : fort.11 fort.10
	diff fort.11 fort.10

fort.11 : read donnees.dat
	$(EXEC_TP8_VERIF) ./read

fort.10 : calcul_exact
	$(EXEC_TP8_VERIF) ./calcul_exact

donnees.dat : exe

exe: poisson
	$(MPIEXEC_TP8)	./poisson

clean :
	rm -f $(OBJS) $(OBJS1) $(OBJS2) $(OBJS3) poisson poisson_non_bloquant calcul_exact read *.mod core donnees.dat fort.10 fort.11
