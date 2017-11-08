/*
 * unanimite.f90 --- TP3 : Diffusions globales et réductions
 *
 * Auteur          : Dimitri LECAS (CNRS/IDRIS - France) <dimitri.lecas@idris.fr>
 *
 */

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int rang, nb_processus;

void mon_germe() {
  srand(time(NULL)+rang);
}

int pile_ou_face() {
  return (int)((double)rand() / ((double)RAND_MAX + 1) * 2);
}

int main(int argc, char *argv[]) {
  int decision = -1;
  int nb_essais_max = 2000;
  int nb_essais = 0;
  int tirage;

  MPI_Init( &argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processus);

  mon_germe();

  while ( (decision != 0) && (decision != nb_processus) &&
	  (nb_essais <= nb_essais_max) ) {
    nb_essais++;
    tirage = pile_ou_face();
    MPI_Allreduce(&tirage,&decision,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }

  if (rang == 0) {
    if (nb_essais > nb_essais_max)
      printf("On a fait %d essais sans faire l'unanimite !\n", nb_essais_max);
    else
      printf("On a fait %d essais avant de faire l'unanimite !\n", nb_essais);
  }
  MPI_Finalize();
  return 0;
}
