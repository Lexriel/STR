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






  mon_germe();




    tirage = pile_ou_face();



  if (rang == 0) {
    if (nb_essais > nb_essais_max)
      printf("On a fait %d essais sans faire l'unanimite !\n", nb_essais_max);
    else
      printf("On a fait %d essais avant de faire l'unanimite !\n", nb_essais);
  }


  return 0;
}
