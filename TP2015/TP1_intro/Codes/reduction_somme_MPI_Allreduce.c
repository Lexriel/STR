#include <stdlib.h>
#include <stdio.h>
#include <math.h> // ne pas oublier '-lm' a la compilation 
#include <time.h>
#include "mpi.h"



int main (int argc, char * argv[]) {
  int my_rank; /* mon rang dans MPI_COMM_WORLD */
  int nproc, somme, orig, dest, i, receive;
  MPI_Status status;


  /* Initialisation 
   */
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


  if (my_rank == 0){
    printf("Nb_procs = %i\n", nproc);
  }

  /* Generation du nombre aleatoire : 
   */
  srand48(clock() + my_rank);
  somme = ((double) 111) * drand48();
  printf("my_rank=%d, initial value is %d\n", my_rank, somme);

  /* Calcul de la somme : 
   */
  MPI_Allreduce(&somme, &receive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  printf("Processus %d : somme = %d\n", my_rank, receive);

  /* Terminaison : 
   */
  MPI_Finalize();
}


