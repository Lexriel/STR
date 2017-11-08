#include <stdlib.h>
#include <stdio.h>
#include <math.h> // ne pas oublier '-lm' a la compilation 
#include <time.h>
#include "mpi.h"



int main (int argc, char * argv[]) {
  int my_rank; /* mon rang dans MPI_COMM_WORLD */
  int nproc, somme, orig, dest, i, receive;
  MPI_Status status;
  int nb_etapes = 0;

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
  /*** Ou avec rand()/srand(): 
   * srand() ; // ou : srand(time(NULL)); 
   * somme = rand()%100; 
   */

  printf("my_rank=%d, initial value is %d\n", my_rank, somme);

  /* Calcul de la somme : 
   */
  /* Il y a N etapes, avec : $ 2^{N-1} < nprocs <= 2^N $  */
  nb_etapes = (int) ceil(log(nproc)/log(2));
  if (my_rank == 0) {printf("Nb etapes : %i\n", nb_etapes);}
  for (i=0; i < nb_etapes; ++i) {
    /* Etape 'i' : difference entre source et destination 1<<i 
       (i.e. "1 decale de i", qui vaut $2^i$ ) */
    /* Autre facon de faire : mask=1; puis mask<<=1; */
    
    //    if (my_rank == 0) {printf("i = %i\n", i);}
    
    if ((my_rank >> i) & 1) {
      /* Je dois emettre : */
      dest = my_rank - (1<<i);
      MPI_Send(&somme, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      break; /* Je quitte la boucle 'for' apres mon premier envoi. */
    }
    else {
      /* Je dois recevoir, si mon emetteur a l'etape 'i' existe : */
      if (my_rank + (1<<i) < nproc ) {
	orig = my_rank + (1<<i);
	MPI_Recv(&receive, 1, MPI_INT, orig, 0, MPI_COMM_WORLD,&status);
	somme += receive;
      }
    }
    
  } /* for i */

  if (my_rank==0)
    printf("somme = %d\n", somme);

  /* Terminaison : 
   */
  MPI_Finalize();
}


