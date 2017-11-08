# include <stdio.h> 
# include <mpi.h> 

int main(int argc, char* argv[]) { 
  int rang, Nprocs, val, val_recue, cible, etiquette; 
  double debut, fin; 
  MPI_Status status; 

  MPI_Init(&argc, &argv); 
  debut = MPI_Wtime(); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rang); 
  MPI_Comm_size(MPI_COMM_WORLD, &Nprocs); 

  if (rang == 0){
    val = 512;
    etiquette = 1;
    cible = 3;
    MPI_Send(&val, 1, MPI_INT, cible, etiquette, MPI_COMM_WORLD);
    printf("Processus[%d] : envoie la valeur %d.\n", rang, val);
  }
  else if (rang == 3){
    MPI_Recv(&val_recue, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
	     MPI_COMM_WORLD, &status); 
    printf("Processus[%d] : reçoit la valeur %d.\n", rang, val_recue);
  }

  fin = MPI_Wtime(); 
  printf("Processus[%d] : temps d'exécution = %fs \n", rang, fin-debut); 
  
  MPI_Finalize(); 
  return 0; 
}
