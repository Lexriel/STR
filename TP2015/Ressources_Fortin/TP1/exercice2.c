# include <stdio.h>
# include <string.h>
# include <mpi.h>
# include <unistd.h>

# define SIZE_H_N 50

int main(int argc, char * argv[]){
  int my_rank;              /* rang du processus    */
  int p;                    /* nombre de processus  */
  int source;               /* rang de l'émetteur   */
  int dest;                 /* rang du récepteur    */
  int tag = 0;              /* étiquette du message */
  char message[100];
  MPI_Status status;
  char hostname[SIZE_H_N];

  gethostname(hostname,SIZE_H_N);
  
  /* Initialisation */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  printf("Processus[%d] vous salue.\n", my_rank);

  if (my_rank != 0){
    /* Creation du message */
    sprintf(message, "Coucou du processus #%d depuis %s! :p\n",
	    my_rank, hostname);
    dest = 0;
    printf("Processus[%d] prêt à l'envoi ! ;-)\n", my_rank);
    MPI_Send(message, strlen(message)+1, MPI_CHAR,
	     dest, tag, MPI_COMM_WORLD);
    printf("Processus[%d] a tout lâché ! :D\n", my_rank);
  }
  else{
    for (source = 1; source < p; source++){
      printf("Processus[0] en attente de %d.\n", source);
      MPI_Recv(message, 100, MPI_CHAR, MPI_ANY_SOURCE, tag,
	       MPI_COMM_WORLD, &status);
      printf("Sur %s, le processus #%d a recu le message : %s\n",
	     hostname, my_rank, message);
    }
  }

  /* Desactivation */
  MPI_Finalize();
}
