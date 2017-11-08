#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define MSG_SIZE 100

int main(int argc, char* argv[])
{
    int         my_rank;      /* rang du processus     */
    int         p;            /* nombre de  processus  */
    int         source;       /* rang de l'emetteur */
    int         dest;         /* rang du recepteur     */
    int         tag = 0;      /* etiquette du  message     */
    char        message[MSG_SIZE];  
    MPI_Status  status;        
                               
    /* Initialisation 
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    

    /* Creation du message */
    sprintf(message, "Coucou du processus  %d!",  my_rank);
    dest = (my_rank == p-1 ? 0 : my_rank + 1);

    MPI_Send(message, strlen(message)+1, MPI_CHAR,
	     dest, tag, MPI_COMM_WORLD);

    /* Reception du message */
    source = (my_rank == 0 ? p-1 : my_rank - 1);
    MPI_Recv(message, MSG_SIZE, MPI_CHAR, source, tag,
	     MPI_COMM_WORLD, &status);
    printf("Je suis le processus %d et j'ai recu ce message du processus %d : %s\n", 
	   my_rank, source, message);

    /* Desactivation
     */
    MPI_Finalize();
} 

