#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>

#define SIZE_H_N 50

int main(int argc, char* argv[])
{
    int         my_rank;      /* rang du processus     */
    int         p;            /* nombre de  processus  */
    int         source;       /* rang de l'emetteur */
    int         dest;         /* rang du recepteur     */
    int         tag = 0;      /* etiquette du  message     */
    char        message[100];  
    MPI_Status  status;        
    char hostname[SIZE_H_N] ;              
    
    gethostname(hostname,SIZE_H_N);  

    /* Initialisation  */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);


    if (my_rank != 0) 
    {
        /* Creation du message */
        sprintf(message, "Coucou du processus #%d depuis %s!",  
		my_rank, hostname);
        dest = 0;

        MPI_Send(message, strlen(message)+1, MPI_CHAR,
            dest, tag, MPI_COMM_WORLD);
    } else
    { 
        for (source = 1; source < p; source++) {
            MPI_Recv(message, 100, MPI_CHAR, source, tag,
                MPI_COMM_WORLD, &status);
            printf("Sur %s, le processus #%d a recu le message : %s\n", 
		   hostname, my_rank, message);
        }
    }

    /* Desactivation */
    MPI_Finalize();
} 

