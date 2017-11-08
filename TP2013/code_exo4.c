#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
int source, cible, nb, recu, i, jeton_recu, tag, tag_recu, new_tag;
int nb_proc = 6;
int jeton=5;
double debut, fin;
MPI_Status status;

MPI_Init(&argc,&argv);
debut=MPI_Wtime();
MPI_Comm_rank(MPI_COMM_WORLD, &source);
MPI_Comm_size(MPI_COMM_WORLD, &nb);

/*source et le numero du processus qui envoi*/
if (source == 0)
{
		 
/*cible est le num du processus a qui on envoi*/
	cible=1;
	MPI_Send(&jeton, 1, MPI_INT, cible, 1, MPI_COMM_WORLD);
	/*on met un tag similaire à chaque couple d'envoi et de reception*/
	tag_recu=1;
	printf("Processus[%d] : OK\n", source);
	MPI_Recv(&jeton_recu, 1, MPI_INT, nb_proc-1, tag_recu, MPI_COMM_WORLD, &status);
	fin=MPI_Wtime();
	printf("Le numéro %d a envoye a %d le jeton %d. Temps d'execution = %f \n", source, cible, jeton, fin-debut);
	MPI_Finalize();
}

else if (source > 0 && source < nb_proc-1) 
{
	cible=(source+1) % nb_proc;
	tag=1;
	new_tag=1;
	MPI_Recv(&jeton_recu, 1, MPI_INT, source-1, tag, MPI_COMM_WORLD, &status);
	printf("Processus[%d] : OK\n", source);
	MPI_Send(&jeton_recu, 1, MPI_INT, cible, new_tag, MPI_COMM_WORLD);
	fin=MPI_Wtime();
	printf("Le numéro %d a envoye a %d le jeton %d. Temps d'execution = %f \n", source, cible, jeton, fin-debut);
	MPI_Finalize();
}
else if (source == (nb_proc -1))
{
	cible=0;
	tag=1;
	new_tag=1;
	MPI_Recv(&jeton_recu, 1, MPI_INT, source-1, tag, MPI_COMM_WORLD, &status);
	printf("Processus[%d] : OK\n", source);
	MPI_Send(&jeton_recu, 1, MPI_INT, cible, new_tag, MPI_COMM_WORLD);
	fin=MPI_Wtime();
	printf("Le numéro %d a envoye a %d le jeton %d. Temps d'execution = %f \n", source, cible, jeton, fin-debut);
	MPI_Finalize();
}


MPI_Finalize();
return 0;
}
