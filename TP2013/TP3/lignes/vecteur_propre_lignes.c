# include <mpi.h>
# include "outils.h"
# define N 16

int main(int argc, char *argv[])
{
  int lignes, colonnes, i;
  int rang, Nprocs, maitre;
  double a, debut, fin, norm, epsilon, test;
  double *A, *x, *y, *loc_ligneA;


  /* Initialisation de MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
  maitre = 0;

  lignes   = (argc < 2) ? N : atoi(argv[1]);
  colonnes = (argc < 2) ? N : atoi(argv[1]);
  epsilon  = (argc < 3) ? 0.001 : atof(argv[2]);

  if (lignes != Nprocs)
    {
      printf("ERREUR : le nombre de lignes doit être égal au nombre de processus.\n");
      exit(EXIT_FAILURE);
    }

  /* Allocation de x, y et ligneA commune à tous les processus */
  x = (double *) malloc( lignes * sizeof(double) );
  loc_ligneA = (double *) malloc( colonnes * sizeof(double) );
  y = (double *) calloc( lignes, sizeof(double) );

  if (rang == maitre)
    {
      debut = MPI_Wtime();

      /* seul le processus de rang maitre travaille sur A */
      A = (double *) malloc( lignes * colonnes * sizeof(double) );

      random_stochastic_matrix(A, lignes, colonnes);
      printf("Le maitre a calculé la matrice A :\n");
      display_matrix(A, lignes, colonnes);

      random_vector(y, lignes);
      normalize_vector(y, lignes);
      printf("Le maitre a calculé le vecteur x (nommé y pour l'instant) :\n\n");
      display_vector_line(y, lignes);
    }
  
  /* Le maitre envoie une ligne A[rang] sur ligneA pour chaque rang */
  MPI_Scatter(A, colonnes, MPI_DOUBLE, loc_ligneA, colonnes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Scatter   : le maitre a envoyé A par lignes aux processus.\n");

  /* Le maitre envoie x (nommé y pour l'instant) à tout le monde */
  MPI_Bcast(y, lignes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Bcast     : le maitre a envoyé y aux processus.\n\n");


  /* ============================================================

              BOUCLE effectuant A*A*...*A*x itérativement

     ============================================================ */

  test = epsilon + 1; /* pour être certain d'effectuer la boucle au-moins une fois */
  i = 0;

  while (test > epsilon)
    {
      /* On copie y dans x pour pouvoir comparer ||y-x|| avec epsilon à la boucle suivante (test) */
      copy_array(x, y, lignes);
      i++;

      if (rang == maitre)
	{
	  printf("Etape %d :\n---------\n", i);
	  printf("y a été copié dans x.\n");
	}
      
      /* Chaque processus calcule le produit scalaire de sa ligne avec y */
      a = scalar_product(loc_ligneA, y, colonnes);
      if (rang == maitre)
	printf("Chaque processus a calculé un coefficient de A*y.\n");

      /* Tous les processus s'envoient leur a et le stockent dans leur y */
      /* MPI_Allgather remplace un MPI_Gather suivi d'un MPI_Bcast. */
      MPI_Allgather(&a, 1, MPI_DOUBLE, y, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      if (rang == maitre)
	printf("MPI_Allgather : les processus ont tous rempli y avec les coefficients calculés par chaque processus.\n");

      normalize_vector(y, lignes);
      if (rang == maitre)
	printf("Chaque processus a normé y.\n");
	

      test = difference_norm(y, x, lignes);
      if (rang == maitre)
	printf("\n");
    }


  if (rang == maitre)
    {
      fin = MPI_Wtime();
      printf("\n");
      printf("La suite A*A*...*A*x converge vers :\n");
      display_vector_line(y, lignes);

      printf("Temps d'exécution du processus maitre : %f s.\n\n", fin-debut);
      free(A);
    }


  free(loc_ligneA);
  free(x);
  free(y);
  
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}