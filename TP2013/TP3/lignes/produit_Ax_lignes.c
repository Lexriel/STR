# include <mpi.h>
# include "outils.h"
# define N 16

int main(int argc, char *argv[])
{
  int lignes, colonnes, i;
  int rang, Nprocs, maitre;
  double a, debut, fin;
  double *A, *x, *res, *loc_ligneA;
  
  /* Initialisation de MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);
  MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
  maitre = 0;

  lignes   = (argc < 2) ? N : atoi(argv[1]);
  colonnes = (argc < 2) ? N : atoi(argv[1]);

  if (lignes != Nprocs)
    {
      printf("ERREUR : le nombre de lignes doit être égal au nombre de processus.\n");
      exit(EXIT_FAILURE);
    }

  /* Allocation de x et ligneA commune à tous les processus */
  x = (double *) malloc( lignes * sizeof(double) );
  loc_ligneA = (double *) malloc( colonnes * sizeof(double) );

  if (rang == 0)
    {
      debut = MPI_Wtime();

      /* seul le processus de rang maitre travaille sur A */
      A = (double *) malloc( lignes * colonnes * sizeof(double) );
      res = (double *) malloc( lignes * sizeof(double) );

      random_stochastic_matrix(A, lignes, colonnes);
      printf("Le maitre a calculé la matrice A :\n");
      display_matrix(A, lignes, colonnes);

      random_vector(x, lignes);
      printf("Le maitre a calculé le vecteur x.\n\n");
      /* display_vector(x, lignes); */
    }
  
  /* Le maitre envoie une ligne A[rang] sur ligneA pour chaque rang */
  MPI_Scatter(A, colonnes, MPI_DOUBLE, loc_ligneA, colonnes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Scatter : le maitre a envoyé A par lignes aux processus.\n");

  /* controle OK
     for (i=0; i<Nprocs; i++) if (rang==i) display_matrix(ligneA, 1, colonnes);  */

  /* Le maitre envoie x à tout le monde */
  MPI_Bcast(x, lignes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Bcast   : le maitre a envoyé x aux processus.\n");

  /* Chaque processus calcul le produit scalaire de sa ligne avec x */
  a = scalar_product(loc_ligneA, x, colonnes);
  if (rang == maitre)
    printf("Chaque processus a calculé un coefficient de res.\n");

  /* controle OK
     printf("res = %.3f \n", res); */

  /* Tout le monde envoie son a dans res au maitre */
  MPI_Gather(&a, 1, MPI_DOUBLE, res, 1, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Gather  : le maitre a reçu tous les coefficients calculés par les processus.\n\n");


  if (rang == maitre)
    {
      printf("A*x :\n");
      display_vector_line(res, lignes);
      free(res);
      free(A);
      
      fin = MPI_Wtime();
      printf("Temps d'exécution du processus maitre : %f s.\n\n", fin-debut);
    }


  free(loc_ligneA);
  free(x);
  
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
