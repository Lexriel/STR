# include <mpi.h>
# include "outils.h"
# define N 16

int main(int argc, char *argv[])
{
  int lignes, colonnes, i;
  int rang, Nprocs, maitre;
  double a, debut, fin, norm, epsilon, test, xi;
  double *A, *tA, *x, *y, *loc_colonneA, *res;
  
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

  /* Allocation de x et colonneA commune à tous les processus */
  x   = (double *) malloc( lignes * sizeof(double) );
  res = (double *) malloc( lignes * sizeof(double) );
  loc_colonneA = (double *) malloc( colonnes * sizeof(double) );

  if (rang == 0)
    {
      debut = MPI_Wtime();

      /* seul le processus de rang maitre travaille sur A */
      A   = (double *) malloc( lignes * colonnes * sizeof(double) );
      tA  = (double *) malloc( colonnes * lignes * sizeof(double) );
      y   = (double *) malloc( lignes * sizeof(double) );

      random_stochastic_matrix_and_its_transpose(A, tA, lignes, colonnes);
      printf("Le maitre a calculé la matrice A :\n");
      display_matrix(A, lignes, colonnes);
      /* En fait on peut ne calculer que tA, mais bon, c'est un choix... */
      printf("Le maitre a calculé la matrice tA :\n");
      display_matrix(tA, colonnes, lignes);

      random_vector(x, lignes);
      printf("Le maitre a calculé le vecteur x.\n\n");
      /* display_vector(x, lignes); */
    }
  
  /* Le maitre envoie une ligne tA[rang] (soit une colonne de A) sur loc_colonneA pour chaque rang */
  MPI_Scatter(tA, lignes, MPI_DOUBLE, loc_colonneA, lignes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Scatter : le maitre a envoyé A par colonnes aux processus.\n");

  /* controle OK
     for (i=0; i<Nprocs; i++) if (rang==i) display_matrix(colonneA, 1, colonnes);  */

  /* Le maitre envoie chaque coefficient de x à un processus différent */
  MPI_Scatter(x, 1, MPI_DOUBLE, &xi, 1, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Scatter : le maitre a envoyé x par morceaux aux processus.\n");

  /* Chaque processus calcule son produit xi*loc_colonneA */
  real_vector_product(res, xi, loc_colonneA, lignes);
  if (rang == maitre)
    printf("Chaque processus a calculé son produit res=xi*loc_colonneA.\n");

  /* controle OK
     printf("res = %.3f \n", res); */

  /* Tout le monde envoie son loc_colonneA à maitre pour les additionner, la somme sera le produit A*x */
  MPI_Reduce(res, y, colonnes, MPI_DOUBLE, MPI_SUM, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Reduce  : le maitre a reçu tous les loc_colonnes envoyés par les processus et les somme pour obtenir A*x.\n\n");


  if (rang == maitre)
    {
      printf("A*x :\n");
      display_vector_line(y, lignes);
      free(y);
      free(A);
      free(tA);

      fin = MPI_Wtime();
      printf("Temps d'exécution du processus maitre : %f s.\n\n", fin-debut);
    }


  free(loc_colonneA);
  free(x);
  free(res);
  
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
