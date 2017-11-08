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
  epsilon  = (argc < 3) ? 0.001 : atof(argv[2]);

  if (lignes != Nprocs)
    {
      printf("ERREUR : le nombre de lignes doit être égal au nombre de processus.\n");
      exit(EXIT_FAILURE);
    }

  /* Allocation de x, y et colonneA commune à tous les processus */
  x   = (double *) malloc( lignes * sizeof(double) );
  y   = (double *) calloc( lignes, sizeof(double) );
  res =  (double *) malloc( lignes * sizeof(double) );
  loc_colonneA = (double *) malloc( colonnes * sizeof(double) );

  if (rang == maitre)
    {
      debut = MPI_Wtime();

      /* seul le processus de rang maitre travaille sur A et tA */
      A  = (double *) malloc( lignes * colonnes * sizeof(double) );
      tA = (double *) malloc( colonnes * lignes * sizeof(double) );

      random_stochastic_matrix_and_its_transpose(A, tA, lignes, colonnes);
      printf("Le maitre a calculé la matrice A :\n");
      display_matrix(A, lignes, colonnes);
      /* En fait on peut ne calculer que tA, mais bon, c'est un choix... */
      printf("Le maitre a calculé la matrice tA :\n");
      display_matrix(tA, colonnes, lignes);

      random_vector(y, lignes);
      normalize_vector(y, lignes);
      printf("Le maitre a calculé le vecteur x (nommé y pour l'instant) :\n");
      display_vector_line(y, lignes);
    }
  
  /* Le maitre envoie une ligne tA[rang] (soit une colonne de A) sur loc_colonneA pour chaque rang */
  MPI_Scatter(tA, lignes, MPI_DOUBLE, loc_colonneA, lignes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Scatter : le maitre a envoyé A par colonnes aux processus.\n");


  /* Le maitre envoie y sur chaque processus, seul intérêt : éviter de renvoyer
     les morceaux xi à maitre à chaque étape pour qu'il calcule test */
  MPI_Bcast(y, lignes, MPI_DOUBLE, maitre, MPI_COMM_WORLD);
  if (rang == maitre)
    printf("MPI_Bcast   : le maitre a envoyé y aux processus.\n\n");


  /* ============================================================

              BOUCLE effectuant A*A*...*A*x itérativement

     ============================================================ */

  test = epsilon + 1; /* pour être certain d'effectuer la boucle au-moins une fois */
  i = 0;

  while (test > epsilon)
    {
      /* On copie y dans x pour pouvoir comparer ||y-x|| avec epsilon à la boucle suivante (test) */
      copy_array(x, y, lignes);
      xi = x[rang];
      i++;

      if (rang == maitre)
	{
	  printf("Etape %d :\n---------\n", i);
	  printf("y a été copié dans x.\n");
	}

      /* Chaque processus calcule son produit xi*loc_colonneA */
      real_vector_product(res, xi, loc_colonneA, lignes);
      
      if (rang == maitre)
	printf("Chaque processus a calculé son produit res=xi*loc_colonneA.\n");


      /* Tout le monde envoie son loc_colonneA à tout le monde pour les additionner, la somme sera le produit A*x */
      MPI_Allreduce(res, y, lignes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (rang == maitre)
	printf("MPI_AllReduce : tout le monde a reçu tous les res envoyés par les processus et les somme pour obtenir y=A*x.\n");

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
      free(tA);
    }


  free(loc_colonneA);
  free(x);
  free(y);
  free(res);
  
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
