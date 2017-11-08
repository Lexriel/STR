# include <mpi.h>
# include "outils.h"
# define N 16
const int faux = 0;

int main(int argc, char *argv[])
{
  int Nprocs, rang, pseudo_rangL, pseudo_rangC, i, lignes, colonnes, maitre, loc_dim1, loc_dim2;
  double a, debut, fin, norm, epsilon, test, xi;
  double *loc_x, *loc_y, *loc_blocA, *loc_res, *y;

  /* Variables pour créer de nouveaux communicateurs : carte2D, carteColonnesBlocs */
  MPI_Comm carte2D, commBlocsColonnes, commBlocsLignes;
  const int Ndims = 2;
  int coordonnees[Ndims];
  int dims[Ndims];
  dims[0] = dims[1] = 4;

  MPI_Init( &argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
  if (dims[0]*dims[1] != Nprocs)
    {
      fprintf(stderr, "Ce n'est pas le nombre de processus nécessaire !\n");
      exit(EXIT_FAILURE);
    }

  debut = MPI_Wtime();
  lignes   = (argc < 2) ? N : atoi(argv[1]);
  colonnes = (argc < 2) ? N : atoi(argv[1]);
  loc_dim1 = lignes/dims[0];
  loc_dim2 = lignes/dims[1];
  maitre = 0;

  /* Déclarations pour la création de carte2D */
  /* le rang de carte2D sera le même que dans MPI_COMM_WORLD grâce à reordonne = faux */
  int reordonne = faux;
  /* spécifie si la grille est périodique, non nécessaire ici */
  int periode[Ndims];
  for (i=0; i<Ndims; i++)
    periode[i] = faux;

  /* Création du communicateur commcart2D : MPI_COMM_WORLD avec une topologie cartésienne 2D */
  MPI_Cart_create(MPI_COMM_WORLD, Ndims, dims, periode, reordonne, &carte2D);
  MPI_Comm_rank(carte2D, &rang); 
  MPI_Cart_coords(carte2D, rang, Ndims, coordonnees);

  /* Subdivision de la grille 2D a l'aide de MPI_COMM_SPLIT */
  /* commBlocsColonnes représente des communicateurs identifiés par leur "colonne" coordonnees[1] 
     et les processus sont identifiés par leur "ligne" coordonnees[0] en tant que 
     pseudo-rang dans ce communicateur */
  MPI_Comm_split(carte2D, coordonnees[1], coordonnees[0], &commBlocsColonnes);
  /* Inversement, on crée*/
  MPI_Comm_split(carte2D, coordonnees[0], coordonnees[1], &commBlocsLignes);
  pseudo_rangL = coordonnees[1];
  pseudo_rangC = coordonnees[0];

  /* Allocation de loc_x et loc_blocA commune à tous les processus */
  loc_x     = (double *) malloc( loc_dim2 * sizeof(double) );
  loc_y     = (double *) malloc( colonnes * sizeof(double) );
  y         = (double *) malloc( colonnes * sizeof(double) );
  loc_res   = (double *) malloc( loc_dim2 * sizeof(double) );
  loc_blocA = (double *) malloc( loc_dim1 * loc_dim2 * sizeof(double) );

  /* Chaque processus crée un morceau d'une matrice A qui sera stochastique */
  random_stochastic_submatrix(loc_blocA, loc_dim1, loc_dim2, loc_dim2);
  if (rang == maitre)
    printf("Chaque processus a calculé son sous-bloc de A.\n");

  /* Chaque maitre local de chaque communicateur commBlocsColonnes (colonnes de blocs) construit
     un morceau de x qu'il doit transmettre aux autres membres de commBlocsColonnes */
  if (coordonnees[0] == maitre)
    random_vector(loc_x, loc_dim2);
  MPI_Bcast(loc_x, loc_dim2, MPI_DOUBLE, maitre, commBlocsColonnes);

  /* Chaque processus calcule son produit loc_res = loc_blocA * loc_x */
  produit_matrice_vecteur(loc_res, loc_blocA, loc_x, loc_dim1, loc_dim2);

  /* Chaque maitre local de chaque communicateur commBlocsLignes (lignes de blocs) doit transmettre
     son loc_res aux autres membres de commBlocsLignes pour l'additionner, les processus de chaque 
     communicateur commBlocsLigne disposera du coordonnees[0]-ième morceau de y = A*x. */
  MPI_Allreduce(loc_res, loc_y, loc_dim2, MPI_DOUBLE, MPI_SUM, commBlocsLignes);


  /* On peut pour vérifier, le résultat, rassembler les loc_y dans y dans chaque communicateur
     commBlocsColonnes */
  MPI_Gather(loc_y, loc_dim2, MPI_DOUBLE, y, loc_dim2, MPI_DOUBLE, maitre, commBlocsColonnes);

  if (rang == maitre)
    {
      printf("A*x :\n");
      display_vector_line(y, lignes);
      fin = MPI_Wtime();
      printf("Temps d'exécution : %f s.\n\n", fin-debut);
    }

  MPI_Comm_free(&commBlocsColonnes);
  MPI_Comm_free(&commBlocsLignes);
  MPI_Comm_free(&carte2D);

  free(loc_x);
  free(loc_y);
  free(loc_res);
  free(y);
  free(loc_blocA);

  MPI_Finalize();
  return 0;
}
