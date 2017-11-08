# include <mpi.h>
# include "outils.h"
# define N 16
const int faux = 0;
int taille = 1;

int main(int argc, char *argv[])
{
  int Nprocs, rang, pseudo_rangL, pseudo_rangC, i, lignes, colonnes, maitre, loc_dim1, loc_dim2;
  double a, debut, fin, norme, epsilon, test, xi, norme_carre, loc_somme;
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
  epsilon  = (argc < 3) ? 0.001 : atof(argv[2]);
  loc_dim1 = lignes/dims[0];
  loc_dim2 = lignes/dims[1];
  maitre = 0;

  srand(time(NULL));

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
  printf("Processus[%d] : loc_x = %f %f %f %f\n", rang, loc_x[0], loc_x[1], loc_x[2], loc_x[3]);





  /* Reste à gérer le problème de la norme, les processus calculent la somme des carrés de leur 
     loc_x, puis communiquent avec les processus de leur communicateur commBlocsLignes leur somme
     partielle pour la sommer, et ensuite tout le monde récupère la somme totale : la norme carrée. */
  loc_somme = somme_carres(loc_x, loc_dim2);
  MPI_Reduce_scatter(&loc_somme, &norme_carre, &taille, MPI_DOUBLE, MPI_SUM, commBlocsLignes);
  norme = sqrt(norme_carre);
  //printf("Processus[%d] : loc_somme = %f, norme = %f\n", rang, loc_somme, norme);
  
  /* On norme x, donc on fait ce qu'il faut à loc_x */
  pseudo_normalize_vector(loc_x, loc_dim2, norme);

  /* On crée loc_y comme copie de loc_x */
  copy_array(loc_y, loc_x, loc_dim2);
  test = epsilon; // + 1; /* pour être certain d'effectuer la boucle au-moins une fois */
  i = 0;

  while (test > epsilon)
    {
      /* On copie loc_x dans loc_y pour pouvoir comparer ||y-x|| avec epsilon à la boucle suivante (test) */
      copy_array(loc_y, loc_x, loc_dim2);
      i++;
      
      if (rang == maitre)
	{
	  printf("Etape %d :\n---------\n", i);
	  printf("y a été copié dans x.\n");
	}

      /* Chaque processus calcule son produit loc_res = loc_blocA * loc_y */
      produit_matrice_vecteur(loc_res, loc_blocA, loc_y, loc_dim1, loc_dim2);
      
      /* Chaque maitre local de chaque communicateur commBlocsLignes (lignes de blocs) doit transmettre
	 son loc_res aux autres membres de commBlocsLignes pour l'additionner, les processus de chaque 
	 communicateur commBlocsLigne disposera du coordonnees[0]-ième morceau de y = A*x. */
      MPI_Allreduce(loc_res, loc_y, loc_dim2, MPI_DOUBLE, MPI_SUM, commBlocsLignes);

      loc_somme = somme_carres(loc_y, loc_dim2);
      MPI_Reduce_scatter(&loc_somme, &norme_carre, &taille, MPI_DOUBLE, MPI_SUM, commBlocsLignes);
      norme = sqrt(norme_carre);
      pseudo_normalize_vector(loc_y, loc_dim2, norme);
      
      /* Il faut calculer ||y-x|| */
      difference_vectors(loc_res, loc_y, loc_x, loc_dim2);
      loc_somme = somme_carres(loc_res, loc_dim2);
      MPI_Reduce_scatter(&loc_somme, &norme_carre, &taille, MPI_DOUBLE, MPI_SUM, commBlocsLignes);
      test = sqrt(norme_carre);
      
    }


  /* On peut pour vérifier, le résultat, rassembler les loc_y dans y dans chaque communicateur
     commBlocsColonnes */
  MPI_Gather(loc_y, loc_dim2, MPI_DOUBLE, y, loc_dim2, MPI_DOUBLE, maitre, commBlocsColonnes);
  

  if (rang == maitre)
    {
      printf("A*...*A*x :\n");
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
