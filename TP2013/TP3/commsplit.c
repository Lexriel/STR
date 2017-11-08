/*
 * commsplit.f90 --- Subdiviser une grille 2D avec MPI_COMM_SPLIT
 *
 *
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
const int faux = 0;

int main(int argc, char *argv[]) {
  int Nprocs, rang, i, maitre = 0;
  double v[4];
  double w;

  /* Variables pour créer de nouveaux communicateurs : carte2D, carte1D */
  MPI_Comm carte2D, carte1D;
  const int Ndims = 2;
  int coordonnees[Ndims];
  int dims[Ndims];
  dims[0] = dims[1] = 4;

  MPI_Init( &argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
  if (dims[0]*dims[1] != Nprocs)
    {
      fprintf(stderr, "Ce n'est pas le nombre de processus necessaire !");
      exit(EXIT_FAILURE);
    }

  /* Déclarations peu importantes pour la création de carte2D */
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

  /* Initialisation du vecteur V et du scalaire w dans la "colonne" maitre */
  w = v[0] = v[1] = v[2] = v[3] = 0;
  if (coordonnees[0] == maitre){
    v[0] = rang; v[1] = rang; v[2] = rang; v[3] = rang; w = rang;}

  /* Subdivision de la grille 2D a l'aide de MPI_COMM_SPLIT */
  /* carte1D représente des communicateurs identifiés par leur "colonne" coordonnees[1] 
     et les processus sont identifiés par leur "ligne" coordonnees[0] en tant que 
     pseudo-rang dans ce communicateur */
  MPI_Comm_split(carte2D, coordonnees[1], coordonnees[0], &carte1D);

  /* Les processus de la ligne 2 diffusent seélectivement
   * le vecteur V aux processus de leur ligne */
  MPI_Scatter(v, 1, MPI_DOUBLE, &w, 1, MPI_DOUBLE, maitre, carte1D);

  printf("Processus[%2d] : coordonnées (%d,%d) , w = %.2f\n", rang, coordonnees[0], coordonnees[1], w);

  MPI_Comm_free(&carte1D);
  MPI_Comm_free(&carte2D);

  MPI_Finalize();
  return 0;
}
