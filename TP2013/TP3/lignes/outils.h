# ifndef _OUTILS_H_
# define _OUTILS_H_

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

/* Fonction créant une matrice stochastique aléatoire mat */
void random_stochastic_matrix(double *mat, int lignes, int colonnes);


/* Fonction créant un vecteur aléatoire vect */
void random_vector(double *vect, int lignes);


/* Fonction d'affichage de matrice */
void display_matrix(double *mat, int lignes, int colonnes);


/* Fonction d'affichage d'un vecteur en colonnes */
void display_vector(double *vect, double lignes);


/* Fonction d'affichage d'un vecteur en lignes */
void display_vector_line(double *vect, double lignes);


/* Fonction effectuant un produit scalaire de deux vecteurs,
   ici mat est en réalité une ligne de mat */
double scalar_product(double *mat, double *vect, int n);


/* Fonction donnant la norme euclidienne de vect */
double euclidean_norm(double *vect, int n);


/* Fonction normalisant un vecteur vect selon la norme euclidienne */
void normalize_vector(double *vect, int n);


/* Fonction calculant la norme de la différence de deux vecteur v1 et v2 */
double difference_norm(double *v1, double *v2, int n);


/* Fonction copiant un tableau src dans un tableau dest */
void copy_array(double *dest, double *src, int n);


/* Fonction calculant le produit matrice-vecteur entre mat et vect, le stocke dans res */
void produit_matrice_vecteur(double *res, double *mat, double *vect, int lignes, int colonnes);


/* Fonction créant une matrice stochastique aléatoire M et sa transposée tM */
void random_stochastic_matrix_and_its_transpose(double *M, double *tM, int lignes, int colonnes);


/* Fonction multipliant un vecteur v par un double a */
void real_vector_product(double *res, double a, double *v, int n);


# endif
