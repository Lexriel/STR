# include "outils.h"

/* Fonction créant une matrice stochastique aléatoire mat */
void random_stochastic_matrix(double *mat, int lignes, int colonnes)
{
  int i, j, indice;
  double somme;

  for (i=0; i<lignes; i++)
    {
      somme = 0;
      for (j=0; j<colonnes; j++)
	{
	  indice = i*lignes + j;
	  mat[indice] = (double) (rand() % 10000);
	  somme += mat[indice];
	}
      for (j=0; j<colonnes; j++)
	mat[i*lignes + j] /= somme;
    }
}


/* Fonction créant un vecteur aléatoire vect */
void random_vector(double *vect, int lignes)
{
  int i;

  for (i=0; i<lignes; i++)
    *vect++ = (rand() % 1000) / 100.;
}


/* Fonction d'affichage de matrice */
void display_matrix(double *mat, int lignes, int colonnes)
{
  int i, j;

  for (i=0; i<lignes; i++)
    {
      for (j=0; j<colonnes; j++)
	printf("%.4f ", *mat++);
      printf("\n");
    }
  printf("\n");
}


/* Fonction d'affichage d'un vecteur en colonnes */
void display_vector(double *vect, double lignes)
{
  display_matrix(vect, lignes, 1);
}


/* Fonction d'affichage d'un vecteur en lignes */
void display_vector_line(double *vect, double lignes)
{
  display_matrix(vect, 1, lignes);
}


/* Fonction effectuant un produit scalaire de deux vecteurs,
   ici mat est en réalité une ligne de mat */
double scalar_product(double *mat, double *vect, int n)
{
  int i;
  double res;

  res = 0;

  for (i=0; i<n; i++)
    res += (*mat++) * (*vect++);

  return res;
}


/* Fonction donnant la norme euclidienne de vect */
double euclidean_norm(double *vect, int n)
{
  int i;
  double res;
  res = 0;
  
  for (i=0; i<n; i++)
    res += (*vect) * (*vect++);

  return sqrt(res);
}


/* Fonction normalisant un vecteur vect selon la norme euclidienne */
void normalize_vector(double *vect, int n)
{
  int i;
  double norm;
  norm = euclidean_norm(vect, n);

  for (i=0; i<n; i++)
    *vect++ /= norm;
}


/* Fonction calculant la norme de la différence de deux vecteur v1 et v2 */
double difference_norm(double *v1, double *v2, int n)
{
  int i;
  double res, tmp;
  res = 0;
  
  for (i=0; i<n; i++)
    {
      tmp = (*v1) * (*v1++) - (*v2) * (*v2++);
      tmp = fabs(tmp);
      res += tmp;
    }

  return sqrt(res);
}


/* Fonction copiant un tableau src dans un tableau dest */
void copy_array(double *dest, double *src, int n)
{
  int i;

  for (i=0; i<n; i++)
    *dest++ = *src++;
}


/* Fonction calculant le produit matrice-vecteur entre mat et vect, le stocke dans res */
void produit_matrice_vecteur(double *res, double *mat, double *vect, int lignes, int colonnes)
{
  double somme;
  int i, j, indice;

  for (j=0; j<colonnes; j++)
    {
      somme = 0;
      for (i=0; i<lignes; i++)
	somme += (*mat++) * vect[i];
      (*res++) = somme;
    }
}


/* Fonction créant une matrice stochastique aléatoire M et sa transposée tM */
void random_stochastic_matrix_and_its_transpose(double *M, double *tM, int lignes, int colonnes)
{
  int i, j, indice;
  double somme;

  for (i=0; i<lignes; i++)
    {
      somme = 0;
      for (j=0; j<colonnes; j++)
	{
	  indice = i*lignes + j;
	  M[indice] = (double) (rand() % 10000);
	  somme += M[indice];
	}
      for (j=0; j<colonnes; j++)
	{
	  indice = i*lignes + j;
	  M[indice] /= somme;
	  tM[j*colonnes + i] = M[indice];
        }
    }
}

/* Fonction multipliant un vecteur v par un double a et affectant le résultat dans res */
void real_vector_product(double *res, double a, double *v, int n)
{
  int i;

  for (i=0; i<n; i++)
    *res++ = (*v++) * a;
}


/* Fonction créant une sous-matrice d'une matrice stochastique aléatoire mat */
void random_stochastic_submatrix(double *mat, int lignes, int colonnes, int proportion)
{
  int i, j, indice;
  double somme;

  for (i=0; i<lignes; i++)
    {
      somme = 0;
      for (j=0; j<colonnes; j++)
	{
	  indice = i*lignes + j;
	  mat[indice] = (double) (rand() % 10000);
	  somme += mat[indice];
	}
      for (j=0; j<colonnes; j++)
	mat[i*lignes + j] /= (proportion*somme);
    }
}


/* Fonction calculant la somme des carrés d'un vecteur v */
double somme_carres(double *v, int n)
{
  int i, somme;
  somme = 0;

  for (i=0; i<n; i++)
    somme += (*v) * (*v++);

  return somme;
}


/* Fonction normalisant juste une partie d'un vecteur sans connaitre le vecteur entier ni sa norme */
void pseudo_normalize_vector(double *vect, int n, double norme)
{
  int i;

  for (i=0; i<n; i++)
    *vect++ /= norme;
}


/* Fonction calculant la différence en valeur absolue de deux vecteurs v1 et v2 */
void difference_vectors(double *res, double *v1, double *v2, int n)
{
  int i;
  double tmp;
  res = 0;
  
  for (i=0; i<n; i++)
    {
      tmp = (*v1++) - (*v2++);
      tmp = fabs(tmp);
      *res++ = tmp;
    }
}
