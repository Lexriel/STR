/*
 * Universit� Pierre et Marie Curie
 * Calcul de convolution sur une image.
 *
 * Compilation, ex�cution : voir fichier Makefile
 * Version parallele avec Scatter
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>   /* pour le rint */
#include <string.h> /* pour le memcpy */

#include <time.h>   /* chronometrage */
#include <mpi.h>

#include "rasterfile.h"

#define MAX(a,b) ((a>b) ? a : b)


/** 
 * \struct Raster
 * Structure d�crivant une image au format Sun Raster
 */

typedef struct {
  struct rasterfile file;  ///< Ent�te image Sun Raster
  unsigned char rouge[256],vert[256],bleu[256];  ///< Palette de couleur
  unsigned char *data;    ///< Pointeur vers l'image
} Raster;


/**
 * Cette procedure convertit un entier LINUX en un entier SUN 
 *
 * \param i pointeur vers l'entier � convertir
 */

void swap(int *i) {
  unsigned char s[4],*n;
  memcpy(s,i,4);
  n=(unsigned char *)i;
  n[0]=s[3];
  n[1]=s[2];
  n[2]=s[1];
  n[3]=s[0];
}

/**
 * \brief Lecture d'une image au format Sun RASTERFILE.
 *
 * Au retour de cette fonction, la structure r est remplie
 * avec les donn�es li�e � l'image. Le champ r.file contient
 * les informations de l'entete de l'image (dimension, codage, etc).
 * Le champ r.data est un pointeur, allou� par la fonction
 * lire_rasterfile() et qui contient l'image. Cette espace doit
 * �tre lib�r� apr�s usage.
 *
 * \param nom nom du fichier image
 * \param r structure Raster qui contient l'image
 *  charg�e en m�moire
 */

void lire_rasterfile(char *nom, Raster *r) {
  FILE *f;
  int i;
    
  if( (f=fopen( nom, "r"))==NULL) {
    fprintf(stderr,"erreur a la lecture du fichier %s\n", nom);
    exit(1);
  }
  fread( &(r->file), sizeof(struct rasterfile), 1, f);    
  swap(&(r->file.ras_magic));
  swap(&(r->file.ras_width));
  swap(&(r->file.ras_height));
  swap(&(r->file.ras_depth));
  swap(&(r->file.ras_length));
  swap(&(r->file.ras_type));
  swap(&(r->file.ras_maptype));
  swap(&(r->file.ras_maplength));
    
  if ((r->file.ras_depth != 8) ||  (r->file.ras_type != RT_STANDARD) ||
      (r->file.ras_maptype != RMT_EQUAL_RGB)) {
    fprintf(stderr,"palette non adaptee\n");
    exit(1);
  }
    
  /* composante de la palette */
  fread(&(r->rouge),r->file.ras_maplength/3,1,f);
  fread(&(r->vert), r->file.ras_maplength/3,1,f);
  fread(&(r->bleu), r->file.ras_maplength/3,1,f);
    
  if ((r->data=malloc(r->file.ras_width*r->file.ras_height))==NULL){
    fprintf(stderr,"erreur allocation memoire\n");
    exit(1);
  }
  fread(r->data,r->file.ras_width*r->file.ras_height,1,f);
  fclose(f);
}

/**
 * Sauve une image au format Sun Rasterfile
 */

void sauve_rasterfile(char *nom, Raster *r)     {
  FILE *f;
  int i;
  
  if( (f=fopen( nom, "w"))==NULL) {
    fprintf(stderr,"erreur a l'ecriture du fichier %s\n", nom);
    exit(1);
  }
    
  swap(&(r->file.ras_magic));
  swap(&(r->file.ras_width));
  swap(&(r->file.ras_height));
  swap(&(r->file.ras_depth));
  swap(&(r->file.ras_length));
  swap(&(r->file.ras_type));
  swap(&(r->file.ras_maptype));
  swap(&(r->file.ras_maplength));
    
  fwrite(&(r->file),sizeof(struct rasterfile),1,f);
  /* composante de la palette */
  fwrite(&(r->rouge),256,1,f);
  fwrite(&(r->vert),256,1,f);
  fwrite(&(r->bleu),256,1,f);
  /* pour le reconvertir pour la taille de l'image */
  swap(&(r->file.ras_width));
  swap(&(r->file.ras_height));
  fwrite(r->data,r->file.ras_width*r->file.ras_height,1,f); 
  fclose(f);
}

/**
 * R�alise une division d'entiers plus pr�cise que
 * l'op�rateur '/' : renvoie le quotient entier le plus proche
 * du quotient r�el. 
 * Remarque : la fonction rint provient de la librairie 
 * math�matique.
 */

unsigned char division(int numerateur,int denominateur) {
  
  if (denominateur != 0)
    return (unsigned char) rint((double)numerateur/(double)denominateur); 
  else 
    return 0;
}

static int ordre( unsigned char *a, unsigned char *b) {
  return (*a-*b);
}


typedef enum {
  CONVOL_MOYENNE1, ///< Filtre moyenneur
  CONVOL_MOYENNE2, ///< Filtre moyenneur central
  CONVOL_CONTOUR1, ///< Laplacien
  CONVOL_CONTOUR2, ///< Max gradient
  CONVOL_MEDIAN    ///< Filtre m�dian
} filtre_t;

/**
 * R�alise une op�ration de convolution avec un noyau pr�d�fini sur
 * un point.
 *
 * \param choix type de noyau pour la convolution :
 *  - CONVOL_MOYENNE1 : filtre moyenneur
 *  - CONVOL_MOYENNE2 : filtre moyenneur avec un poid central plus fort
 *  - CONVOL_CONTOUR1 : filtre extracteur de contours (laplacien)
 *  - CONVOL_CONTOUR2 : filtre extracteur de contours (max des composantes du gradient)
 *  - CONVOL_MEDIAN : filtre m�dian (les 9 valeurs sont tri�es et la valeur
 *     m�diane est retourn�e).
 * \param NO,N,NE,O,CO,E,SO,S,SE: les valeurs des 9 points
 *  concern�s pour le calcul de la convolution (cette derni�re est
 *  formellement une combinaison lin�aire de ces 9 valeurs).
 * \return la valeur de convolution.
 */

unsigned char filtre( filtre_t choix, 
		      unsigned char NO, unsigned char N,unsigned char NE, 
		      unsigned char O,unsigned char CO, unsigned char E, 
		      unsigned char SO,unsigned char S,unsigned char SE) {
  int numerateur,denominateur;
  unsigned char tab[] = {NO,N,NE,O,CO,E,SO,S,SE};

  switch (choix)
    {
    case CONVOL_MOYENNE1:
	  /* filtre moyenneur */
	  numerateur = (int)NO + (int)N + (int)NE + (int)O + (int)CO + 
	    (int)E + (int)SO + (int)S + (int)SE;
	  denominateur = 9;
	  return division(numerateur,denominateur); 

    case CONVOL_MOYENNE2:
	  /* filtre moyenneur */
	  numerateur = (int)NO + (int)N + (int)NE + (int)O + 4*(int)CO +
	    (int)E + (int)SO + (int)S + (int)SE;
	  denominateur = 12;
	  return division(numerateur,denominateur);	

    case CONVOL_CONTOUR1:
	  /* extraction de contours */
	  numerateur = -(int)N - (int)O + 4*(int)CO - (int)E - (int)S;
	  /* numerateur = -(int)NO -(int)N - (int)NE - (int)O + 8*(int)CO -
		 (int)E - (int)SO - (int)S - (int)SE;
	  */
	  return ((4*abs(numerateur) > 255) ? 255 :  4*abs(numerateur));

    case CONVOL_MEDIAN:
	  /* filtre non lineaire : tri rapide sur la brillance */
	  qsort( tab, 9, sizeof(unsigned char), (int (*) (const void *,const void *))ordre);
	  return tab[4];

    case CONVOL_CONTOUR2:
	  /* extraction de contours */
	  numerateur = MAX(abs(CO-E),abs(CO-S));
	  return ((4*numerateur > 255) ? 255 :  4*numerateur);
    default:
	  printf("\nERREUR : Filtre inconnu !\n\n");
	  exit(1);
    }
}

/**
 * Convolution d'une image par un filtre pr�d�fini
 * \param choix choix du filtre (voir la fonction filtre())
 * \param tab pointeur vers l'image
 * \param nbl, nbc dimension de l'image
 *
 * \sa filtre()
 */

int convolution( filtre_t choix, unsigned char tab[],int nbl,int nbc) {
  int i,j;
  unsigned char *tmp;
  
  tmp = (unsigned char*) malloc(sizeof(unsigned char) *nbc*nbl);
  if (tmp == NULL) {
    printf("Erreur dans l'allocation de tmp dans convolution \n");
    return 1;
  }
  /* on laisse tomber les bords */
  for(i=1 ; i<nbl-1 ; i++)
    for(j=1 ; j<nbc-1 ; j++)
      tmp[i*nbc+j] = filtre(
			    choix,
			    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
			    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
			    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
  
  /* Recopie de l'image apres traitement dans l'image initiale,
   * On remarquera que la premiere, la derniere ligne, la premiere
   * et la derniere colonne ne sont copi�es (ce qui force a faire
   * la copie ligne par ligne). */
  for( i=1; i<nbl-1; i++)
    memcpy( tab+nbc*i+1, tmp+nbc*i+1, (nbc-2)*sizeof(unsigned char));

  free(tmp);   
}


/**
 * Interface utilisateur
 */

static char usage [] = "Usage : %s <nom image SunRaster> [0|1|2|3|4] <nbiter>\n";

/*
 * Partie principale
 */

int main(int argc, char *argv[]) {

  /* Variables se rapportant a l'image elle-meme */
  Raster r;
  int    w, h;	/* nombre de lignes et de colonnes de l'image */

  /* Variables liees au traitement de l'image */
  int 	 filtre;		/* numero du filtre */
  int 	 nbiter;		/* nombre d'iterations */

  /* Variables liees au chronometrage */
  double  debut, fin;  

  /* Variables de boucle */
  int 	i,j;

  /* Variable pour la parallelisation */
  unsigned char *local_data;
  int local_h;
  int p, np;
  MPI_Status status;
#define ROOT 0


  MPI_Init( &argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &p);
  MPI_Comm_size( MPI_COMM_WORLD, &np);

  if (argc != 4) {
    fprintf( stderr, usage, argv[0]);
    return 1;
  }
      
  /* Saisie des param�tres */
  filtre = atoi(argv[2]);
  nbiter = atoi(argv[3]);
        
  /* Lecture du fichier Raster */
  if( p == ROOT) {
    lire_rasterfile( argv[1], &r);
    h = r.file.ras_height;
    w = r.file.ras_width;

    if ( h % np != 0 ) {
      printf("La hauteur de l'image %d n'est pas divisible",h);
      printf("par le nombre de processus %d \n",np);
      exit(1);
    }
  }

  /* debut du chronometrage */
  debut = MPI_Wtime();            

  /* Diffusion de nbl et nbc, seulement connu du proc
   * qui charge l'image */
  MPI_Bcast( &w, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &h, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Le tampon est �ventuellement augment� de deux lignes 
   * pour le calcul de la convolution au bord : */
  local_h = h/np + (p != 0) + (p != np-1);

  /* Allocation memoire du tampon : */
  local_data = (unsigned char *)malloc( sizeof(unsigned char)*w*local_h);
  if (local_data == NULL){
    fprintf(stderr, "Erreur pour l'allocation memoire de 'local_data'.\n");
    exit(1);
  }

  /* Repartition de l'image sur tous les processus: */
  MPI_Scatter( r.data, w*h/np, MPI_UNSIGNED_CHAR,
	       local_data + w*(p!=0), w*h/np, MPI_UNSIGNED_CHAR,
	       ROOT, MPI_COMM_WORLD);
    
  /* La convolution a proprement parler */
  for(i=0 ; i < nbiter ; i++) {
    double t_com_start, t_com_end, t_comp_start, t_comp_end;

    /* Echange des lignes adjacentes */
    t_com_start = MPI_Wtime();
    if( p > 0) {
      MPI_Send( local_data + w, w, MPI_UNSIGNED_CHAR, p-1, 0, MPI_COMM_WORLD);
      MPI_Recv( local_data, w, MPI_UNSIGNED_CHAR, p-1, 0, MPI_COMM_WORLD, &status);
    }
    if( p < np-1) {
      MPI_Send( local_data + w*(local_h - 2), w, MPI_UNSIGNED_CHAR, p+1, 0, MPI_COMM_WORLD);
      MPI_Recv( local_data + w*(local_h - 1), w, MPI_UNSIGNED_CHAR, p+1, 0, 
		MPI_COMM_WORLD, &status);
    }
    t_com_end = MPI_Wtime();

    t_comp_start = MPI_Wtime();
    convolution( filtre, local_data, local_h, w);
    t_comp_end = MPI_Wtime();
    
    if (p==1){
      printf("Proc 1: i=%d t_com = %g t_comp = %g \n", i, t_com_end - t_com_start, t_comp_end - t_comp_start);
    }
    MPI_Barrier(MPI_COMM_WORLD); // sinon les resultats sont moins "clairs"...
  }


  MPI_Gather( local_data + w*(p!=0), w*h/np, MPI_UNSIGNED_CHAR,
	      r.data, w*h/np, MPI_UNSIGNED_CHAR,
	      ROOT, MPI_COMM_WORLD);

  /* Liberation memoire du tampon : */
  free(local_data);

  if( p == ROOT) {
    fin = MPI_Wtime();
    printf("Temps total de calcul : %g seconde(s) \n", fin - debut);
    /* Sauvegarde du fichier Raster */
    sauve_rasterfile("post-convolution-par.ras",&r);
  }


  MPI_Finalize();
  return 0;
}

