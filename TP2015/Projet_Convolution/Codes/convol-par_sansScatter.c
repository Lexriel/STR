/*
 * Universite Pierre et Marie Curie
 * Calcul de convolution sur une image.
 *
 * Compilation, execution : voir fichier Makefile
 * Version parallele sans Scatter
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>   /* pour le rint */
#include <string.h> /* pour le memcpy */

#include <time.h>   /* chronometrage */
#include <mpi.h>

#include "rasterfile.h"

#define MAX(a,b) ((a>b) ? a : b)
#define STOP -1


/** 
 * \struct Raster
 * Structure d�crivant une image au format Sun Raster
 */

typedef struct {
  struct rasterfile file;  ///< Entete image Sun Raster
  unsigned char rouge[256],vert[256],bleu[256];  ///< Palette de couleur
  unsigned char *data;    ///< Pointeur vers l'image
} Raster;


/**
 * Cette procedure convertit un entier LINUX en un entier SUN 
 *
 * \param i pointeur vers l'entier à convertir
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
 * avec les données liée à l'image. Le champ r.file contient
 * les informations de l'entete de l'image (dimension, codage, etc).
 * Le champ r.data est un pointeur, alloué par la fonction
 * lire_rasterfile() et qui contient l'image. Cette espace doit
 * être libéré après usage.
 *
 * \param nom nom du fichier image
 * \param r structure Raster qui contient l'image
 *  chargée en mémoire
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
 * Réalise une division d'entiers plus précise que
 * l'opérateur '/' : renvoie le quotient entier le plus proche
 * du quotient réel. 
 * Remarque : la fonction rint provient de la librairie 
 * mathématique.
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
  CONVOL_MEDIAN    ///< Filtre médian
} filtre_t;

/**
 * Réalise une opération de convolution avec un noyau prédéfini sur
 * un point.
 *
 * \param choix type de noyau pour la convolution :
 *  - CONVOL_MOYENNE1 : filtre moyenneur
 *  - CONVOL_MOYENNE2 : filtre moyenneur avec un poid central plus fort
 *  - CONVOL_CONTOUR1 : filtre extracteur de contours (laplacien)
 *  - CONVOL_CONTOUR2 : filtre extracteur de contours (max des composantes du gradient)
 *  - CONVOL_MEDIAN : filtre médian (les 9 valeurs sont triées et la valeur
 *     médiane est retournée).
 * \param NO,N,NE,O,CO,E,SO,S,SE: les valeurs des 9 points
 *  concernés pour le calcul de la convolution (cette dernière est
 *  formellement une combinaison linéaire de ces 9 valeurs).
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
 * Convolution d'une image par un filtre prédéfini
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
   * et la derniere colonne ne sont copiées (ce qui force a faire
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
  Raster 		r;
  int		nbl, nbc;	/* nombre de lignes et de colonnes de l'image */

  /* Variables liees au traitement de l'image */
  int 		filtre;		/* numero du filtre */
  int 		nbiter;		/* nombre d'iterations */

  /* Variables de la grille locale */
  unsigned char  *region;		/* grille locale */
  int		region_nb_ligne;/* nombre de lignes de la grille locale */
  int		nb_ligne_trait;	/* nombre de lignes de la grille locale a traiter */
  int		ligne_suiv_env;	/* numero de la ligne a envoyer au processus suivant */
  int		ligne_prec_env;	/* numero de la ligne a envoyer au processus precedent */
  int		ligne_suiv_recv;/* numero de la ligne a recevoir du processus suivant */
  int		ligne_prec_recv;/* numero de la ligne a recevoir du processus precedent */

  /* Variables liees au chronometrage */
  double          debut, fin;  

  /* Variables liees aux communications */

  int             rang;           /* rang du processus */
  int             p;              /* nombre de processus dans l'application */
  int             source;         /* numero du processus du maitre */
  int 		precedent, suivant;/* numero des processus de qui l'on doit recevoir */
  int             tag;            /* etiquette */
  MPI_Status      status;

  /* Variables de boucle */
  int 	i,j;


  /* Initialisation de MPI */

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);

  if (argc != 4) {
    printf("Usage : %s <nom generique fichier SunRaster> [0|1|2|3|4] <nbiter>\n",argv[0]);
    exit(1);
  }
  
  tag = 0;
  source = 0;
	
  filtre = atoi(argv[2]);
  nbiter = atoi(argv[3]);
	
  if ( rang == source ) {
    /* Lecture du fichier Raster */
    lire_rasterfile(argv[1],&r);
    nbl = r.file.ras_height;
    nbc = r.file.ras_width;
	  
    if ( nbl % p != 0 ) {
      printf("La hauteur de l'image %d n'est pas divisible",nbl);
      printf("par le nombre de processus %d \n",p);
      exit(1);
    }

  }

  /* debut du chronometrage */
  debut = MPI_Wtime();            

  /* Diffusion des variables nbl et nbc */
  MPI_Bcast(&nbl, 1, MPI_INT, source, MPI_COMM_WORLD);
  MPI_Bcast(&nbc, 1, MPI_INT, source, MPI_COMM_WORLD);
	
  region_nb_ligne = nbl / p;
  nb_ligne_trait  = region_nb_ligne + 2;
	
  ligne_suiv_env = region_nb_ligne;
  ligne_suiv_recv = region_nb_ligne + 1;
  ligne_prec_env = 1;
  ligne_prec_recv = 0;
  if ( rang ==  source ) {
    ligne_prec_env = STOP;
    ligne_prec_recv = STOP;
  } else if (rang == (p-1) ) {
    ligne_suiv_env = STOP;
    ligne_suiv_recv = STOP;
  }
	
  /* Allocation memoire de la region locale : */
  region = (unsigned char *)calloc( nbc * nb_ligne_trait,
				    sizeof(unsigned char));
  if ( region == NULL) {
    printf("Erreur dans allocation de region \n");
    exit(1);
  }
	

  /* Ma remarque : ici on envoie la grille locale a chaque processus
   * directement avec les bords des voisins necessaires a la premiere
   * iteration. 
   * C'est peut etre plus efficace (moins de comms) mais ca complique le code 
   * (cf. corps boucle "for(i=0 ; i < nbiter ; i++)") 
   * Pour une version où on envoie pas les bords avec la grille locale,  
   * voir [my_]convol-par2.c (avec MPI_Gather/MPI_Scatter). */ 

  if( rang == source) {
    /* Envoi de la grille locale a chaque processus */
    for (i = 1; i < p-1; i++) {
      MPI_Send(&r.data[(i*region_nb_ligne-1) * nbc],
	       nbc * (region_nb_ligne+2), MPI_UNSIGNED_CHAR, i, tag,
	       MPI_COMM_WORLD);
    }

    if (p != 1 ) {
      /* Traitement specifique pour le processus de rang p-1 */
      MPI_Send(&r.data[((p-1)*region_nb_ligne-1) * nbc],
	       nbc * (region_nb_ligne+1), MPI_UNSIGNED_CHAR, p-1, tag,
	       MPI_COMM_WORLD);
    }
	  
    /* Mise a jour de sa propre grille locale */
    memcpy(region+nbc,r.data, nbc * (region_nb_ligne + 1)* sizeof(unsigned char));
  } else if (rang == (p-1) ) {
    /* Reception de la grille locale */
    MPI_Recv(region, (region_nb_ligne+1) *nbc, MPI_UNSIGNED_CHAR, MPI_ANY_TAG,
	     source, MPI_COMM_WORLD, &status);
  } else 
    /* Reception de la grille locale */
    MPI_Recv(region, (region_nb_ligne+2) *nbc, MPI_UNSIGNED_CHAR, MPI_ANY_TAG,
	     source, MPI_COMM_WORLD, &status);
	
  /* Calcul des processus suivant et precedent */
  precedent = ( rang == 0 ? STOP :  rang-1);
  suivant = ( rang == (p-1) ? STOP : rang+1 );
	
  /* Convolution a proprement parler */
  for(i=0 ; i < nbiter ; i++) {
	  
    if ( rang == source ){
      convolution(filtre,&region[nbc],region_nb_ligne+1,nbc);
    }
    else {
      if ( rang == (p-1) ){
	convolution(filtre,region,region_nb_ligne+1,nbc);
      }
      else {
	convolution(filtre,region,nb_ligne_trait,nbc);
      }
    }
	  
    if ( i != ( nbiter - 1 ) ) {

      /* Envoi au processus precedent de sa premiere ligne  */		
      if (precedent != STOP )
	MPI_Send(region+ligne_prec_env *nbc, nbc, MPI_UNSIGNED_CHAR,
		 precedent, tag, MPI_COMM_WORLD);
		  		  
      /* Envoi au processus suivant de sa deniere ligne  */
      if (suivant != STOP )
	MPI_Send(region+ligne_suiv_env * nbc, nbc, MPI_UNSIGNED_CHAR,
		 suivant, tag, MPI_COMM_WORLD);
		
		
      /* Reception du processus precedent de sa ligne 0  */
      if (precedent != STOP )
	MPI_Recv(region+ligne_prec_recv * nbc, nbc, MPI_UNSIGNED_CHAR,
		 precedent, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		
		
      /* Reception du processus suivant de sa ligne nb_ligne_trait  */
      if (suivant != STOP )
	MPI_Recv(region+ligne_suiv_recv * nbc, nbc, MPI_UNSIGNED_CHAR,
		 suivant, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		
    }
	  
  } /* for i */


  /* Reception des grille locales */
  MPI_Gather(&region[nbc], region_nb_ligne*nbc, MPI_UNSIGNED_CHAR, r.data,
	     region_nb_ligne*nbc, MPI_UNSIGNED_CHAR, source, MPI_COMM_WORLD);
	
  /* Liberation memoire de la region locale : */	
  free(region);	
	
  MPI_Finalize();
	
  if( rang == source ) {
    fin = MPI_Wtime();
    printf("Temps total de calcul : %g seconde(s) \n", fin - debut);
    /* Sauvegarde du fichier Raster */
    sauve_rasterfile("post-convolution-par_sansScatter",&r);
  }
   
}

