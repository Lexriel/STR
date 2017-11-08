/*
 * Universit� Pierre et Marie Curie
 * Calcul de l'ensemble de Mandelbrot,
 * Version Parallele avec equilibrage statique de la charge
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>	/* chronometrage */
#include <string.h>     /* pour memset */
#include <math.h>
#include <sys/time.h>

#include "rasterfile.h"

#include <mpi.h>

// #define USE_MPI_GATHER



char info[] = "\
Usage:\n\
      mandel dimx dimy xmin ymin xmax ymax prof\n\
\n\
      dimx,dimy : dimensions de l'image a generer\n\
      xmin,ymin,xmax,ymax : domaine a calculer dans le plan complexe\n\
      prof : nombre maximale d'iteration\n\
\n\
Quelques exemples d'execution\n\
      mandel 800 800 0.35 0.355 0.353 0.358 200\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 500\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 300\n\
      mandel 800 800 -1.48478 0.00006 -1.48440 0.00044 100\n\
      mandel 800 800 -1.5 -0.1 -1.3 0.1 10000\n\
";



/**
 * Convertion entier (4 octets) LINUX en un entier SUN
 * @param i entier � convertir
 * @return entier converti
 */

int swap(int i) {
  int init = i; 
  int conv;
  unsigned char *o, *d;
	  
  o = ( (unsigned char *) &init) + 3; 
  d = (unsigned char *) &conv;
  
  *d++ = *o--;
  *d++ = *o--;
  *d++ = *o--;
  *d++ = *o--;
  
  return conv;
}


/*** 
 * Par Francois-Xavier MOREL (M2 SAR, oct2009): 
 */

unsigned char power_composante(int i, int p) {
  unsigned char o;
  double iD=(double) i;

  iD/=255.0;
  iD=pow(iD,p);
  iD*=255;
  o=(unsigned char) iD;
  return o;
}

unsigned char cos_composante(int i, double freq) {
  unsigned char o;
  double iD=(double) i;
  iD=cos(iD/255.0*2*M_PI*freq);
  iD+=1;
  iD*=128;
  
  o=(unsigned char) iD;
  return o;
}

/*** 
 * Choix du coloriage : definir une (et une seule) des constantes
 * ci-dessous :  
 */
//#define ORIGINAL_COLOR
#define COS_COLOR 

#ifdef ORIGINAL_COLOR
#define COMPOSANTE_ROUGE(i)    ((i)/2)
#define COMPOSANTE_VERT(i)     ((i)%190)
#define COMPOSANTE_BLEU(i)     (((i)%120) * 2)
#endif /* #ifdef ORIGINAL_COLOR */
#ifdef COS_COLOR
#define COMPOSANTE_ROUGE(i)    cos_composante(i,13.0)
#define COMPOSANTE_VERT(i)     cos_composante(i,5.0)
#define COMPOSANTE_BLEU(i)     cos_composante(i+10,7.0)
#endif /* #ifdef COS_COLOR */


/**
 *  Sauvegarde le tableau de donn�es au format rasterfile
 *  8 bits avec une palette de 256 niveaux de gris du blanc (valeur 0)
 *  vers le noir (255)
 *    @param nom Nom de l'image
 *    @param largeur largeur de l'image
 *    @param hauteur hauteur de l'image
 *    @param p pointeur vers tampon contenant l'image
 */

void sauver_rasterfile( char *nom, int largeur, int hauteur, unsigned char *p) {
  FILE *fd;
  struct rasterfile file;
  int i;
  unsigned char o;

  if ( (fd=fopen(nom, "w")) == NULL ) {
	printf("erreur dans la creation du fichier %s \n",nom);
	exit(1);
  }

  file.ras_magic  = swap(RAS_MAGIC);	
  file.ras_width  = swap(largeur);	  /* largeur en pixels de l'image */
  file.ras_height = swap(hauteur);         /* hauteur en pixels de l'image */
  file.ras_depth  = swap(8);	          /* profondeur de chaque pixel (1, 8 ou 24 )   */
  file.ras_length = swap(largeur*hauteur); /* taille de l'image en nb de bytes		*/
  file.ras_type    = swap(RT_STANDARD);	  /* type de fichier */
  file.ras_maptype = swap(RMT_EQUAL_RGB);
  file.ras_maplength = swap(256*3);

  fwrite(&file, sizeof(struct rasterfile), 1, fd); 
  
  /* Palette de couleurs : composante rouge */
  i = 256;
  while( i--) {
    o = COMPOSANTE_ROUGE(i);
    fwrite( &o, sizeof(unsigned char), 1, fd);
  }

  /* Palette de couleurs : composante verte */
  i = 256;
  while( i--) {
    o = COMPOSANTE_VERT(i);
    fwrite( &o, sizeof(unsigned char), 1, fd);
  }

  /* Palette de couleurs : composante bleu */
  i = 256;
  while( i--) {
    o = COMPOSANTE_BLEU(i);
    fwrite( &o, sizeof(unsigned char), 1, fd);
  }

  // pour verifier l'ordre des lignes dans l'image : 
  //fwrite( p, largeur*hauteur/3, sizeof(unsigned char), fd);
  
  // pour voir la couleur du '0' :
  // memset (p, 0, largeur*hauteur);
  
  fwrite( p, largeur*hauteur, sizeof(unsigned char), fd);
  fclose( fd);
}



/**
 * 	Etant donne les coordonnees d'un point a+ib dans le plan
 * 	complexe, retourne la couleur correspondante estimant
 *	a quel distance de l'ensemble de mandelbrot le point est:
 *	En definisant la suite  ainsi (Z=x+iy):
<pre>
 		Z(0) = 0
 		Z(n+1) = Z(n) * Z(n) - (a +ib)
</pre>
 *  
 *	le nombre d'iterations que la suite met pour diverger
 *	correspond a une couleur dans la palette des couleurs
 *
 *  @param a,b : coordonnees dans le plan complexe du point
 *               a calculer
 *  @param prof : nombre maximale d'iteration
 *  @return couleur du point calcule
 */

unsigned char xy2color(double a, double b, int prof) {
  double x, y, temp, x2, y2;
  int i;

  x = y = 0.;
  for( i=0; i<prof; i++) {
    /* garder la valeur pr�c�dente de x qui va etre ecrase */
    temp = x;
    /* nouvelles valeurs de x et y */
    x2 = x*x;
    y2 = y*y;
    x = x2 - y2 + a;
    y = 2*temp*y + b;
    if( x2 + y2 >= 4.0) break;
  }
  return (i==prof)?255:(int)((i%255)); 
}

/* 
 * Partie principale
 */

int main(int argc, char *argv[]) {
  /* Domaine de calcul dans le plan complexe */
  double xmin, ymin;
  double xmax, ymax;
  /* Dimension de l'image */
  int w,h;
  /* Pas d'incrementation */
  double xinc, yinc;
  /* Profondeur d'iteration */
  int prof;
  /* Image resultat */
  unsigned char	*ima, *pima;
  /* Variables intermediaires */
  int  i, j;
  double x, y;
  /* Chronometrage */
  double debut, fin;

  /* Variables liees a la parallelisation */
#define R_MASTER     0   /* Rang du processeur maitre 
			  * (attention : code ci-dessous pas correct si R_MASTER != 0 */
  int np, p;
  int local_h;
  unsigned char *local_ima;
  MPI_Status status;

  /* Initialisation MPI */
  MPI_Init( &argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &p);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if( np < 1 ) {
    /* seul le maitre parle */
    if( p == R_MASTER) fprintf( stderr, "Necessite au moins 2 processus.\n");
    MPI_Finalize();
    return 0;
  }

  if( argc == 1) {
    if( p == R_MASTER) fprintf( stderr, "%s\n", info);
  }

  /* Valeurs par defaut de la fractale */
  xmin = -2; ymin = -2;
  xmax =  2; ymax =  2;
  w = h = 800;
  prof = 10000;
  
  /* Recuperation des parametres */
  if( argc > 1) w    = atoi(argv[1]);
  if( argc > 2) h    = atoi(argv[2]);
  if( argc > 3) xmin = atof(argv[3]);
  if( argc > 4) ymin = atof(argv[4]);
  if( argc > 5) xmax = atof(argv[5]);
  if( argc > 6) ymax = atof(argv[6]);
  if( argc > 7) prof = atoi(argv[7]);

  /* Calcul des pas d'incrementation */
  xinc = (xmax - xmin) / (w-1);
  yinc = (ymax - ymin) / (h-1);
  
  /* debut du chronometrage */
  debut = MPI_Wtime();
  
  if( p == R_MASTER) {
    /* affichage parametres pour verification */
    fprintf( stderr, "Domaine: {[%lg,%lg]x[%lg,%lg]}\n", xmin, ymin, xmax, ymax);
    fprintf( stderr, "Increment : %lg %lg\n", xinc, yinc);
    fprintf( stderr, "Prof: %d\n",  prof);
    fprintf( stderr, "Dim image: %dx%d\n", w, h);
    fprintf( stderr, "#ligne/proc: %d\n", h/np);
    fprintf( stderr, "nproc: %d\n", np);
  } 

  /* Chaque proccessus calcule sa portion de taille w*nlines 
   * dans un tampon local */
  
  if (h % np != 0){
    /* seul le maitre parle */
    if( p == R_MASTER) fprintf( stderr, "Le nombre de processus ne divise pas le nombre de lignes.\n");
    MPI_Finalize();
    return 0;
  }
  local_h = h / np;
  pima = local_ima = malloc( w * local_h * sizeof( unsigned char));
  
  y = ymin + p * local_h * yinc;
  for ( i = 0; i < local_h; i++) {	
    x = xmin;
    for (j = 0; j < w; j++) {
      *pima++ = xy2color( x, y, prof); 
      x += xinc;
    }
    y += yinc;
  }

  /* Allocation memoire du tableau resultat */  
  if( p == R_MASTER) {
    ima = (unsigned char *)malloc( w*h*sizeof(unsigned char));
    if( ima == NULL) {
      fprintf( stderr, "Erreur allocation memoire du tableau \n");
      MPI_Finalize();
      return 0;
    }
  }



#ifdef USE_MPI_GATHER

  /* Attention : si le MPI_Gather implique une sychronisation de tous les processus (d�pend de l'impl�mentation 
   * et de la taille des messages), tous les processus vont terminer en meme temps : � �viter donc... */

  MPI_Gather( local_ima, w * local_h, MPI_UNSIGNED_CHAR,
	      ima, w * local_h, MPI_UNSIGNED_CHAR,
  	      R_MASTER, MPI_COMM_WORLD);

#else

  if( p == R_MASTER) {
    int source;
    MPI_Status status;

    memcpy( ima, local_ima, w * local_h);

#if 0
    /*** 1ere version : naive ***/
    /* Probleme de cette version : on empeche un processus 
     * qui a termine plus tot de nous le signaler. */ 
    for( source = 1; source < np; source ++){
      MPI_Recv( ima + w * local_h * source, w * local_h, 
		MPI_UNSIGNED_CHAR, source, 0 /* tag */, MPI_COMM_WORLD, &status);
      printf("Re�u %d lignes commen�ant � la ligne %d du processus %d\n",
	     local_h + (source == np-1 ? h%np : 0), 
	     local_h * source, 
	     source);
    } /* for source */
#endif 

    /*** 2ieme version : amelioree ***/
    /* Cette version permet de r�cup�rer en premier les r�sultats des processus 
     * qui ont termin� avant les autres (sans tenir compte de leur rang). 
     * On a besoin de MPI_Probe() pour conna�tre l'�metteur. 
     * (puis MPI_Recv() pour r�cup�rer le message correspondant). */
    for( source = 1; source < np; source ++){
      int s = 0;
      MPI_Probe(MPI_ANY_SOURCE, 0 /* tag */, MPI_COMM_WORLD, &status);
      s = status.MPI_SOURCE;
      MPI_Recv( ima + w * local_h * s, w * local_h, 
		MPI_UNSIGNED_CHAR, s, 0 /* tag */, MPI_COMM_WORLD, &status);
      printf("Re�u %d lignes commen�ant � la ligne %d du processus %d\n",
	     local_h + (s == np-1 ? h%np : 0), 
	     local_h * s, 
	     s);
    } /* for source */


  } 
  else {
    MPI_Send( local_ima, w * local_h, MPI_UNSIGNED_CHAR, R_MASTER, 0, 
	      MPI_COMM_WORLD);
  }

#endif 

  if( p == R_MASTER) {    
    /* fin du chronometrage */
    fin = MPI_Wtime();
    fprintf( stderr, "Temps total de calcul : %g sec\n", fin - debut);

    /* Sauvegarde de la grille dans le fichier resultat "mandel.ras" */
    sauver_rasterfile( "mandel.ras", w, h, ima);    
  } 
  else {
    fin = MPI_Wtime();
    fprintf( stderr, "Proc. %d: %g sec\n", p, fin - debut);
  }
  MPI_Finalize();  
  return 0;
}
