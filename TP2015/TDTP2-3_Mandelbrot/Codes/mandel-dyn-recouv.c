/*
 * Université Pierre et Marie Curie
 * Calcul de l'ensemble de Mandelbrot,
 * Version Parallele avec equilibrage dynamique de la charge
 * et recouvrement communications/calcul
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




char info[] = "\
Usage:\n\
      mandel dimx dimy xmin ymin xmax ymax prof nlines rmaster\n\
\n\
      dimx,dimy : dimensions de l'image a generer\n\
      xmin,ymin,xmax,ymax : domaine a calculer dans le plan complexe\n\
      prof : nombre maximale d'iteration\n\
      nlines : nombre de ligne par bloc (grain)\n\
      rmaster : rang du maître\n\
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
 * @param i entier à convertir
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
 *  Sauvegarde le tableau de données au format rasterfile
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
    /* garder la valeur précédente de x qui va etre ecrase */
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
  /* Nombre de lignes par bloc (grain) */
  int nlines;

  /* Variables liees a la parallelisation */
#define TAG_REQ  0   /* Etiquette de nos messages */
#define TAG_DATA 1
#define TAG_END  2
  int np, p;
  int r_master;      /* rang du maitre */
  MPI_Status status;

  /* Initialisation MPI */
  MPI_Init( &argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &p);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if( np < 2 ) {
    /* comme le maitre n'est pas encore défini, seul le processus #0 parle */
    if( p == 0) fprintf( stderr, "Necessite au moins 2 processus.\n");
    MPI_Finalize();
    return 0;
  }

  /* comme le maitre n'est pas encore défini, seul le processus #0 parle */
  if( argc == 1 && p == 0) fprintf( stderr, "%s\n", info);
  
  /* Valeurs par defaut de la fractale */
  xmin = -2; ymin = -2;
  xmax =  2; ymax =  2;
  w = h = 800;
  prof = 200;
  nlines = 1;
  r_master = 0;

  /* Recuperation des parametres */
  if( argc > 1) w    = atoi(argv[1]);
  if( argc > 2) h    = atoi(argv[2]);
  if( argc > 3) xmin = atof(argv[3]);
  if( argc > 4) ymin = atof(argv[4]);
  if( argc > 5) xmax = atof(argv[5]);
  if( argc > 6) ymax = atof(argv[6]);
  if( argc > 7) prof = atoi(argv[7]);
  if( argc > 8) nlines = atoi(argv[8]);
  if( argc > 9) r_master = atoi(argv[9]);

  /* Calcul des pas d'incrementation */
  xinc = (xmax - xmin) / (w-1);
  yinc = (ymax - ymin) / (h-1);




  /*********************** Processus Maitre (serveur) *************************/ 
  /* Le processus Maitre est serveur et 
   * distribue les requetes aux clients
   */
  if( p == r_master) {
    /* Chronometrage */
    double debut, fin;
    /* Divers */
    int dest, source, bloc, nbloc, i;
    int nb_directly_sent = 0;
    unsigned char *buf, *pbuf;
    MPI_Status status;

    /* debut du chronometrage */
    debut = MPI_Wtime();

    /* affichage parametres pour verification */
    fprintf( stderr, "Domaine: {[%lg,%lg]x[%lg,%lg]}\n", xmin, ymin, xmax, ymax);
    fprintf( stderr, "Increment : %lg %lg\n", xinc, yinc);
    fprintf( stderr, "Prof: %d\n",  prof);
    fprintf( stderr, "Dim image: %dx%d\n", w, h);
    fprintf( stderr, "Grain: %d lignes\n", nlines);
    fprintf( stderr, "Nproc: %d\n", np);
    fprintf( stderr, "Root: %d\n", r_master);

    buf = (unsigned char *)malloc( w*h*sizeof(unsigned char));
    if( buf == NULL) {
      fprintf( stderr, "Erreur allocation memoire du tableau \n");

      /* Fin des hostilites */
      for( dest=0; dest<np; dest++){
	if( dest != r_master){
	  MPI_Send( &dest, 0, MPI_INT, dest, TAG_END, MPI_COMM_WORLD);      
	}
      }
      MPI_Finalize();
      return 0;
    }
    

    /***** On commence par envoyer les premieres taches */
    nbloc = h / nlines;
    /*     printf("nbloc = %d\n", nbloc);  */
    if (nbloc < np){ 
      fprintf(stderr, "Le code actuel ne gère pas les cas où : nbloc < np \n"); // cf. envoi messages de terminaison

      /* Fin des hostilites */
      for( dest=0; dest<np; dest++){
	if( dest != r_master){
	  MPI_Send( &dest, 0, MPI_INT, dest, TAG_END, MPI_COMM_WORLD);      
	}
      }
      MPI_Finalize();
      return 0;
    }
    if (h % nlines != 0){ 
      fprintf(stderr, "Le code actuel ne gère pas les cas où : 'h' n'est pas divisible par 'nlines' \n");

      /* Fin des hostilites */
      for( dest=0; dest<np; dest++){
	if( dest != r_master){
	  MPI_Send( &dest, 0, MPI_INT, dest, TAG_END, MPI_COMM_WORLD);      
	}
      }
      MPI_Finalize();
      return 0;
    }

    for ( dest = 0, bloc = 0; 
	  (dest < np ) && (bloc < nbloc);
	  dest++) {
      if( dest != r_master){
	/* 	printf("Processus %d : envoi d'une requete pour bloc=%d à dest=%d\n", p, bloc, dest);  */
	MPI_Send( &bloc, 1, MPI_INT, dest, TAG_REQ, MPI_COMM_WORLD);
	bloc++;
      }
    } /* for dest,bloc */


    /***** Puis on envoie les secondes taches sans attendre*/
    for ( dest = 0;
	  (dest < np ) && (bloc < nbloc);
	  dest++) {
      if( dest != r_master){ /* juste pour faire np-1 tours de boucle... */
	/* 	printf("Processus %d : envoi d'une requete pour bloc=%d à dest=%d (tag=%d)\n", p, bloc, dest, TAG_REQ); */
	MPI_Send( &bloc, 1, MPI_INT, dest, TAG_REQ, MPI_COMM_WORLD);
	bloc++;
      }
    } /* for dest,bloc */
    nb_directly_sent = bloc;

    
    /***** On se place en attente des résultats pour envoyer
     ***** d'autres requetes : il reste nbloc bloc a calculer */
    for ( ; bloc < nbloc; bloc++) {
      /* Reception d'un numéro de bloc */
      MPI_Recv( &i, 1, MPI_INT, MPI_ANY_SOURCE, TAG_REQ, MPI_COMM_WORLD, &status);
      source = status.MPI_SOURCE;
      
      /* Reception du bloc de donnée */
      MPI_Recv( buf + i * w * nlines, w * nlines, MPI_UNSIGNED_CHAR, source,
		TAG_DATA, MPI_COMM_WORLD, &status);

      /* Envoi d'une nouvelle tache */
      /*       printf("Processus %d : envoi d'une requete pour bloc=%d à dest=%d (tag=%d)\n", p, bloc, source, TAG_REQ); */
      MPI_Send( &bloc, 1, MPI_INT, source, TAG_REQ, MPI_COMM_WORLD); /* a priori, la réception correspondante est 
									déjà lancée (cf. code "esclave") donc pas 
									besoin d'émission non bloquante... */
    }
    
 
    /***** Reception des derniers blocs : on attend le meme nombre de
     ***** message que ceux envoyés dans la premiere boucle (cf. 'nb_directly_sent') */
    for ( dest = 0; dest < nb_directly_sent ; dest++ ) {
      /* Reception d'un numéro de bloc */
      MPI_Recv( &i, 1, MPI_INT, MPI_ANY_SOURCE, TAG_REQ, MPI_COMM_WORLD, &status);
      source = status.MPI_SOURCE;

      /* Reception du bloc de donnée */
      MPI_Recv( buf + i * w * nlines, w * nlines, MPI_UNSIGNED_CHAR, source, 
		TAG_DATA, MPI_COMM_WORLD, &status); 

      /* Envoi de la terminaison */
      /*       printf("Processus %d : envoi d'une terminaison à dest=%d (tag=%d)\n", p, source, TAG_END);  */
      /*       MPI_Send( &dest /\* non utilisé en réception *\/, 1, MPI_INT,  */
      /* 		source, TAG_END, MPI_COMM_WORLD); */
      MPI_Send( NULL, 0 /* message de taille nulle */, MPI_INT,
      		source, TAG_END, MPI_COMM_WORLD);
    }

    /* fin du chronometrage */    
    fin = MPI_Wtime();

    fprintf( stderr, "Temps total de calcul : %g sec\n", fin - debut);  
    fprintf( stdout, "%g\n", fin - debut);    

    /* Sauvegarde de la grille dans le fichier resultat "mandel.ras" */
    /*     printf("Début sauvegarde.\n"); */
    sauver_rasterfile( "mandel.ras", w, h, buf);
    /*     printf("Fin sauvegarde.\n"); */
  } 
  else {
    /*********************** Processus Esclave (clients) **********************/ 

    /* Variables intermediaires */
    int  i, j;
    double x, y;
    int local_h;
    unsigned char *local_buf, *pbuf;

    MPI_Request requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
    MPI_Status statuses[2] = {0};
    int blocs[2] = {0, 0};
    int n;

    local_h = nlines;
    local_buf = malloc( w * local_h * sizeof( unsigned char));

    /*** Technique dite du "double buffering" : ***/
    MPI_Irecv(blocs, 1, MPI_INT, r_master, MPI_ANY_TAG, MPI_COMM_WORLD, requests);
    MPI_Irecv(blocs+1, 1, MPI_INT, r_master, MPI_ANY_TAG, MPI_COMM_WORLD, requests+1);

    /* Les clients sont en attente de requetes */
    n = 0;
    for(;;) {
      int bloc = 0;

      /* En attente d'une requete */
      /*       printf("Processus %d : attente d'une requete de r_master=%d\n", p, r_master);  */
      MPI_Wait(requests+n, statuses+n);	  
      /*       printf("Processus %d : réception d'un message d'étiquette=%d\n", p, statuses[n].MPI_TAG); */
      /*       printf("Processus %d : réception d'un message d'étiquette=%s\n", p, */
      /* 	     (statuses[n].MPI_TAG == TAG_REQ ? "TAG_REQ" : */
      /* 	      (statuses[n].MPI_TAG == TAG_DATA ? "TAG_DATA" : */
      /* 	       (statuses[n].MPI_TAG == TAG_END ? "TAG_END" : "inconnue!")))); */


      /* Cas d'une fin de tache */
      if( statuses[n].MPI_TAG == TAG_END){
	/* on annule l'autre réception non bloquante en cours : */
	MPI_Cancel(requests+(1-n));
	break;
      }
      else {
	/* on sauve le numero de bloc : */
	bloc = blocs[n];
	/* on relance une autre réception : */
	MPI_Irecv(blocs+n, 1, MPI_INT, r_master, MPI_ANY_TAG, MPI_COMM_WORLD, requests+n);
      }

      /* Sinon on recoit une requete de calcul d'un bloc : */
      /*       printf("Processus %d : bloc = %d\n", p, bloc); */

      /* le point de depart depend du bloc a calculer  */
      y = ymin + yinc*bloc*local_h;
      pbuf = local_buf;

      /* Chaque proccessus calcule sa portion de taille w*nlines 
       * dans un tampon local */          
      for ( i = 0; i < local_h; i++) {	
	x = xmin;
	for (j = 0; j < w; j++) {
	  *pbuf++ = xy2color( x, y, prof); 
	  x += xinc;
	}
	y += yinc;
      }
      

      /* Envoi du bloc */          
      /*       printf("Processus %d : envoi bloc=%d et données à dest=%d\n", p, bloc, r_master);   */
      MPI_Send( &bloc, 1, MPI_INT, r_master, TAG_REQ, MPI_COMM_WORLD);
      MPI_Send( local_buf, w*local_h, MPI_UNSIGNED_CHAR, r_master, TAG_DATA, MPI_COMM_WORLD);
      /*       printf("Processus %d : FIN DE envoi bloc=%d et données à dest=%d\n", p, bloc, r_master);   */

      /* mise à jour de 'n' pour la prochaine itération : */
      n = 1-n;
    } /* for */ 

    free( local_buf);
  }

  MPI_Finalize();  
  return 0;
}
