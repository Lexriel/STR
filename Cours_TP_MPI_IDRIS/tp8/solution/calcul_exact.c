#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  FILE * fichier;
  int ntx,nty;
  double * u_exact;
  int iter,iterx,itery;
  double hx,hy;
  double x,y;

  /* Lecture de ntx et nty */
  fichier = fopen("poisson.data","r");
  fscanf(fichier,"%d", &ntx);
  fscanf(fichier,"%d", &nty);
  fclose(fichier);

  /* Allocation dynamique u_exact */
  u_exact = malloc(ntx*nty*sizeof(double));

  /* Initialisation de u_exact */
  for (iter=0; iter<ntx*nty; iter++)
    u_exact[iter] = 0.;

  /* Pas */
  hx = 1./(ntx+1);
  hy = 1./(nty+1);

  /* Calcul de la solution exacte */
  for (iterx=1; iterx<ntx+1; iterx++) {
    for (itery=1; itery<nty+1; itery++) {
      x = iterx*hx;
      y = itery*hy;
      u_exact[ (iterx-1)*nty + itery-1] = x*y*(x-1)*(y-1);
    }
  }

  fichier = fopen("fort.10","w");
  for (iter=0; iter<ntx*nty; iter++)
    fprintf(fichier, "%12.5e\n", u_exact[iter]);
  fclose(fichier);
  return 0;
}

