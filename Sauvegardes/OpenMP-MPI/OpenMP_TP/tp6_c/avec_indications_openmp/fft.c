////////////////////////////////////////////////////////////////////////////////
// fft.c --- TP6 :  transformée de Fourier multiple réelle-complexe 3D
//
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Mon Aug 26 15:51:03 2013
// Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Dern. mod. le   : Mon Aug 26 15:51:03 2013
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

// Interface des routines Fortran de JMFFTW
extern void scfftm_(int* isign, int* n, int* lot, double* scale, double* x,         int* ldx, double complex* y, int* ldy, double* table, double* work, int* isys);
extern void csfftm_(int* isign, int* n, int* lot, double* scale, double complex* x, int* ldx, double* y,         int* ldy, double* table, double* work, int* isys);

// Initialisation aleatoire d'un tableau
void random_number(double* array, int size) {
  for(int i=0; i<size; i++) {
    // Generation d'un nombre dans l'intervalle [0, 1[
    double r = (double)rand() / (double)(RAND_MAX - 1);
    array[i] = r;
  }
}

int main() {

  int nx=128, ny=200, nz=5000, ldx=nx+1, ldy=nx/2+1;
  double x[ldx*ny*nz], xx[ldx*ny*nz];
  double table[100+2*nx];

#ifdef _OPENMP
  int nb_taches;
#pragma omp parallel
  { nb_taches = omp_get_num_threads(); }
  fprintf(stdout, "\n\n   Execution fft en parallele avec %d threads\n", nb_taches);
#endif // _OPENMP
 
  // Initialisation d'un tableau x.
  random_number(&x[0], ldx*ny*nz); 
  memcpy(xx, x, ldx*ny*nz*sizeof(double)); // xx = x

  // Temps CPU de calcul initial.
  clock_t t_cpu_0 = clock();

  // Temps elapsed de reference.
  struct timeval t_elapsed_0;
  gettimeofday(&t_elapsed_0, NULL);

  // Calcul de la FFT (si parallele, on affecte une tranche d'éléments par tache).
  #pragma omp ......................................................
  {
    int z_debut, z_fin, z_tranche;

#ifndef _OPENMP

    // Initialisation pour le cas d'une exécution séquentielle.
    z_tranche = nz;
    z_debut   = 0;
    z_fin     = nz-1;

#else

    // Détermination du rang d'une tache.
    int rang = ...

    // Détermination du nombre totale de taches.
    // int nb_taches = ...

    // Calcul du nombre d'éléments et des indices de début et de fin
    // d'une tranche dans la direction z.
    int reste = nz % nb_taches;
    z_tranche = nz / nb_taches;
    if (rang >= reste) {
      z_debut   = rang*z_tranche+reste;
      z_fin     = (rang+1)*z_tranche+reste-1;
    } else {
      z_tranche = z_tranche+1;
      z_debut   = rang*z_tranche;
      z_fin     = (rang+1)*z_tranche-1;
    }

#endif // _OPENMP


    // Allocation dynamique de mémoire des tableaux de travail et temporaires.
    double*         work   = (double*)        malloc(((2*nx+4)*ny*z_tranche)* sizeof(double));
    double*         temp_x = (double*)        malloc((ldx*ny*z_tranche)     * sizeof(double));
    double complex* temp_y = (double complex*)malloc((ldy*ny*z_tranche)     * sizeof(double complex));

    int code = 0;

    // Définition de constantes pour appeler les routines Fortran
    int    zero = 0, one  = 1, minus_one = -1;
    double d_one = 1.0;
    double inv_nx = 1.0/nx;
    int    lot = ny*z_tranche;
    
    // Initialisation de la table des coefficients trigonometriques par une seule tache.
#pragma omp ......................................................
    scfftm_(&zero, &nx, &lot, &d_one, x, &ldx, (double complex*)x, &ldy, table, work, &code);
    
    // Calcul de la FFT directe en parallèle.
    memcpy(temp_x, x+(z_debut*ldx*ny), (ldx*ny*z_tranche)*sizeof(double)); // temp_x(:,:,:)=x(:,:,z_debut:z_fin)
    scfftm_(&one, &nx, &lot, &d_one, temp_x, &ldx, temp_y, &ldy, table, work, &code);
  
    // Calcul de la FFT inverse en parallèle.    
    csfftm_(&minus_one, &nx, &lot, &inv_nx, temp_y, &ldy, temp_x, &ldx, table, work, &code);
    memcpy(x+(z_debut*ldx*ny), temp_x, (ldx*ny*z_tranche)*sizeof(double)); // x(:,:,z_debut:z_fin)=temp_x(:,:,:);

    // Désallocation des tableaux dynamiques.
    free(work);
    free(temp_x);
    free(temp_y);

  } // omp end parallel
  
  // Temps elapsed final
  struct timeval t_elapsed_1;
  gettimeofday(&t_elapsed_1, NULL);
  double t_elapsed = (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;
  
  // Temps CPU de calcul final
  clock_t t_cpu_1 = clock();
  double t_cpu = (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;

  // Vérification du résultat.
  // ecart = maxval(abs(x(1:nx,1:ny,1:nz) - xx(1:nx,1:ny,1:nz)))/real(nx*ny*nz,kind=8)
  double ecart = 0;
  for(int k=0; k<ldx*ny*nz; k++) {
    double curr = fabs(x[k] - xx[k]);
    if (curr > ecart)
      ecart = curr; 
  }
  ecart /= nx*ny*nz;

  // Impression du resultat
  fprintf(stdout, "\n\n"
	  "   Ecart FFT |Directe - Inverse| : %10.3E\n"
	  "   Temps elapsed              : %10.3E sec.\n"
	  "   Temps CPU                  : %10.3E sec.\n",
	  ecart, t_elapsed, t_cpu
	  );

  return EXIT_SUCCESS;

}
