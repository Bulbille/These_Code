////////////////////////////////////////////////////////////////////////////////
// principal.c --- TP7 : résolution d'un système linéaire général multiple
//             --- avec la méthode du Bi_CGSTAB.
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Fri Aug 02 10:57:12 2013
// Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Dern. mod. le   : Mon Aug 26 14:57:55 2013
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

// Interface de la routine Fortran
void bi_cgstab_(int* n, int* m, double* a, int* lda, double* pr, double* b, int* ldb, double* x, int* ldx, double* norme, int* nb_iter_min, int* nb_iter_max);

// Initialisation aleatoire d'un tableau
void random_number(double* array, int size) {
  for(int i=0; i<size; i++) {
    // Generation d'un nombre dans l'intervalle [0, 1[
    double r = (double)rand() / (double)(RAND_MAX - 1);
    array[i] = r;
  }
}

int main() {

  int n=VAL_N, m=VAL_M, lda=VAL_LDA, ldb=lda, ldx=lda;
  double a[n][lda]; // column-major order (Fortran)
  double b[m][ldb];
  double x[m][ldx];
  double pr[n];
  double norme;
  int nb_iter_min, nb_iter_max;

#ifdef _OPENMP
  int nb_taches;
#pragma omp parallel
  { nb_taches = omp_get_num_threads(); }
  fprintf(stdout, "\n\n   Execution bi_cgstab en parallele avec %d threads\n", nb_taches);
#endif // _OPENMP

  // On initialise la matrice et le second membre
  random_number(&a[0][0],n*lda);
  random_number(&b[0][0],m*ldb);

  // On muscle la diagonale principale
  for(int i=0; i<n; i++) a[i][i] += 100;

  // On choisi le preconditionnement Jacobi=Diagonale(a)
  for(int i=0;i<n; i++) pr[i] = a[i][i];

  // Solution initiale.
  for(int i=0; i<m; i++) 
    for(int j=0; j<ldx; j++) 
      x[i][j] = 1;

  // Temps CPU de calcul initial.
  clock_t t_cpu_0 = clock();

  // Temps elapsed de reference.
  struct timeval t_elapsed_0;
  gettimeofday(&t_elapsed_0, NULL);

  // Resolution du systeme lineaire multiple Ax=b en parallèle.
  #pragma omp parallel
  {
    bi_cgstab_(&n, &m, &a[0][0], &lda, &pr[0], &b[0][0], &ldb, &x[0][0], &ldx, &norme, &nb_iter_min, &nb_iter_max);
  }

  // Temps elapsed final
  struct timeval t_elapsed_1;
  gettimeofday(&t_elapsed_1, NULL);
  double t_elapsed = (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;
  
  // Temps CPU de calcul final
  clock_t t_cpu_1 = clock();
  double t_cpu = (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;

  // Impression du resultat.
  fprintf(stdout, "\n\n"
	  "   Nb iterations min. et max. : %4d %4d\n"
	  "   Norme                      : %10.3E\n"
	  "   Temps elapsed              : %10.3E sec.\n"
	  "   Temps CPU                  : %10.3E sec.\n",
	  nb_iter_min, nb_iter_max, norme, t_elapsed, t_cpu
	  );

  return EXIT_SUCCESS;

}
