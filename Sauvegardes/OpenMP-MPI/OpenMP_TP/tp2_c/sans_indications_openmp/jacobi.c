////////////////////////////////////////////////////////////////////////////////
// jacobi.c --- TP2 : resolution d'un systeme lineaire par la methode de jacobi
//
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Mon Aug 26 14:56:32 2013
// Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Dern. mod. le   : Mon Aug 26 14:56:32 2013
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>


// Dimension par defaut de la taille des matrices
#ifndef VAL_N
#define VAL_N 1201
#endif
#ifndef VAL_D
#define VAL_D 800
#endif

// Initialisation aleatoire d'un tableau
void random_number(double* array, int size) {
  for(int i=0; i<size; i++) {
    // Generation d'un nombre dans l'intervalle [0, 1[
    double r = (double)rand() / (double)(RAND_MAX - 1);
    array[i] = r;
  }
}

int main() {
  
  int n=VAL_N, diag=VAL_D;
  int i, j, iteration=0;
  double a[n][n];
  double x[n], x_courant[n], b[n];
  double norme;

  struct timeval t_elapsed_0, t_elapsed_1;
  double t_elapsed;

  clock_t t_cpu_0, t_cpu_1;
  double  t_cpu;


  // Initialisation de la matrice et du second membre
  srand(421); // pour des resultats reproductibles
  random_number(&a[0][0], n*n);
  random_number(&b[0], n);

  // On muscle la diagonale principale de la matrice
  for(int i=0; i<n; i++) {
    a[i][i] += diag;
  }

  // Solution initiale
  for(int i=0; i<n; i++) {
    x[i] = 1.0;
  }

  // Temps CPU de calcul initial.
  t_cpu_0 = clock();

  // Temps elapsed de reference.
  gettimeofday(&t_elapsed_0, NULL);

  // Resolution par la methode de Jacobi
  while(1) {
    iteration++;
    
    for(int i=0; i<n; i++) {
      x_courant[i] = 0;
      for(j=0; j<i; j++) {
	x_courant[i] += a[j][i]*x[j];
      }
      for(j=i+1; j<n; j++) {
	x_courant[i] += a[j][i]*x[j];
      }
      x_courant[i] = (b[i] - x_courant[i])/a[i][i];
    }
    
    // Test de convergence
    {
      double absmax = 0;
      for(int i=0; i<n; i++) {
	double curr = fabs(x[i] - x_courant[i]);
	if (curr > absmax)
	  absmax = curr; 
      }
      norme = absmax / n;
    }
    if( (norme <= DBL_EPSILON) || (iteration >= n) ) break;
    
    // Copie de x_courant dans x
    memcpy(x, x_courant, n*sizeof(double));
  }
  
  // Temps elapsed final
  gettimeofday(&t_elapsed_1, NULL);
  t_elapsed = (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;
  
  // Temps CPU de calcul final
  t_cpu_1 = clock();
  t_cpu = (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;
  
  // Impression du resultat
  fprintf(stdout, "\n\n"
	  "   Taille du systeme   : %5d\n"
	  "   Iterations          : %4d\n"
	  "   Norme               : %10.3E\n"
	  "   Temps elapsed       : %10.3E sec.\n"
	  "   Temps CPU           : %10.3E sec.\n",
	  n, iteration, norme, t_elapsed, t_cpu
	  );
  
  return EXIT_SUCCESS;
}
