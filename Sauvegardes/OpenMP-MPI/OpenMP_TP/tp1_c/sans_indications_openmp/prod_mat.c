////////////////////////////////////////////////////////////////////////////////
// prod_mat.c --- TP1 : produit de matrices : C = A * B
//
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Mon Aug 26 14:56:32 2013
// Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Dern. mod. le   : Mon Aug 26 14:56:32 2013
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// Dimension par defaut de la taille des matrices
#ifndef VAL_M
#define VAL_M 701
#endif
#ifndef VAL_N
#define VAL_N 801
#endif

int main() {
  
  int m=VAL_M, n=VAL_N;
  double a[m][n], b[n][m], c[m][m];

  struct timeval t_elapsed_0, t_elapsed_1;
  double t_elapsed;

  clock_t t_cpu_0, t_cpu_1;
  double  t_cpu;


  // Temps CPU de calcul initial.
  t_cpu_0 = clock();

  // Temps elapsed de reference.
  gettimeofday(&t_elapsed_0, NULL);
    
  // Initialisation des matrices A, B et C.
  for(int i=0; i < m; i++) {
    for(int j=0; j < n; j++) {
      a[i][j] = (i+1)+(j+1);
    }
  }
    
  for(int i=0; i < n; i++) {
    for(int j=0; j < m; j++) {
      b[i][j] = (i+1)-(j+1);
    }
  }

  for(int i=0; i < m; i++) {
    for(int j=0; j < m; j++) {
      c[i][j] = 0;
    }
  }

  // Produit de matrices
  for(int i=0; i < m; i++) {
    for(int j=0; j < m; j++) {
      for(int k=0; k < n; k++) {
	c[i][j] = c[i][j] + a[i][k] * b[k][j];
      }
    }
  }
    
  // Temps elapsed final
  gettimeofday(&t_elapsed_1, NULL);
  t_elapsed = (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;
  
  // Temps CPU de calcul final
  t_cpu_1 = clock();
  t_cpu = (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;
  
  // Impression du resultat.
  fprintf(stdout, "\n\n"
	  "   Valeurs de m et n   : %5d %5d\n"
	  "   Temps elapsed       : %10.3E sec.\n"
	  "   Temps CPU           : %10.3E sec.\n"
	  "   Resultat partiel    : %10.3E %10.3E ... %10.3E %10.3E\n\n", 
	  m, n, t_elapsed, t_cpu, c[1][1], c[2][2], c[m-3][m-3], c[m-2][m-2]
	  );
  
  return EXIT_SUCCESS;
}

