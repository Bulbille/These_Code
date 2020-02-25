////////////////////////////////////////////////////////////////////////////////
// dependance.c --- TP9 : boucles avec dépendances
//
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Mon Aug 26 13:53:54 2013
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

int main() {

  const int nx=VAL_NX, ny=VAL_NY;
  double V[ny][nx], W[ny][nx];

#ifdef _OPENMP
  int nb_taches;
#pragma omp parallel
  // Affichage du nombre de taches actives
  { nb_taches = omp_get_num_threads(); }
  fprintf(stdout, "\n\n   Execution dependance en parallele avec %d threads\n", nb_taches);
#endif // _OPENMP 

  // Intialisation des matrices V et W
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++) {
      V[j][i] = (i+j)/(double)(nx+ny);
      W[j][i] = (i+j)/(double)(nx+ny);
    }
  }

  // Temps CPU de calcul initial.
  clock_t t_cpu_0 = clock();

  // Temps elapsed de reference.
  struct timeval t_elapsed_0;
  gettimeofday(&t_elapsed_0, NULL);

  // Allocation et initialisation du tableau tab_sync
  // qui gere la synchronisation entre les differents threads
#ifdef _OPENMP
  int* tab_sync = malloc(nb_taches*sizeof(int));
  for(int k=0; k<nb_taches; k++) tab_sync[k] = 1;
#endif // _OPENMP

  // Boucles avec dependance
#pragma omp parallel shared(nb_taches, tab_sync)
  {
    
#ifdef _OPENMP
    int rang = omp_get_thread_num();
#endif // _OPENMP
    
    for(int j=1; j<ny; j++) {
      
#ifdef _OPENMP
      // Synchronisation des threads
      if (rang != 0) {
	while(1) {
#pragma omp flush(tab_sync) 
	  if (tab_sync[rang-1] >= tab_sync[rang]+1) break;
	}
#pragma omp flush(tab_sync, V)
      }
#endif // _OPENMP
      
#pragma omp for nowait schedule(static)
      for(int i=1; i<nx; i++) {
	V[j][i] = (V[j][i] + V[j][i-1] + V[j-1][i])/3;
      }
      
#ifdef _OPENMP
      // Mis a jour du tableau tab_sync
      tab_sync[rang]++;
#pragma omp flush(tab_sync,V)
#endif // _OPENMP
      
    } // end for

  } // omp end parallel

  // Temps elapsed final
  struct timeval t_elapsed_1;
  gettimeofday(&t_elapsed_1, NULL);
  double t_elapsed = (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;
  
  // Temps CPU de calcul final
  clock_t t_cpu_1 = clock();
  double t_cpu = (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;

#ifdef _OPENMP
  free(tab_sync);
#endif // _OPENMP

  // Verification de la justesse de la parallelisation
  // Ne pas modifier cette partie SVP
  for(int j=1; j<ny; j++) {
    for(int i=1; i<nx; i++) {
      W[j][i] = (W[j][i] + W[j][i-1] + W[j-1][i]) / 3;
    }
  }
  double norme = 0.0;
  for(int j=1; j<ny; j++) {
    for(int i=1; i<nx; i++) {
      norme += (V[j][i] - W[j][i]) * (V[j][i] - W[j][i]);
    }
  }

  // Impression du resultat.
  fprintf(stdout, "\n\n"
	  "   Valeurs de nx et ny  : %8d %8d\n"
	  "   Temps elapsed        : %10.3E sec.\n"
	  "   Temps CPU            : %10.3E sec.\n"
	  "   Norme (PB si /= 0)   : %10.3E\n",
	  nx, ny, t_elapsed, t_cpu, norme
	  );

  return EXIT_SUCCESS;

}
