#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

int main()
{
  #pragma omp parallel
  {
    int rang = 0;

#ifdef _OPENMP
    #pragma omp single
    {
      int nb_taches = omp_get_num_threads();
      fprintf(stdout, "\n\n   Execution de hello_world en parallele avec %d threads\n", nb_taches);
    } // omp end single

    rang = omp_get_thread_num();
#endif // _OPENMP

    fprintf(stdout, "   Bonjour depuis le thread de rang %d\n", rang);
  } // omp end parallel

  return EXIT_SUCCESS;
}
