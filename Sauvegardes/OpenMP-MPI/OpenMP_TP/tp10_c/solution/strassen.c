#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

void strassen_multiply(int nn, double AA[][nn], double BB[][nn], double CC[][nn]);

void random_number(double* array, int size);
void matmul(double* restrict AA, double* restrict BB, double* restrict CC, int nn);
void matcopy(double* AA, int AAstride, double* CC, int CCstride, int nn);
void matadd(double* AA, int AAstride, double* BB, int BBstride, double* CC, int CCstride, int nn);
void matsub(double* AA, int AAstride, double* BB, int BBstride, double* CC, int CCstride, int nn);
double calcul_erreur(double* restrict AA, double* restrict BB, double* restrict CC, int nn);

#define THRESHOLD 128

int main()
{
  const int n = 2048;
  double A[n][n], B[n][n], C[n][n];

  struct timeval t_elapsed_0, t_elapsed_1;
  double t_elapsed;

  clock_t t_cpu_0, t_cpu_1;
  double  t_cpu;

#ifdef _OPENMP
  int nb_taches;
  #pragma omp parallel
  { nb_taches = omp_get_num_threads(); }
  fprintf(stdout, "\n\n   Execution strassen en parallele avec %d threads\n", nb_taches);
#endif // _OPENMP

  // Initialisation des matrices A et B
  srand(421); // pour des résultats reproductibles
  random_number(&A[0][0], n * n);
  random_number(&B[0][0], n * n);

  // Temps CPU de calcul initial.
  t_cpu_0 = clock();

  // Temps elapsed de reference.
  gettimeofday(&t_elapsed_0, NULL);

  #pragma omp parallel
  {
    #pragma omp single
    {
      strassen_multiply(n, A, B, C);
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
    "   Erreur        : %10.3E\n"
    "   Temps elapsed : %10.3E sec.\n"
    "   Temps CPU     : %10.3E sec.\n",
    calcul_erreur(&A[0][0], &B[0][0], &C[0][0], n), t_elapsed, t_cpu
  );

  return EXIT_SUCCESS;
}

void strassen_multiply(int nn, double AA[][nn], double BB[][nn], double CC[][nn])
{
  if ((nn & 1) || nn < THRESHOLD) {
    matmul(&AA[0][0], &BB[0][0], &CC[0][0], nn);
  } else {
    int nn2 = nn / 2;

    double Q1[nn2][nn2], Q2[nn2][nn2], Q3[nn2][nn2], Q4[nn2][nn2], Q5[nn2][nn2], Q6[nn2][nn2], Q7[nn2][nn2];

    #pragma omp task shared(AA,BB,Q1,nn,nn2) // final(nn < 128) mergeable
    {
      double A11pA22[nn2][nn2], B11pB22[nn2][nn2];
      matadd(&AA[0][0], nn, &AA[nn2][nn2], nn, &A11pA22[0][0], nn2, nn2);
      matadd(&BB[0][0], nn, &BB[nn2][nn2], nn, &B11pB22[0][0], nn2, nn2);

      strassen_multiply(nn2, A11pA22, B11pB22, Q1);
    }
    #pragma omp task shared(AA,BB,Q2,nn,nn2) // final(nn < 128) mergeable
    {
      double A21pA22[nn2][nn2], B11[nn2][nn2];
      matadd(&AA[nn2][0], nn, &AA[nn2][nn2], nn, &A21pA22[0][0], nn2, nn2);
      matcopy(&BB[0][0], nn, &B11[0][0], nn2, nn2);

      strassen_multiply(nn2, A21pA22, B11, Q2);
    }
    #pragma omp task shared(AA,BB,Q3,nn,nn2) // final(nn < 128) mergeable
    {
      double A11[nn2][nn2], B12mB22[nn2][nn2];
      matcopy(&AA[0][0], nn, &A11[0][0], nn2, nn2);
      matsub(&BB[0][nn2], nn, &BB[nn2][nn2], nn, &B12mB22[0][0], nn2, nn2);

      strassen_multiply(nn2, A11, B12mB22, Q3);
    }
    #pragma omp task shared(AA,BB,Q4,nn,nn2) // final(nn < 128) mergeable
    {
      double A22[nn2][nn2], mB11pB21[nn2][nn2];
      matcopy(&AA[nn2][nn2], nn, &A22[0][0], nn2, nn2);
      matsub(&BB[nn2][0], nn, &BB[0][0], nn, &mB11pB21[0][0], nn2, nn2);

      strassen_multiply(nn2, A22, mB11pB21, Q4);
    }
    #pragma omp task shared(AA,BB,Q5,nn,nn2) // final(nn < 128) mergeable
    {
      double A11pA12[nn2][nn2], B22[nn2][nn2];
      matadd(&AA[0][0], nn, &AA[0][nn2], nn, &A11pA12[0][0], nn2, nn2);
      matcopy(&BB[nn2][nn2], nn, &B22[0][0], nn2, nn2);

      strassen_multiply(nn2, A11pA12, B22, Q5);
    }
    #pragma omp task shared(AA,BB,Q6,nn,nn2) // final(nn < 128) mergeable
    {
      double mA11pA21[nn2][nn2], B11pB12[nn2][nn2];
      matsub(&AA[nn2][0], nn, &AA[0][0], nn, &mA11pA21[0][0], nn2, nn2);
      matadd(&BB[0][0], nn, &BB[0][nn2], nn, &B11pB12[0][0], nn2, nn2);

      strassen_multiply(nn2, mA11pA21, B11pB12, Q6);
    }
    #pragma omp task shared(AA,BB,Q7,nn,nn2) // final(nn < 128) mergeable
    {
      double A12mA22[nn2][nn2], B21pB22[nn2][nn2];
      matsub(&AA[0][nn2], nn, &AA[nn2][nn2], nn, &A12mA22[0][0], nn2, nn2);
      matadd(&BB[nn2][0], nn, &BB[nn2][nn2], nn, &B21pB22[0][0], nn2, nn2);

      strassen_multiply(nn2, A12mA22, B21pB22, Q7);
    }
    #pragma omp taskwait

    // C11 = Q1+Q4-Q5+Q7
    matadd(&Q1[0][0], nn2, &Q4[0][0], nn2, &CC[0][0], nn, nn2); // C11 = Q1+Q4
    matsub(&CC[0][0], nn, &Q5[0][0], nn2, &CC[0][0], nn, nn2); // C11 = C11-Q5
    matadd(&CC[0][0], nn, &Q7[0][0], nn2, &CC[0][0], nn, nn2); // C11 = C11+Q7
    // C21 = Q2+Q4
    matadd(&Q2[0][0], nn2, &Q4[0][0], nn2, &CC[nn2][0], nn, nn2);
    // C12 = Q3+Q5
    matadd(&Q3[0][0], nn2, &Q5[0][0], nn2, &CC[0][nn2], nn, nn2);
    // C22 = Q1+Q3-Q2+Q6
    matadd(&Q1[0][0], nn2, &Q3[0][0], nn2, &CC[nn2][nn2], nn, nn2); // C22 = Q1+Q3
    matsub(&CC[nn2][nn2], nn, &Q2[0][0], nn2, &CC[nn2][nn2], nn, nn2); // C22 = C22-Q2
    matadd(&CC[nn2][nn2], nn, &Q6[0][0], nn2, &CC[nn2][nn2], nn, nn2); // C22 = C22+Q6
  }
}

// Initialisation aleatoire d'un tableau
void random_number(double* array, int size)
{
  for (int i = 0; i < size; i++) {
    // Generation d'un nombre dans l'intervalle [0, 1[
    array[i] = (double)rand() / (double)(RAND_MAX - 1);
  }
}

void matmul(double* restrict AA, double* restrict BB, double* restrict CC, int nn)
{
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < nn; j++) {
      CC[i * nn + j] = 0.0;

      for (int k = 0; k < nn; k++) {
        CC[i * nn + j] += AA[i * nn + k] * BB[k * nn + j];
      }
    }
  }
}

void matcopy(double* AA, int AAstride, double* CC, int CCstride, int nn)
{
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < nn; j++) {
      CC[i * CCstride + j] = AA[i * AAstride + j];
    }
  }
}

void matadd(double* AA, int AAstride, double* BB, int BBstride, double* CC, int CCstride, int nn)
{
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < nn; j++) {
      CC[i * CCstride + j] = AA[i * AAstride + j] + BB[i * BBstride + j];
    }
  }
}

void matsub(double* AA, int AAstride, double* BB, int BBstride, double* CC, int CCstride, int nn)
{
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < nn; j++) {
      CC[i * CCstride + j] = AA[i * AAstride + j] - BB[i * BBstride + j];
    }
  }
}
