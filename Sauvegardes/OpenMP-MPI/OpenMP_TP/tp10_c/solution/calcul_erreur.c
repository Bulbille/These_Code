#include <math.h>

static void matmul(double* restrict AA, double* restrict BB, double* restrict CC, int nn)
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

double calcul_erreur(double* restrict AA, double* restrict BB, double* restrict CC, int nn)
{
  double DD[nn][nn];
  matmul(AA, BB, &DD[0][0], nn); // D = A * B
 
  double erreur = 0.0;
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < nn; j++) {
      double e = CC[i * nn + j] - DD[i][j];
      erreur += e * e;
    }
  }
  erreur = sqrt(erreur) / nn;

  return erreur;
}
