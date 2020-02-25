////////////////////////////////////////////////////////////////////////////////
// bi_cgstab.c --- --- TP7 : méthode du bi-gradient conjugué stabilisée
//             --- --- (Réf. Templates for the Solution of Linear Systems. p.27)
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Mon Aug 26 14:56:32 2013
// Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Dern. mod. le   : Mon Aug 26 14:56:32 2013
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <mkl_cblas.h>

#define max(a,b)	       \
  ({ __typeof__ (a) _a = (a);	\
    __typeof__ (b) _b = (b);	\
    _a > _b ? _a : _b; })

#define min(a,b)	       \
  ({ __typeof__ (a) _a = (a);	\
    __typeof__ (b) _b = (b);	\
    _a < _b ? _a : _b; })

void bi_cgstab(/* in  */ int n, int m, double* a, int lda, double* pr, double* b, int ldb, /* inout */ double* x, int ldx, 
	       /* out */ double* norme, int* nb_iter_min, int* nb_iter_max) {

  double r[n], r_tilde[n], p[n], p_hat[n], s[n], s_hat[n], v[n], t[n];
  double norme_locale, rho, rho_old, alpha, beta, omega;

  *nb_iter_min=n; *nb_iter_max=0; *norme=0;
  
  // Resolution de m systemes lineaires independants.
  for(int j=0; j<m; j++) {
    
    // Initialisation du compteur sur le nombre d'iterations.
    int iteration = 0;
    double* b_j = b + j*ldb; // b(1:n, j)
    double* x_j = x + j*ldx; // x(1:n, j)
    
    // Ecart initiale.
    // r(:) = b(1:n,j) - matmul(a(1:n,1:n),x(1:n,j))
    memcpy(&r[0], b_j, n* sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1.0, a, lda, x_j, 1, 1.0, r, 1);

    memcpy(&r_tilde[0], &r[0], n* sizeof(double)); // r_tilde(:) = r(:);

    // Resolution par la methode du bi-gradient conjugue stabilisee.
    while(1) {
      iteration++;

      rho = cblas_ddot(n, r_tilde, 1, r, 1);

      if (rho == 0) {
	fprintf(stdout, " La procedure Bi-CGSTAB a echoue, rho=0\n");
	break;
      }
      if (iteration == 1) {
	memcpy(&p[0], &r[0], n* sizeof(double)); // p(:) = r(:);
      } else {
	beta = (rho / rho_old) * (alpha / omega);

	for(int k=0; k<n; k++) {
	  p[k] = r[k] + beta * (p[k] - omega * v[k]);
	}

      }

      for(int k=0; k<n; k++) {
        p_hat[k] = p[k] / pr[k];
      }

      // v(:) = matmul(a(1:n,1:n), p_hat(:))
      cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, a, lda, p_hat, 1, 0.0, v, 1);

      alpha = rho / cblas_ddot(n, r_tilde, 1, v, 1);
      for(int k=0; k<n; k++) {
	s[k] = r[k] - alpha * v[k];
      }

      // Test de convergence 1
      norme_locale = fabs(s[cblas_idamax(n, s, 1)])/n; // maxval( abs( s(:) ) )/real(n,kind=8)
      
      if((norme_locale <= 10*DBL_EPSILON) || (iteration >= n)) {
	for(int k=0; k<n; k++) {
	  x_j[k] += alpha * p_hat[k];
	}
	break;
      }

      for(int k=0; k<n; k++) {
	s_hat[k] = s[k] / pr[k];
      }

      // t(:) = matmul(a(1:n,1:n), s_hat(:))
      cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, a, lda, s_hat, 1, 0.0, t, 1);

      omega = cblas_ddot(n, t, 1, s, 1) / cblas_ddot(n, t, 1, t, 1);
      if (omega == 0) {
	fprintf(stdout, " La procedure Bi-CGSTAB a echoue, omega=0\n");
	break;
      }
      for(int k=0; k<n; k++) {
	x_j[k] += alpha * p_hat[k] + omega * s_hat[k];
      }
        
      for(int k=0; k<n; k++) {
	r[k] = s[k] - omega * t[k];
      }

      // Test de convergence 2
      norme_locale = fabs(s[cblas_idamax(n, s, 1)])/n; // maxval( abs( s(:) ) )/real(n,kind=8)
      if((norme_locale <= 10*DBL_EPSILON) || (iteration >= n)) break;

      rho_old = rho;
    } // while(1)

    {
      *nb_iter_min = min(*nb_iter_min, iteration);
      *nb_iter_max = max(*nb_iter_max, iteration);
      *norme       = max(*norme, norme_locale);
    }

  } // for
  
}
