#include <iostream> //cout
#include <math.h>
#include <gsl/gsl_fft_complex.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define ABS(z,i)  pow(pow(REAL(z,i),2)+pow(IMAG(z,i),2),0.5)
#define TOT(z,i,l)  pow(pow(z[i],2)+pow(z[l-i],2),0.5)

#include <fstream> //files
#include <complex>
#include <gsl/gsl_fft_real.h>
#define J dcomplex(0.0,1.0)
typedef std::complex<double> dcomplex;
using namespace std;

/**** Définitions des fonctions ****/
double getEquilibrium(int l,double chim);
bool EsosKaw(int* array, int x, int sens, double kbeta,double f);
void FnCorrel(double* cor, double* val, int tmax, int tcor);
double IntCor(double* cor, int tmax);
void SpCorrel(double*cor,double*val,int LX,double mbar);

/***********************************/
double gaussian(double x,double moy,double mean,int max,int L){
    double a = max*(2.*x/L-1);
    return exp(-pow(x-moy,2)/2/pow(mean,2));
    return sin(moy*x);
    return 3;
    return (x==moy) ? 1 : 0 ;
    return (abs(a-moy)<mean) ? 1 : 0;
    return exp(-mean*x*x);
}
int main(int argc,char* argv[]){
    int L = 2<<7;
    double max = 4;
    double moy = 80,mean = 20;
//    double* hauteur = new double[2*L];
//    for(int x=0;x<L;x++){
//        REAL(hauteur,x) = gaussian(x,moy,mean);
////        REAL(hauteur,x) = exp(-2*mean*abs(x));
//        IMAG(hauteur,x) = 0;
//    }
//
//    ofstream fmag("plot",std::ofstream::out);
//    gsl_fft_complex_radix2_forward (hauteur, 1, L);
//    for(int x=0;x<L;x++)
//        fmag <<  x << " "  <<  ABS(hauteur,x) << " " << REAL(hauteur,x) << " " << IMAG(hauteur,x) << endl;
//    fmag.close();

    double* haumai = new double[L];
    ofstream vmag("vrai",std::ofstream::out);
    for(int x=0;x<L;x++){
        haumai[x] = gaussian(x,moy,mean,max,L);
        vmag << x << " " << haumai[x] << endl;
    }
    vmag.close();
    ofstream rmag("fft_real",std::ofstream::out);
    gsl_fft_real_radix2_transform(haumai,1,L);
    for(int x=0;x<L/2;x++){
        rmag << x << " " << TOT(haumai,x,L)<< " " << haumai[x] << " " << haumai[L-x]  << endl;
        //rmag <<  x*max/L << " "  <<  TOT(hauteur,x,L) << endl;
    }
    rmag.close();

return 0;
}
