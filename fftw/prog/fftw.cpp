//Basic stuff
#include <iostream> //standard Input-Output
#include <fstream>  //standard Files
#include <sys/file.h> //Files lock
#include <iomanip>  //setprecision
#include <sstream>  //stringstream
#include <string>   //strings
//Benchamrk
#include <chrono>
//Mathematics
#include <cmath>
#include <complex>

using namespace std;
complex<double> im = complex<double>(0,1);
double FORW = 1.0;
double BACK = -1.0;

void fftw(complex<double>* in, complex<double>* out,int N,double SIG){

    for(int q = 0; q<N ; q++){
        complex<double> sum = complex<double>(0,0);
                double k    = 2*M_PI*q/N;

        for(int n = 0; n<N ; n++){
            sum += exp(SIG * im * k * static_cast<double>(n)) * in[n] ;
        }
//        cout << q << " " << sum << " " << norm(sum) << "\n";
        out[q] = (SIG > 0 ) ? sum : sum/static_cast<double>(N);
    }
}



int main(){
    int LX = 15;
    complex<double>* test = new complex<double>[LX];
    for(int i=0;i<LX;i++)
        test[i] = exp(-static_cast<double>(i)/2);

    complex<double>* res = new complex<double>[LX];
    complex<double>* tmp = new complex<double>[LX];

    fftw(test,tmp,LX,FORW);
    fftw(tmp,res,LX,BACK);

    for(int i=0;i<LX;i++)
        cout << norm(test[i]) << " " << norm(tmp[i]) << " " << norm(res[i]) << "\n";

return 0;
}
