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
#include <random>
#include <gsl/gsl_histogram.h>


using namespace std;
//Propriétés de l'Hamiltonien
const double T_C    = 2./log(1.+sqrt(2)), J = 1;
double       Beta,ttc,H,F;
string prefix = ".",suffix="";
string algo;
//Propriétés de la grille
int       LX    = 100, LY = 500, N = LX*LY;
const int T_EQ  = 1e3;
int       T_PHO = LX*1.0e2;
const int T_DONNEE  = T_PHO*1e3;
//const long int T_MAX = T_PHO*1.2e8;
//Générateur nombres aléatoires
default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

/***********************************/
/**** Définitions des fonctions ****/
int modulo(int a, int b); 
double modulo_d(double a, int b); 
void parametres(int argc, char* argv[]);
void ecriture(int* grille,int* histogramme,int* flow,int time,double* mag,double energie, double* fourier);
bool isOnlyDouble(const char* str);
void fftw(int* in, double* out,int N,double SIG); complex<double> im = complex<double>(0,1);
void fftw(double* in, double* out,int N,double SIG); 
void ft_correl(int* in,double* out,int N,double SIG);

bool Esos(int* array, int x, int sens,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+sens,LX)],
        hn  = array[modulo(x-sens,LX)],
        hnp = array[modulo(x+2*sens,LX)];

    double sos    = abs(hn-(hx-1)) +  abs((hx-1)-(hxp+1)) +   abs((hxp+1)-hnp) - 
                   (abs(hn-hx)     +  abs(hx-hxp)         +   abs(hxp-hnp)      );
    double dg    = 3*(1+hxp-hx)+hn-hnp;
    double champ = Champ*( abs(hx-1) - abs(hx) + abs(hxp+1) - abs(hxp) );
    double shear = - J*F*sens*(hx+hxp)/2;
//    double D_e   = (hxp < LY and hx > -LY) * exp(-kbeta*(dg+shear+champ));
    double D_e   = (hxp < LY and hx >= -LY) * exp(-kbeta*(sos+shear+champ));
    return (rand_01(generator) <  D_e ) ? true : false;
}
/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);
    std::uniform_int_distribution<int> rand_lx(0,LX-1);
    std::uniform_int_distribution<int> rand_int01(0,1);
    generator.seed(time(NULL));
    int* histo  = new int[(2*LY+1)];
    int* flux   = new int[(2*LY+1)];
    double* mag = new double[(2*LY+1)]; 
    double* fourier = new double[LX];
    int* system = new int[LX]; 
    double energie  = 0;

    int photo_suiv=0; int nb_photos=0; int donnee_suiv = 0;

    for(int x = 0;x<LX;x++)     system[x] = fourier[x] = 0;
    for(int y = 0;y<=2*LY;y++)  flux[y]   = mag[y]     = histo[y] = 0;

    for(long int t = 0; true ; t++){
        int tirage  = rand_lx(generator);
        int sens    = -1+2*rand_int01(generator);

        if(Esos(system,tirage,sens,Beta,H)){
//            flux[LY+system[tirage]] += sens;
            system[tirage] -= 1;
            system[modulo(tirage+sens,LX)] += 1;
        }

        if(photo_suiv <= t && (t-1)%T_PHO == 0) {
            photo_suiv += T_PHO; nb_photos++;
            
            double* mag_t  = new double[2*LY+1]; 
            double e_par   = 0;
            for(int y = -LY;y<=LY;y++){
                mag_t[LY+y] = 0;
                for(int x = 0; x<LX;x++)
                    if(system[x] > y) mag_t[LY+y]-- ;
                    else if(system[x] < y) mag_t[LY+y]++;
            }
            for(int x = 0; x<LX;x++){
                e_par += J * abs(system[x]-system[modulo(x+1,LX)]) + H*abs(system[x]);
            }
            energie  = energie*(nb_photos-1)/nb_photos + e_par/nb_photos/static_cast<double>(LX);
            for(int y = -LY;y<=LY;y++){
                mag[LY+y]  = mag[LY+y]*(nb_photos-1)/nb_photos + mag_t[LY+y]/nb_photos/static_cast<double>(LX);
            }

            double* fourier_par = new double[LX];
  //          ft_correl(system,fourier_tmp,LX,1);
            ft_correl(system,fourier_par,LX,1);

            for(int x=0;x<LX;x++){
                histo[LY+system[x]]++;
                fourier[x] = fourier[x]*(nb_photos-1)/nb_photos + fourier_par[x]/nb_photos;
            }

            if(donnee_suiv <=t){
                donnee_suiv+=T_DONNEE;
                ecriture(system,histo,flux,t,mag,energie,fourier);
            }
            delete[] mag_t;
        }
    }
//    auto debut = std::chrono::high_resolution_clock::now();
//    auto fin = std::chrono::high_resolution_clock::now();
//    auto duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
//    std::cout << "Duree : " << duree/1e9 <<  std::endl;
return 0;
}
/***************************************************************/
/************** Fin du main ************************************/
/***************************************************************/

/*** Modulo **/
int modulo(int a, int b){ 
    while(a < 0 or a >= b){ 
        if(a<0) a+= b;
        else if(a >= b) a-=b;
    }   
    return a;
}

double modulo_d(double a, int b){ 
    while(a < 0 or a >= b){ 
        if(a<0) a+= b;
        else if(a >= b) a-=b;
    }   
    return a;
}

/*** Parametres ***/
void parametres(int argc, char* argv[]){
    Beta=1./T_C;
    ttc = 1;
    F   = 0;
    H   = 0;
    bool taille = false;
    if(argc > 1){
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "-t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    Beta = 1/(ttc*T_C);
                }
                else{ cout << "Le paramètre temperature n'est pas un double"; abort();}
            }
            else if(arg == "-h"){
                if(isOnlyDouble(argv[i+1])){
                    H = atof(argv[i+1]);
                }   
                else{ cout << "Le paramètre h n'est pas un double"; abort();}
            }   
            else if(arg == "-prefix"){
                prefix = argv[i+1];
            }   
            else if(arg == "-f"){
                if(isOnlyDouble(argv[i+1])){
                    F = atof(argv[i+1]);
                }   
                else{ cout << "Le paramètre F n'est pas un double"; abort();}
            }   
            else if(arg == "-lx"){
                if(isOnlyDouble(argv[i+1])){
                    LX = round(atof(argv[i+1]));
                    taille = true;
                }   
                else{ cout << "Le paramètre F n'est pas un double"; abort();}
            }   
            else if(arg == "-suffix"){
                suffix = argv[i+1];
            }   
        }   
    }   
    T_PHO *= LX;
    if(system(("mkdir -p " + prefix).c_str())) {;}
    algo = prefix + "/";
    {
        stringstream stream;
        stream << fixed << setprecision(2) << ttc ;
        algo += "chip_T-" + stream.str();
    }
    if(fabs(H) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(4) << H;
        algo += "-h-" + stream.str();
    }
    if(fabs(F) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(2) << F;
        algo += "-f-" + stream.str();
    }
    if(taille){
        stringstream stream;
        stream << fixed << setprecision(3) << LX;
        algo += "-lx-" + stream.str();
    }
}


/************* ÉCRITURE DANS FICHIER **********/
void ecriture(int* grille,int* histogramme,int* flow,int time, double* mag, double energie, double* fourier){
    /*********** CREATION DOSSIER POUR RESULTATS ********/
    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    string str = algo;

    ofstream result(str.c_str());
    for(int x=0;x<LX;x++) result << x << " " << grille[x] << "\n" << x+0.999 << " " << grille[x] << "\n" ;
    result.close();

    str = algo+"-histo";
    ofstream fhisto(str.c_str());
    for(int y=-LY;y<=LY;y++) fhisto << y << " " << histogramme[LY+y] << "\n" ;
    fhisto.close();

    str = algo+"-flux";
    ofstream fflux(str.c_str());
    for(int y=-LY;y<=LY;y++) fflux << y/static_cast<double>(LY) << " " << flow[LY+y]/static_cast<double>(time) << "\n" ;
    fflux.close();

    str = algo+"-mag";
    ofstream fmag(str.c_str());
    for(int y=-LY;y<=LY;y++) fmag << y/static_cast<double>(LY) << " " << mag[LY+y] << "\n" ;
    fmag.close();

    str = algo+"-fourier";
    ofstream ffour(str.c_str());
    for(int x=0;x<LX;x++) ffour << x << " " << fourier[x] << "\n" ;
    ffour.close();

    str = prefix+"/energie";
//    ofstream fenergie(str.c_str(), ios::app);
//    fenergie << ttc << " " << H << " " << F << " " << LX << " " << energie/static_cast<double>(LX) << "\n";
    FILE* fenergie;
    fenergie = fopen(str.c_str(),"a+");
    flock(fileno(fenergie),LOCK_EX);
    fprintf(fenergie,"%f %f %f %i %f \n",ttc,H,F,LX, energie/static_cast<double>(LX) );

    fclose(fenergie);
}

bool isOnlyDouble(const char* str){
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

    return true;
}


void fftw(int* in,double* out,int N,double SIG){

    for(int q = 0; q<N ; q++){
        complex<double> sum = complex<double>(0,0);
        double k    = 2*M_PI*q/N;

        for(int n = 0; n<N ; n++){
            sum += exp(SIG * im * k * static_cast<double>(n)) * static_cast<double>(in[n]) ;
        }
        //        cout << q << " " << sum << " " << norm(sum) << "\n";
        out[q] = (SIG == 1 ) ? norm(sum) : norm(sum)/N ;
    }
}
void fftw(double* in,double* out,int N, double SIG){

    for(int q = 0; q<N ; q++){
        complex<double> sum = complex<double>(0,0);
        double k    = 2*M_PI*q/N;

        for(int n = 0; n<N ; n++){
            sum += exp(SIG * im * k * static_cast<double>(n)) * static_cast<double>(in[n]) ;
        }
        //        cout << q << " " << sum << " " << norm(sum) << "\n";
        out[q] = (SIG == 1 ) ? norm(sum) : norm(sum)/N ;
    }
}

void ft_correl(int* in,double* out,int N,double SIG){
    
    for(int q = 0; q<N ; q++){
        complex<double> sum = complex<double>(0,0);
        double k    = 2*M_PI*q/N;

        for(int n = 0; n<N ; n++){
            for(int m = 0; m<N; m++){
                if( in[n]*in[m] > 0)
                    sum += exp(SIG * im * k * static_cast<double>(n-m));
                else
                    sum -= exp(SIG * im * k * static_cast<double>(n-m));
            }
        }
        //        cout << q << " " << sum << " " << norm(sum) << "\n";
        out[q] = (SIG == 1 ) ? norm(sum) : norm(sum)/ (N^2) ;
    }
}
