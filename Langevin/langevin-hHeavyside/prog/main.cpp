#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <random>
#include <string>
#include <gsl/gsl_histogram.h>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <algorithm> // sort() array


using namespace std;

void parametres(int argc, char* argv[]),
     ecriture(double grille[],double temps, double moyenne,double correl[],gsl_histogram histogramme),
     cfl(double* delta_t, double grille[]);
bool isOnlyDouble(const char* str);
template<typename Type> Type sign(Type a);
template<typename Type> Type modulo(Type a, int b);

std::default_random_engine generator;
std::normal_distribution<double> gaussienne(0.0,1.0);
clock_t cputime = clock();

// Paramètres thermodynamiques
const double    T_C          = 2./log(1+sqrt(2));
double          BETA         = 1/T_C;
double          ttc          = 1;
double          H            = 0;
double          cisaillement = 0;
const double    sigma        = 1;
// Paramètres de la simulation
const double    dx           = 0.1;
const int       LX           = 50;
double          dt           ; // dt est calculé à chaque pas de temps par cfl()
// Paramètres sauvegarde dans fichiers
const long int  NB_PHOTO     = 1e5;
const long int  T_PHOTO      = 1e2;
const long int  T_MAX        = T_PHOTO*NB_PHOTO;
const long int  T_DONNEE     = T_PHOTO*1e0;
const int       histo_min    = -10;
const int       histo_max    = -histo_min;
const int       histo_bin    = 1e2;
string          prefix       = ".";
string          suffix       = "";

/****************** Main ******************/
int main(int argc, char *argv[]){

    generator.seed(time(NULL));
    /********** OPTIONS DEPUIS LE BASH *********/
    parametres(argc,argv); 
    ttc = static_cast<double>(static_cast<int>(ttc*1000))/1000;

    gsl_histogram * histo = gsl_histogram_alloc (histo_bin);
    gsl_histogram_set_ranges_uniform (histo, histo_min,histo_max);

    /***** INITIALISATION DU SYSTÈME ***********/
    long int photo_suiv = 0; 
    double   nb_photo = 1,
             donnee_suiv = 0;
    double   g_correl[LX],
             grille[LX];
    for(int x=0;x<LX;x++){
        g_correl[x] = 0;
        grille[x] = x*dx;
    }
    cfl(&dt,grille);

    for(double t = 0 ; t<=T_MAX ; t+= dt){
        cfl(&dt,grille);
        //        cout << t << " " << dt << "\n";
        if(dt < 1e-8) abort();
        if(photo_suiv <= static_cast<int>(t) && static_cast<int>(modulo(t,T_PHOTO)) == 0) {
            double moyenne = 0,
                   sign_moyenne = 0;
            double correl[LX],
                   n_correl[LX];

            for(int x=0 ; x<LX ; x++){
                moyenne      += grille[x]/LX;
                sign_moyenne += sign(grille[x]-moyenne)/LX;
                correl[x]     = 0;
                n_correl[x]   = 0;
            }
            for(int x=0 ; x<LX ; x++){
                gsl_histogram_increment (histo, grille[x]-moyenne);
                for(int x2=0;x2<LX-x;x2++){
                    correl[x2] += sign(grille[modulo(x2+x,LX)]) * sign(grille[x]) ;
                    n_correl[x2]++;
                }
            }
            for(int x=0 ; x<LX ; x++){
                g_correl[x] = g_correl[x] * (nb_photo-1.) / nb_photo
                    + (correl[x]/n_correl[x]-sign_moyenne*sign_moyenne) / nb_photo;
            }
            nb_photo++;
            photo_suiv += T_PHOTO;

            if(donnee_suiv <=static_cast<int>(t)){ 
                ecriture(grille,t,moyenne,g_correl,*histo);
                donnee_suiv+=T_DONNEE;
            }
        }

        double grille_tmp[LX]; 
        //**        double alea[LX+1];
        /*        for(int i=0 ; i<LX+1 ; i++)
                  if(i==0) alea[i] = 0;
                  else if(i==LX) alea[i] = 1;
                  else alea[i] = gaussienne(generator);
                  std::sort(std::begin(alea), std::end(alea));
                  for(int i=0; i<LX+1; i++)
                  cout << alea[i] << "\n";
                  cout << "----\n";*/

        for(int x=0 ; x<LX ; x++)
            grille_tmp[x] = grille[x]
                + dt / dx * sigma * (grille[modulo(x+1,LX)] - 2*grille[x] + grille[modulo(x-1,LX)])
        //        + dt / dx * cisaillement * grille[x] * (grille[modulo(x+1,LX)] - grille[modulo(x-1,LX)] ) / 2
                + sqrt(2 * dt / BETA) * gaussienne(generator)
                - dt*H*sign(grille[x]); //step model, expected airy
        //                        - dt*H*grille[x] //linear model, expected gaussian

        for(int x=0 ; x<LX ; x++)
            grille[x] = grille_tmp[x];
    }


    gsl_histogram_free(histo);
    return 0;
}
/*************** FIN  **********************/

/********* Calcul du pas de temps dt *******/
void cfl(double* delta_t, double grille[]){
    // dt doit satisfaire à la condition CFL
    // Pour dh/dt + a*h*dh/dx-b*ddh/ddx = 0
    // dt = 1 / ( 2b/dx^2 + max(|h(j+1)-h(j)|)/dx)
    double max = abs(grille[2]-grille[0]);
    for(int i=0 ; i<LX; i++)
        if( abs(grille[modulo(i+1,LX)] - grille[modulo(i-1,LX)]) > max)
            max = abs(grille[modulo(i+1,LX)] - grille[modulo(i-1,LX)]);

    *delta_t = 1 / ( 2 * sigma / dx / dx + cisaillement * max / dx);
    if( *delta_t < 1e-5 ) *delta_t = 1e-5;
}


/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[]){
    if(argc > 1){ 
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "--t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    BETA = 1/(ttc*T_C);
                }else{ 
                    cout << "Le paramètre temperature n'est pas un double"; abort();
                }
            }else if(arg == "--h"){
                if(isOnlyDouble(argv[i+1]))
                    H = atof(argv[i+1]);
                else{ 
                    cout << "Le paramètre h n'est pas un double"; abort();
                }
            }else if(arg == "--g"){
                if(isOnlyDouble(argv[i+1]))
                    cisaillement = atof(argv[i+1]);
                else{
                    cout << "Le paramètre de cisaillement n'est pas un double"; abort();
                }
            }else if(arg == "--prefix")
                prefix = argv[i+1];
            else if(arg == "--suffix")
                suffix = argv[i+1];
        }
    }
}

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(double grille[],double temps, double moyenne, double correl[],gsl_histogram histogramme){
    /*********** CREATION DOSSIER POUR RESULTATS ********/

    if(system(("mkdir -p " + prefix).c_str())) {;}
    string algo = prefix + "/";
    {   
        stringstream stream;
        stream << fixed << setprecision(2) << ttc ;
        algo += "lan_T-" + stream.str();
    }  
    // algo += "-t-" + to_string(t);

    if(fabs(cisaillement) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(2) << cisaillement; 
        algo += "-g-" + stream.str();
    }
    if(fabs(H) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(2) << H; 
        algo += "-h-" + stream.str();
    }
    /******* ÉCRITURES GRANDEURS THERMO ****/
    string str = prefix + "/thermo" + suffix;
    ifstream test(str.c_str());
    ofstream thermo_res;
    if(static_cast<bool>(test)){
        thermo_res.open(str.c_str(), ios::out | ios::app);
    }else{
        thermo_res.open(str.c_str(), ios::out);
    }

    thermo_res << ttc   << " " << temps << " " << float(clock()-cputime)/CLOCKS_PER_SEC << " " << moyenne << "\n";
    thermo_res.close();

    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    str = algo;
    ofstream result(str.c_str());
    str = algo+"-correl";
    ofstream fcorrel(str.c_str());
    str = algo+"-histo";
    FILE* fhisto;
    fhisto = fopen(str.c_str(),"w");
    for(int x=0; x<LX ; x++){
        result << x*dx << " " << grille[x] << "\n";
        fcorrel << x*dx << " " << correl[x] << "\n";
    }
    gsl_histogram_fprintf(fhisto,&histogramme, "%f","%f");
    result.close();
    fcorrel.close();
    fclose(fhisto);

}

template<typename Type> Type modulo(Type a, int b){
    while(a < 0 || a >= b){
        if(a<0) a+= b;
        else if(a >= b) a-=b;
    }
    return a;
}

template<typename Type> Type sign(Type a){
    return (a>0)-(a<0);
}

bool isOnlyDouble(const char* str){
    char* endptr = 0;
    strtod(str, &endptr);
    return (*endptr == '\0' || endptr != str);
}
