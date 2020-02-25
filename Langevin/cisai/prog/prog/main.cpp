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
     ecriture(double* grille,double temps, double moyenne,gsl_histogram histogramme),
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
const int       LX           = 500;
double          dt           = 0.01; // dt est calculé à chaque pas de temps par cfl()
// Paramètres sauvegarde dans fichiers
const long int  T_PHOTO      = 1e1;
const long int  T_DONNEE     = T_PHOTO*1e3;
const int       histo_min    = -5e2;
const int       histo_max    = -histo_min;
const int       histo_bin    = histo_max*1e1;
string          prefix       = ".";
string          suffix       = "";
string          algo;

/****************** Main ******************/
int main(int argc, char *argv[]){

    generator.seed(time(NULL));
    /********** OPTIONS DEPUIS LE BASH *********/
    parametres(argc,argv); 
    ttc = static_cast<double>(static_cast<int>(ttc*1000))/1000;


    /***** INITIALISATION DU SYSTÈME ***********/
    long int photo_suiv = 0; 
    double   nb_photo = 1,
             donnee_suiv = 0;
    double* grille = new double[LX];
    gsl_histogram * histo = gsl_histogram_alloc (histo_bin);
    gsl_histogram_set_ranges_uniform (histo, histo_min,histo_max);
    for(int x=0;x<LX;x++) grille[x] = 0;

    for(double t = 0 ; true; t+= dt){
//        cfl(&dt,grille);

        if(photo_suiv <= static_cast<int>(t) && static_cast<int>(modulo(t,T_PHOTO)) == 0) {
            double moyenne = 0;

            for(int x=0 ; x<LX ; x++){
                moyenne      += grille[x]/LX;
                gsl_histogram_increment (histo, grille[x]-moyenne);
            }
            
            nb_photo++;
            photo_suiv += T_PHOTO;

            if(donnee_suiv <=static_cast<int>(t)){ 
            cout << t << "\n";
                ecriture(grille,t,moyenne,*histo);
                donnee_suiv+=T_DONNEE;
            }
        }

        double* grille_tmp = new double[LX]; 
        double contrainte = 0;

        for(int x=0 ; x<LX ; x++){
            double gauss = gaussienne(generator);

            grille_tmp[x] = grille[x]
                + dx * sigma * (grille[modulo(x+1,LX)] - 2*grille[x] + grille[modulo(x-1,LX)])
                + dt / dx * cisaillement * sign(grille[x]) * (grille[modulo(x+1,LX)] - grille[modulo(x-1,LX)] ) / 2
                + sqrt(2 * dt / BETA) * gauss
                - H*sign(grille[x]); //step model, expected airy
        //                        - dt*H*grille[x]; //linear model, expected gaussian

            contrainte +=grille_tmp[x];
        }

        for(int x=0 ; x<LX ; x++)
            grille[x] = grille_tmp[x]- 1./LX*(contrainte);
    }

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
    *delta_t = dx / sigma;
}


/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[]){
    if(argc > 1){ 
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "-t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    BETA = 1/(ttc*T_C);
                }else{ 
                    cout << "Le paramètre temperature n'est pas un double"; abort();
                }
            }else if(arg == "-h"){
                if(isOnlyDouble(argv[i+1]))
                    H = atof(argv[i+1]);
                else{ 
                    cout << "Le paramètre h n'est pas un double"; abort();
                }
            }else if(arg == "-g"){
                if(isOnlyDouble(argv[i+1]))
                    cisaillement = atof(argv[i+1]);
                else{
                    cout << "Le paramètre de cisaillement n'est pas un double"; abort();
                }
            }else if(arg == "-prefix")
                prefix = argv[i+1];
            else if(arg == "-suffix")
                suffix = argv[i+1];
        }
    }
    if(system(("mkdir -p " + prefix).c_str())) {;}
    algo = prefix + "/";
    stringstream stream;
    stream << fixed << setprecision(2) << ttc ;
    algo += "lan_T-" + stream.str();
    stream.str("");
    stream << fixed << setprecision(2) << H; 
    algo += "-h-" + stream.str();
    stream.str("");
    stream << fixed << setprecision(2) << cisaillement; 
    algo += "-g-" + stream.str();
}

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(double* grille,double temps, double moyenne,gsl_histogram histogramme){
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
    str = algo+"-histo";
    FILE* fhisto;
    fhisto = fopen(str.c_str(),"w");
    for(int x=0; x<LX ; x++){
        result << x*dx << " " << grille[x] << "\n";
    }
    gsl_histogram_fprintf(fhisto,&histogramme, "%f","%f");
    result.close();
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
