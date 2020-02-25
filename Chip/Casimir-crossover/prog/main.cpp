//Basic stuff
#include <iostream> //standard Input-Output
#include <fstream>  //standard Files
#include <sys/file.h> //Files lock
#include <iomanip>  //setprecision
#include <sstream>  //stringstream
#include <string>   //strings
//Mathematics
#include <cmath>
#include <random>

#include <omp.h>

using namespace std;
//Propriétés de l'Hamiltonien
const double T_C    = 2./log(1.+sqrt(2)), J = 1.0;
double       Beta,ttc,H,F;
string prefix = ".",suffix="";
string algo;
//Propriétés de la grille
int       LX    = 50, LY = 10, N = LX*LY;
int       NLambda = 40, K0 = 0; 
const int T_EQ  = 1e2;
int       T_PHO = 1.0e1;
const double    error = 1e-11;
//const long int T_MAX = T_PHO*1.2e8;
//Générateur nombres aléatoires
default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

/***********************************/
/**** Définitions des fonctions ****/
int modulo(int a, int b); 
double modulo_d(double a, int b); 
void parametres(int argc, char* argv[]);
void ecriture(double energie);
bool isOnlyDouble(const char* str);

bool Esos(int* array, int x, int ajout,double kbeta,double Champ,double couplage){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)];
    int hk,hkm,hkp,y;

    double Hcouplage = 0;
    int tabi [3] = {hx,hxp,hxm};
    for(int i=0;i<3;i++){
        y = tabi[i];
        hk  = (y >= K0) -   (y < K0);
        hkm = (y >= K0-1) - (y < K0+1);
        hkp = (y >= K0+1) - (y < K0-1);
        Hcouplage += hkm*hkp - hk * (hkm + hkp);
    }
    Hcouplage *= (1-couplage);

    double champ = Champ*( abs(hx+ ajout) - abs(hx) );
//    double shear = - J*F*sens*(hx+hxp)/2;

    double sos    = abs(hx + ajout - hxm) + abs(hx+ajout - hxp) - (abs(hx - hxm)  +  abs(hx - hxp) );
    double D_e   = (hxp < LY and hx > -LY) * exp(-kbeta*(sos+champ+couplage));
//    double gauss  = pow(hx + ajout - hxm,2) + pow(hx+ajout - hxp,2) -(pow(hx - hxm,2)  +  pow(hx - hxp,2) );
//    double D_e   = (hx + ajout <= LY and hx + ajout >= -LY) * exp(-kbeta*(gauss+champ+Hcouplage));
    return (rand_01(generator) <  D_e );
}
/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);
    std::uniform_int_distribution<int> rand_lx(0,LX-1);
    generator.seed(time(NULL));
    double Free_energy [2] = {0,0};
    double dh = 1./(NLambda-1);

/*** Boucle pour calculer F(L) et F(2*L) ****/
    for(int longueur = 0 ; longueur <= 1 ; longueur++){
        if(longueur == 1) LY*= 2;

    /**** Boucle pour calculer F(L,lambda) ****/
#pragma omp parallel for schedule(dynamic) 
        for(int iter_lambda = 0; iter_lambda < NLambda; iter_lambda++){

            int* system = new int[LX]; 
            double energie = 0, e_tmp = 1, H1_H0;
            int hk,hkm,hkp,hx;
            int photo_suiv = 0, nb_photos = 0;
            long int t;
            for(int x = 0;x<LX;x++) system[x] = 0;

            for(t = 0; abs(e_tmp-energie) > error ; t++){

                /***** Étape de Monte Carlo *******/
                int tirage  = rand_lx(generator);
                int ajout   = -1+2*round(rand_01(generator));

                if(Esos(system,tirage,ajout,Beta,H,static_cast<double>(iter_lambda)/(NLambda-1)))
                    system[tirage] += ajout;

                /***** Mesure du système *******/
                if(photo_suiv <= t && (t-1)%T_PHO == 0 && t >= T_EQ) {
                    photo_suiv += T_PHO; nb_photos++;

                    //*
                    H1_H0 = 0;
                    for(int x = 0; x<LX;x++){
                        hx = system[x];
                        hk  = (hx >= K0) - (hx < K0);
                        hkm = (hx >= K0-1) - (hx < K0-1);
                        hkp = (hx >= K0+1) - (hx < K0+1);
                        H1_H0 -= hkm*hkp - hk * (hkm + hkp);
                    }
                    e_tmp = energie;
                    energie  = energie*(nb_photos-1)/nb_photos + H1_H0 /static_cast<double>(LX)   /nb_photos;
                    cout << energie << " " << H1_H0 /static_cast<double>(LX)<< "\n";
                }
            }
#pragma omp critical
            {
            if(iter_lambda == 0 or iter_lambda == NLambda-1)
                Free_energy[longueur] += energie;
            else if (iter_lambda %2 == 0 )
                Free_energy[longueur] += 2*energie;
            else
                Free_energy[longueur] += 4*energie;
            }
        /*** Fin boucle sur les lambda ****/
        }
        Free_energy[longueur] *= dh/3.;
    /**** Fin boucle sur les longueurs ****/
    }
    cout << Free_energy[0] << " " << Free_energy[1] << " " << Free_energy[0]- Free_energy[1] << "\n";

    ecriture(Free_energy[0]-Free_energy[1]);
//    auto debut = std::chrono::high_resolution_clock::now();
//    auto fin = std::chrono::high_resolution_clock::now();
//    auto duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
//    std::cout << "Duree : " << duree/1e9 <<  std::endl;
return 0;
}
/***************************************************************/
/************** Fin du main ************************************/
/***************************************************************/

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(double energie){
    /*********** CREATION DOSSIER POUR RESULTATS ********/
    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    string str = algo;
    str = prefix+"/energie";
    FILE* fenergie;
    fenergie = fopen(str.c_str(),"a+");
    flock(fileno(fenergie),LOCK_EX);
    fprintf(fenergie,"%.2f %f %i %f \n",ttc,H,LY, energie);

    fclose(fenergie);
}

bool isOnlyDouble(const char* str){
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

    return true;
}

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
    bool taille_x = false;
    bool taille_y = false;
    if(argc > 1){
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "-t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    Beta = 1/(ttc*T_C);
                }
            }
            else if(arg == "-h"){
                if(isOnlyDouble(argv[i+1]))
                    H = atof(argv[i+1]);
            }
            else if(arg == "-prefix"){
                prefix = argv[i+1];
            }
            else if(arg == "-f"){
                if(isOnlyDouble(argv[i+1]))
                    F = atof(argv[i+1]);
            }
            else if(arg == "-lx"){
                if(isOnlyDouble(argv[i+1])){
                    LX = round(atof(argv[i+1]));
                    taille_x = true;
                }
            }
            else if(arg == "-ly"){
                if(isOnlyDouble(argv[i+1])){
                    LY = round(atof(argv[i+1]));
                    taille_y = true;
                }
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
    if(taille_x){
        stringstream stream;
        stream << fixed << setprecision(3) << LX;
        algo += "-lx-" + stream.str();
    }
    if(taille_y){
        stringstream stream;
        stream << fixed << setprecision(3) << LY;
        algo += "-ly-" + stream.str();
    }

}
