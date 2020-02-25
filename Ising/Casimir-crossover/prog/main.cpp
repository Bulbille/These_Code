#include "init.h"
#include "math.h"
#include "dynamique.h"
#include <gsl/gsl_histogram.h>

void parametres(int argc, char* argv[]);
void ecriture(int grille[][L_Y], int* histogramme,double magy[]);

double BETA;
double ttc = 1;
string prefix = ".",suffix="";
string algo;
double H=H0, w=w0;
const double NLambda = 40
const double error = 1e-11;

std::default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

//0 → L-1 : périodique
//1 → L-2 : fixe
std::uniform_int_distribution<int> rand_lx(0,L_X-1);
std::uniform_int_distribution<int> rand_ly(1,L_Y-2);

clock_t cputime = clock();


/****************** Main ******************/
int main(int argc, char *argv[]){

    generator.seed(time(NULL));
    parametres(argc,argv); 
    double Free_energy [2] = {0,0};
    double dh = 1./(NLambda-1);


    /********************************************/
    /****** DYNAMIQUE DE KAWASAKI ***************/
    /*******************************************/
for(int longueur = 0 ; longueur <= 1 ; longueur++){
    if(longueur == 1) LY*= 2;

#pragma omp parallel for schedule(dynamic) 
    for(int iter_lambda = 0; iter_lambda < NLambda; iter_lambda++){
        double energie = 0, e_tmp = 1, H1_H0; 

        int grille[L_X][L_Y], liens[L_X * L_Y * 5];
        double acceptance_rate[L_X * L_Y * 2];
        generation(grille);
        long int photo_suiv=0,nb_photos=0,t; 
        double donnee_suiv=0,delta_t;
        int nb_liens;
        maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens, static_cast<double>(iter_lambda)/(NLambda-1));

        for(t = 0; abs(e_tmp-energie) > error ; t+= delta_t){

            maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens, static_cast<double>(iter_lambda)/(NLambda-1));
            int tirage = tirage_lien(acceptance_rate, nb_liens);
            grille[liens[tirage*4+0]][liens[tirage*4+1]] *= -1; 
            grille[liens[tirage*4+2]][liens[tirage*4+3]] *= -1;

            if(photo_suiv <= t && (static_cast<int>(t)-1)%T_PHO == 0 && t >= T_EQ) { 
                photo_suiv += T_PHO; nb_photos++; 
                H1_H0 = 0;
                for(int x = 0; x<LX;x++){
                    int hk,hkm,hkp;
                    hk  = grille[x][K0+0];
                    hkm = grille[x][K0-1];
                    hkp = grille[x][K0+1];
                    H1_H0 -= hkm*hkp - hk * (hkm + hkp);
                }
                e_tmp = energie;
                energie  = energie*(nb_photos-1)/nb_photos + H1_H0 /static_cast<double>(LX)   /nb_photos;
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


    /*************** FIN  **********************/
    return 0;
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
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
}

void ecriture(double energie){
    /*********** CREATION DOSSIER POUR RESULTATS ********/
    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    str = prefix+"/energie";
    FILE* fenergie;
    fenergie = fopen(str.c_str(),"a+");
    flock(fileno(fenergie),LOCK_EX);
    fprintf(fenergie,"%.2f %f %i %f \n",ttc,H,LY, energie);

    fclose(fenergie);
}
