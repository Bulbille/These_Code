#include "init.h"
#include "math.h"
#include "dynamique.h"
#include <gsl/gsl_histogram.h>

void parametres(int argc, char* argv[]);
void ecriture(int grille[][L_Y], long int temps, int* histogramme,double magy[]);

double BETA;
double ttc = 1;
string prefix = ".",suffix="";
string algo;
double H=H0;
double w=w0;

std::default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

//0 → L-1 : périodique
//1 → L-2 : fixe
std::uniform_int_distribution<int> rand_lx(0,L_X-1);
std::uniform_int_distribution<int> rand_ly(1,L_Y-2);

clock_t cputime = clock();


/****************** Main ******************/
int main(int argc, char *argv[]){

    int grille[L_X][L_Y];
    generator.seed(time(NULL));
    /********** OPTIONS DEPUIS LE BASH *********/
    /* Le programme peut accepter deux options : */
    /* --beta xxx */
    /* --reprendre fichier */
    /* La première option permet de faire la simulation avec une température donnée */
    /* La deuxième permet de reprendre la simulation à partir d'un fichier */
    /**/    parametres(argc,argv); /**/
    /**/    generation(grille);
    /**/    ttc = static_cast<double>(static_cast<int>(ttc*1000))/1000;
    /******************************************/

    /******* VARIABLES MICRO *********/
    int* histo = new int[L_Y+1];
    for(int y=0; y<L_Y+1; y++) histo[y] = 0;

    /**** VARIABLES INTERMÉDIAIRES *****/
    int nb_liens;
    int liens[L_X * L_Y * 5];
    double acceptance_rate[L_X * L_Y * 2];
    double delta_t=0;
    maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens);
    long int photo_suiv=0; int nb_photos=0;double donnee_suiv=0;

    double magy[L_Y]; 
    for(int y=0;y<L_Y;y++) magy[y]=0;

    /********************************************/
    /****** DYNAMIQUE DE KAWASAKI ***************/
    /*******************************************/

    double Hprime = H;
    H = H0;
    double wprime = w;
    w = w0;
    for(double t = 0;t <= T_CHAMP; t+= delta_t){
        maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens);
        int tirage = tirage_lien(acceptance_rate, nb_liens);
        int x1 = liens[tirage*4+0];
        int y1 = liens[tirage*4+1];
        int x2 = liens[tirage*4+2];
        int y2 = liens[tirage*4+3];
        grille[x1][y1] *= -1;
        grille[x2][y2] *= -1;
    }
    H=Hprime;
    w=wprime;

    for(double t = 0; true ; t+= delta_t){
        /******************************************************/
        /*********** CALCULS GRANDEURS THERMO *****************/
        /** moyenne_n = moyenne_(n-1) * (n-1) / n + a_n / n ***/
        /******************************************************/

        int modulo_photo = static_cast<int>(modulo_d(delta_t+t,T_PHOTO));
        if(photo_suiv <= static_cast<int>(t) && modulo_photo == 0) {
            photo_suiv += T_PHOTO; nb_photos++;
            // L'interface est par construction à L_X/2
            /*
            for(int x = 0;x<L_X ; x++){
                int mag_x = 0;
                for(int y = 0; y<L_Y; y++)
                    mag_x += grille[x][y];
                histo[L_Y/2+mag_x]++;
            }

            /*** Grandeurs thermo vectorielles ***
            double p_magy[L_Y];
            for(int y=0;y<L_Y;y++){
                p_magy[y] = 0; 
                // Magnetisation et énergies de bandes
                for(int x=0;x<L_X;x++){
                    p_magy[y] += grille[x][y];
                }
            }
            for(int y=0;y<L_Y;y++){
                magy[y] = 1.*magy[y]*(nb_photos-1)/nb_photos+1.*p_magy[y]/L_X/nb_photos;
            }*/

            if(donnee_suiv <=static_cast<int>(t)){
                ecriture(grille,t,histo,magy);
                donnee_suiv+=T_DONNEE;
            }

        }
        /****** MONTE CARLO STEP ****/

        maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens);

        int tirage = tirage_lien(acceptance_rate, nb_liens);
        int x1 = liens[tirage*4+0];
        int y1 = liens[tirage*4+1];
        int x2 = liens[tirage*4+2];
        int y2 = liens[tirage*4+3];
        grille[x1][y1] *= -1;
        grille[x2][y2] *= -1;
        //*/
    }

    /*************** FIN  **********************/
    return 0;
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[]){
    BETA=1./T_C;
    if(argc > 1){ 
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "-t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    BETA = 1/(ttc*T_C);
                }
                else{ cout << "Le paramètre temperature n'est pas un double"; abort();}
            }
            else if(arg == "-h"){
                if(isOnlyDouble(argv[i+1])){
                    H = atof(argv[i+1]);
                }
                else{ cout << "Le paramètre h n'est pas un double"; abort();}
            }
            else if(arg == "-f"){
                if(isOnlyDouble(argv[i+1])){
                    w = atof(argv[i+1]);
                }
                else{ cout << "Le paramètre w n'est pas un double"; abort();}
            }
            else if(arg == "-prefix"){
                prefix = argv[i+1];
            }
            else if(arg == "-suffix"){
                suffix = argv[i+1];
            }
        }
    }

    if(system(("mkdir -p " + prefix).c_str())) {;}
    algo = prefix + "/";
    {
        stringstream stream;
        stream << fixed << setprecision(2) << ttc ;
        algo += "kaw_T-" + stream.str();
    }
    if(fabs(H) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(4) << H;
        algo += "-h-" + stream.str();
    }
    if(fabs(w) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(2) << w;
        algo += "-f-" + stream.str();
    }
}

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(int grille[][L_Y], long int temps, int* histogramme,double magy[]){

    /******* ÉCRITURES GRANDEURS THERMO ****/
    string str = prefix + "/thermo_kaw" + suffix;
    FILE* fthermo;
    fthermo = fopen(str.c_str(),"a+");
    flock(fileno(fthermo),LOCK_EX);
    fprintf(fthermo,"%f %f %f \n",ttc,H,w);
    fclose(fthermo);

    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    str = algo+"-time-"+to_string(temps);
    ofstream result(str.c_str());
    str = algo+"-magy";
    ofstream fmagy(str.c_str());
    for(int y=0; y<L_Y;y++)
    {
        fmagy << (1+2*y)/static_cast<double>(2*L_Y) << " " << magy[y] << "\n";
        for(int x=0; x<L_X;x++){
            result << x << " " << y << " " << grille[x][y] << "\n";
        }//*/
    }
    fmagy.close();
    result.close();
    str = algo+"-histo";
    ofstream fhisto(str.c_str());
    for(int y=0;y<=L_Y;y+=2) fhisto << y-L_Y/2 << " " << histogramme[y] << "\n" ;
    fhisto.close();

}
