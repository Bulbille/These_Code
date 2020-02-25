#include <iostream>
#include <fstream>
#include <random>
#include <sys/file.h> //Files lock
#include <sstream>
#include <iomanip>
using namespace std;
/*** Constantes ***/
const double T_C = 2./log(1+sqrt(2)),J=1;
double BETA,H,w,F,ttc;
long int T_MAX = 1e8;
const int TAILLE_X =100,TAILLE_Y = 20;
string prefix = ".",algo;

std::default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);
std::uniform_int_distribution<int> rand_lx(0,TAILLE_X-1);
std::uniform_int_distribution<int> rand_ly(1,TAILLE_Y-2);

void parametres(int argc, char* argv[]);
void ecriture(int grille[][TAILLE_Y], double mag_y[], double e_x[], double e_y[], double magnetisation, double e_int);
void generation(int grille[][TAILLE_Y]);
int modulo(int a, int b);
bool isOnlyDouble(const char* str);

/****************** Main ******************/
int main(int argc, char *argv[]){

    int grille[TAILLE_X][TAILLE_Y];
    generator.seed(time(NULL));
    parametres(argc,argv); 
    generation(grille);

    /******* VARIABLES THERMO *********/
    double mag=0, energie=0;
    double* magy   = new double[TAILLE_Y];
    double* ex     = new double[TAILLE_Y];
    double* ey     = new double[TAILLE_Y];

    double** magypar = new double*[TAILLE_Y];
    double** expar   = new double*[TAILLE_Y];
    double** eypar   = new double*[TAILLE_Y];

    for(int y=0;y<TAILLE_Y;y++){
        magy[y]=0;
        ex[y]=0;
        ey[y]=0;
        magypar[y]   = new double [2]; 
        expar[y]     = new double [2]; 
        eypar[y]     = new double [2]; 
        for(int x=0;x<TAILLE_X;x++){
            magypar[y][0] += grille[x][y]; 
            expar[y][0]   -= grille[x][y]*grille[modulo(x+1, TAILLE_X)][y]; 
            if(y>0)
                eypar[y][0]   -= grille[x][y]*grille[x][y-1];
            if(y<TAILLE_Y-1)
                eypar[y][0]   -= grille[x][y]*grille[x][y+1]; 
        }
        magypar[y][1] = 0;
        expar[y][1]   = 0;
        eypar[y][1]   = 0;    
    }


    /********************************************/
    /*******************************************/
    double t, delta_t = 1; // delta_t = 1 pour Glauber

    for(t = 0; t<T_MAX ; t+= delta_t){

        /****** MONTE CARLO STEP ****/
        /**** Kawasaki **/
        /*** Glauber ****/
        int x = rand_lx(generator);
        int y = rand_ly(generator);
        double e_flip = grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][y-1]+grille[x][y+1];
        e_flip *= 2*J*grille[x][y];
        double mag =  H*grille[x][y]* ( 2*(y < TAILLE_Y/2) -1) ;
        if(rand_01(generator) < exp(-BETA*(e_flip+mag))){
            grille[x][y] *= -1;

            /****** CALCULS GRANDEURS THERMO **************/
            // Énergie verticale
            eypar[y][0] -= 2*grille[x][y]*(grille[x][y-1]+grille[x][y+1]);
            ey[y]       += expar[y][0] * (t-eypar[y][1]);
            eypar[y][1] = t;

            eypar[y+1][0] -= 2*grille[x][y]*grille[x][y+1];
            ey[y+1]       += expar[y+1][0] * (t-eypar[y+1][1]);
            eypar[y+1][1] = t;

            eypar[y-1][0] -= 2*grille[x][y]*grille[x][y-1];
            ey[y-1]       += expar[y-1][0] * (t-eypar[y-1][1]);
            eypar[y-1][1] = t;

            // Énergie horizontale
            expar[y][0] -= 2*grille[x][y]*( grille[modulo(x+1, TAILLE_X)][y]+ grille[modulo(x-1, TAILLE_X)][y]);
            ex[y]       += expar[y][0] * (t-expar[y][1]);
            expar[y][1] = t;

            magypar[y][0] += 2*grille[x][y];
            magy[y]       += magypar[y][0] * (t-magypar[y][1]);
            magypar[y][1] = t;
        }
    }
    /**** Fin simulation MC ***/
    for(int y=0;y<TAILLE_Y;y++){
        ex[y]    += expar[y][0]    * (t-expar[y][1]);      ex[y]   /= static_cast<double>(TAILLE_X)*t;
        ey[y]    += eypar[y][0]    * (t-eypar[y][1]);      ey[y]   /= static_cast<double>(TAILLE_X)*t;
        magy[y]  += magypar[y][0]  * (t-magypar[y][1]);    magy[y] /= static_cast<double>(TAILLE_X)*t;
        cout <<y << " | " <<  " " << ey[y] << " " << eypar[y][0] << "\n";
    }
    ecriture(grille,magy,ex,ey,mag,energie);
    return 0;
}

/*for(int y=0;y<TAILLE_Y;y++){
    for(int x=0;x<TAILLE_X;x++){
        //mag_conf += grille[x][y]; 
        magy[y] += grille[x][y] * delta_t;
        ex[y]   -= grille[x][y]*grille[modulo(x+1, TAILLE_X)][y]*delta_t; 
        if(y < TAILLE_Y-1){ 
            ey[y] -= grille[x][y]*grille[x][y+1]              *delta_t;
        }
    }
}*/

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(int grille[][TAILLE_Y], double mag_y[], double e_x[], double e_y[], double magnetisation, double e_int){
    string str;
    /******* ÉCRITURES GRANDEURS THERMO ****
    str = prefix + "/thermo_kaw";
    FILE* fenergie;
    fenergie = fopen(str.c_str(),"a+");
    flock(fileno(fenergie),LOCK_EX);
    fprintf(fenergie,"%f %f %f %ld \n",ttc,magnetisation,e_int,temps);
    fclose(fenergie);
    //*/
    str = algo;
    ofstream result(str.c_str());
    /* Magnetisation */ 
    str = algo +  "-mag_y";
    ofstream res_mag_y(str.c_str());
    /**** ENERGIE DE BANDE **/
    str = algo + "-energie";
    ofstream energie_bande(str.c_str());

    for(int y=0; y<TAILLE_Y;y++){
        res_mag_y <<  y << " " << mag_y[y] << "\n";
        energie_bande <<  y << " " << e_x[y] << " " << e_y[y]<< "\n";
        for(int x=0; x<TAILLE_X;x++)
            result << x << " " << y << " " << grille[x][y] << "\n";
    }
    result.close();
    res_mag_y.close();
    energie_bande.close();
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[]){
    ttc=1;
    BETA=1./ttc;
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
            else if(arg == "-w"){
                if(isOnlyDouble(argv[i+1]))
                    w = atof(argv[i+1]);
                else{ cout << "Le paramètre w0 n'est pas un double"; abort();}
            }
            else if(arg == "-f"){
                if(isOnlyDouble(argv[i+1]))
                    F = atof(argv[i+1]);
                else{ cout << "Le paramètre F0 n'est pas un double"; abort();}
            }
            else if(arg == "-prefix")
                prefix = argv[i+1];
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
        algo += "-w-" + stream.str();
    }   
    if(fabs(F) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(2) << F;
        algo += "-f-" + stream.str();
    }   
}

void generation(int grille[][TAILLE_Y]){
    /* Configuration initiale */
    for(int j=0;j<TAILLE_Y;j++){
        for(int i=0;i<TAILLE_X;i++){
            if(TAILLE_Y % 2 ==0){
                if(j>=TAILLE_Y/2)
                    grille[i][j] = -1;
                else
                    grille[i][j] = 1;
            }
            else{
                if(j>TAILLE_Y/2)
                    grille[i][j] = -1;
                else if(j == TAILLE_Y/2 && i < TAILLE_X/2)
                    grille[i][j] = -1;
                else
                    grille[i][j] = 1;
            }    
        }
    }

}
/******* Fonction modulo ********/
int modulo(int a, int b){
    while(a < 0 or a >= b){
        if(a<0) a+= b;
        else if(a >= b) a-=b;
    }
    return a;
}
/**** VERIFICATION DOUBLE *****/
bool isOnlyDouble(const char* str)
{
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

    return true;
}

