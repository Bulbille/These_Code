/******************************************************************/
// Plan SemiInfini, dynamique de Kawasaki et de Glauber
// Modèle Chipping, hamiltonien SOS
/******************************************************************/

/*** Calcule l'énergie du systeme libre en fonction de la temperature T **/

//Basic stuff
#include <mpi.h>
#include <sys/file.h> //Files lock
#include <random> 
#include <iostream> //cout
#include <sstream> //stringstream
#include <iomanip> // setprecision
#include <fstream> //files

#define MODEL "CHIM"

//Modules diagonalisation de matrices
using namespace std;
const double T_C    = 2./log(1.+sqrt(2)), J = 1;
const double Hmin =0.01, Hmax = 1;
const int Hn = 20;
double  Beta = 1, ttc;
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 20;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e7;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<0 or hx2>LY) return false;
    double champ = -Champ*(abs(hx2-LY/2)-abs(hx-LY/2));
//    double champ = Champ*(hx2-hx);
    double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

double normspace(int step,double min,double max, double n){
    return (n == 1) ? max : min+(max-min)/(n-1)*1.*step;
}


/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
    ostringstream streamT,streamM;
    string str;

    int Nb_H= Hn/nbprocs;
    double* valeurs = new double[Nb_H];
    double* valTot = new double[Hn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_H; k<(rang+1)*Nb_H ; k++){
        cout << k << endl;
        valeurs[cmpt] = 0;
        double H =  normspace(k,Hmin,Hmax,Hn);
        /******** Initialisation du systeme *****/
        int* system = new int[LX]; 
        for(int x=0;x<LX;x++)
            system[x]  = LY/2;


        // Équilibrage
        int ajout,tirage;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;
            if(EsosGlau(system,tirage,ajout,Beta,H))
                system[tirage]+=ajout;
        }
        //Calcul des valeurs initiales        
        int phi = 0;
        for(int x=0;x<LX;x++){
            phi += abs(system[x]-LY/2);
//            phi += system[x];            
        }  
       
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,H)){
//                phi += ajout;
                phi += abs(system[tirage]+ajout-LY/2)-abs(system[tirage]-LY/2);
                system[tirage]+=ajout;
            }
            valeurs[cmpt] += phi;
        }
        valeurs[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_H,MPI_DOUBLE,valTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/"+MODEL+"L"+to_string(LY);
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Hn;i++)
            fmag << normspace(i,Hmin,Hmax,Hn)<< "\t" << valTot[i] <<  "\n";
        fmag.close();
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] valTot;
return 0;
}
/********** Fin main ***********/

