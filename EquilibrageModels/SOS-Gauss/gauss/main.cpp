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

//Modules diagonalisation de matrices
using namespace std;
const double T_C    = 2./log(1.+sqrt(2)), J = 1;
const double Mmin = 0, Mmax = 3;
const int Mn = 12;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 100;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 2e6;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ);
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

    int Nb_T = Mn/nbprocs;

    /***** Parallélisation sur les champs magnétiques ****/
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        double mu =  normspace(k,Mmin,Mmax,Mn);
        Beta = 1;
        /******** Initialisation du systeme *****/
        //Mise en place du vecteur des positions
        int* system = new int[LX]; 
        //
        //Début à T=0
        for(int x=0;x<LX;x++)
            system[x]  = LY/2;
        // Équilibrage
        int ajout,tirage;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,mu)){
                system[tirage]+=ajout;
            }
        }
        double phi = 0;
        for(int x=0;x<LX;x++){
            phi+=system[x];
        }
        double mag = 0;
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,mu)){
                phi += ajout;
                system[tirage]+=ajout;
            }
            mag+= phi;
            /***** Mesure du système *******/
        }
        /** Fonction d'autocorrélation **/
        mag /=static_cast<double>(LX*T_MAX*LX);
        cout << mu << " " << mag << endl;
    }
    MPI_Finalize();

return 0;
}
/********** Fin main ***********/

bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<0 or hx2>=LY) return false;
    double champ = Champ*(abs(hx2)-abs(hx));
    //double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );
    double sos    = J*( pow(hx2-hxm,2) + pow(hx2-hxp,2) - (pow(hx - hxm,2)  +  pow(hx - hxp,2) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}
