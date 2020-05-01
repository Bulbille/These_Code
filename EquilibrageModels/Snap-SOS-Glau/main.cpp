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
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 40;
const double H = 2;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e5;

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ);
double normspace(int step,double min,double max, double n){
    return (n == 1) ? max : min+(max-min)/(n-1)*1.*step;
}

#define MODEL 'A'


/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);
    string str;

    int* system = new int[LX]; 

    for(int x=0;x<LX;x++)
        system[x]  = LY/2;

    int ajout,tirage;
    for(long int t = 0; t<LX*T_EQ ; t++){
        tirage = rand_lx(generator);
        ajout   = 2*r01int(generator)-1;
        if(EsosGlau(system,tirage,ajout,Beta,H))
            system[tirage]+=ajout;
    }
    str = prefix+"/"+MODEL+to_string(H);
    ofstream feq(str.c_str(),std::ofstream::out);
    for(int x=0;x<LX;x++)
        feq << x << " " <<  system[x] << endl;
    feq.close();

return 0;
}
/********** Fin main ***********/

bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<0 or hx2 >= LY) return false;
#if MODEL == 'A'
    double champ = Champ*ajout;
#elif MODEL == 'B'
    double champ = -Champ*(abs(hx2-LY/2)-abs(hx-LY/2));
#endif
    double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}
