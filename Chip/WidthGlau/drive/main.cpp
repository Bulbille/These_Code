/******************************************************************/
// Plan SemiInfini, dynamique de Kawasaki et de Glauber
// Modèle Chipping, hamiltonien SOS
// Parallélisé sur les cisaillements
/******************************************************************/



//Basic stuff
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
const int       bitLX = 11, LX = 2<<bitLX;
const int LY = 4000;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e6, T_MAX = 1e5, N_MAX = 4;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ);

/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    string str;

    double* ecarts = new double[T_MAX];
    double* valeurs = new double[T_MAX];
    for(int t=0;t<T_MAX;t++)
        ecarts[t] = 0;
    int* system = new int[LX]; 
    int ajout,tirage;
    int hx;
    double phi = 0, phi2 = 0;
       
    for(int n = 0; n < N_MAX;n++){
        for(int x=0;x<LX;x++)
            system[x]  = 0;
        phi = 0; phi2 = 0;
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,0)){
                hx = system[tirage];
                phi += ajout;
                phi2 += 2*hx*ajout+1;
                system[tirage]+=ajout;
            }
            /***** Mesure du système *******/
            if(t%LX == 0){
                ecarts[static_cast<int>(t/LX)] += (phi2/LX)-pow(phi/LX,2);
                valeurs[static_cast<int>(t/LX)] += phi/LX;
            }
        }
    }
    str = prefix+"data"+to_string(N_MAX);
    ofstream fdata(str.c_str(),std::ofstream::out);
    for(int t=0;t<T_MAX;t++)
        fdata << t << " " << sqrt(ecarts[t]/N_MAX) << " " << valeurs[t] << endl;
    fdata.close();
    delete[] system;
    delete[] ecarts;

return 0;
}


bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<=-LY or hx2>=LY) return false;
    double champ = Champ*(abs(hx2)-abs(hx))/2;
    double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}


