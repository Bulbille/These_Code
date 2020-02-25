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
const int LY = 400;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const int  Mun = 100;
const long int T_EQ = 1e6, T_MAX = 1e6, N_MAX = 1;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f);

/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    string str;

    double* ecarts = new double[T_MAX];
    for(int t=0;t<T_MAX;t++)
        ecarts[t] = 0;
    int* system = new int[LX]; 
    for(int x=0;x<LX;x++)
        system[x]  = 0;
    int sens,tirage;
    int hx,hn,hnp,hxp;
    double phi = 0, phi2 = 0;
       
    for(int n = 0; n < N_MAX;n++){
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;

            if(EsosKaw(system,tirage,sens,Beta,0,0)){
                hx = system[tirage];
                hxp  = system[modulo(tirage+sens,LX)]; 
                phi2 += 2*(hxp-hx)+2;
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
            /***** Mesure du système *******/
            if(t%LX == 0){
                ecarts[static_cast<int>(t/LX)] += (phi2/LX);
            }
        }
    }
    str = prefix+"data"+to_string(N_MAX);
    ofstream fdata(str.c_str(),std::ofstream::out);
    for(int t=0;t<T_MAX;t++)
        fdata << t << " " << sqrt(ecarts[t]/N_MAX) << endl;
    fdata.close();
    delete[] system;
    delete[] ecarts;

return 0;
}


bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f){
    int hn  = array[modulo(x-s,LX)],
    hx  = array[modulo(x  ,LX)],
    hxp = array[modulo(x+s,LX)],
    hnp = array[modulo(x+2*s,LX)];
    if(hx+1 >= LY  or hxp-1 <= -LY) return false;
    int ene = abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
    double mag = Champ*(abs(hx+1)-abs(hx)+abs(hxp-1)-abs(hxp));
    double sos = J*ene + f*s+mag/2;

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}
