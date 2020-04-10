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
const double Fmin = 3.7, Fmax = 4.3;
const int Fn = 20;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
const double mu = 0.1;
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 25;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e7;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f);
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

    string str;
    int Nb_F = Fn/nbprocs;

    double* energies = new double[Nb_F];
    double* eneTot = new double[Fn];


    int cmpt = 0;

    /***** Parallélisation sur les champs magnétiques ****/
    for(int k = rang*Nb_F; k<(rang+1)*Nb_F ; k++){
        cout << k << endl;
        double drive =  normspace(k,Fmin,Fmax,Fn);
        /******** Initialisation du systeme *****/
        //Mise en place du vecteur des positions
        int* system = new int[LX]; 
        //
        //Début à T=0
        for(int x=0;x<LX;x++)
            system[x]  = 0;
        // Équilibrage
        int sens,tirage;
        int hx,hn,hxp,hnp;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,mu,drive)){
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
        }
        energies[cmpt] = 0;
        double ene = 0;
        for(int x=0;x<LX;x++){
            ene += abs(system[x]-system[modulo(x+1,LX)]);
        }
        // Simulation

        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,mu,drive)){
                hn  = system[modulo(tirage-sens,LX)];
                hx  = system[modulo(tirage  ,LX)];
                hxp = system[modulo(tirage+sens,LX)];
                hnp = system[modulo(tirage+2*sens,LX)];
                ene += abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
                ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
            energies[cmpt] += ene;
        }
        energies[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(energies ,Nb_F,MPI_DOUBLE,eneTot,Nb_F,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(rang == 0){
//        str = prefix+"/X"+to_string(LX)+"_T"+to_string(ttc)+"_mu"+to_string(mu)+"_F"+to_string(Fmax);
        str = prefix+"/X"+to_string(LX)+"_T"+to_string(ttc)+"_mu"+to_string(mu);
        ofstream fmag(str.c_str(),std::ofstream::app);
        for(int i=0;i<Fn;i++)
            fmag << normspace(i,Fmin,Fmax,Fn) << " " << J*eneTot[i] << endl;
        fmag.close();
    }

    MPI_Finalize();
    delete[] energies;
    delete[] eneTot;

return 0;
}
/********** Fin main ***********/
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
