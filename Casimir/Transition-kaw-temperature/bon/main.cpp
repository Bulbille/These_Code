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
const double Tmin = 0.4, Tmax = 0.5;
const int Tn = 1;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
int lambdaN = 2;
const int LMax = 300, LY = 30;
int rangee = LY/2;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e6,T_COR = 1e4;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f, double couplage, int k0);
double normspace(int step,double min,double max, double n){
    return (n == 1) ? max : min+(max-min)/(n-1)*1.*step;
}
int sign(double a){
    return (a == 0) ? 0 : ( (a > 0 ) ? 1 : -1 );
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
    while(lambdaN%nbprocs != 0) lambdaN++;

    int Nb_L = lambdaN/nbprocs;
    /***** Boucle sur les températures ****/
    for(int kt = 0; kt < Tn ; kt++){
        Beta = 1./(ttc*normspace(kt,Tmin,Tmax,Tn));
        cout << kt << endl;

        /*** Parallélisation sur le découplage lambda ***/
        double* FrTot = new double[lambdaN];
        double* FrEneL = new double[Nb_L];
        int cmpt = 0;
        for(int k = rang*Nb_L; k<(rang+1)*Nb_L ; k++){
            double lambda =  normspace(k,0,1,lambdaN);
            /******** Initialisation du systeme *****/
            int system [LX] = {}; 
            for(int n=0;n<N;n++)
                system[rand_lx(generator)]++;
            int sens,tirage;
            int hx,hn,hxp,hnp,y;
            /// Équilibrage
            for(long int t = 0; t<LX*T_EQ ; t++){
                tirage = rand_lx(generator);
                sens   = 2*r01int(generator)-1;
                if(EsosKaw(system,tirage,sens,Beta,0,0,lambda,rangee)){
                    system[tirage]++;
                    system[modulo(tirage+sens,LX)]--;
                }
            }
            /// Simulation
            double ene0 = 0,diffEne = 0;
            for(int x=0;x<LX;x++){
                y = system[x];
//                ene0 += abs(y-system[modulo(x,LX)]);
                diffEne += sign(rangee-1-y)*sign(rangee+1-y) - sign(rangee-y) * ( sign(rangee-1-y)+sign(rangee+1-y) );
            }
            FrEneL[cmpt] = diffEne;

            for(long int t = 0; t<LX*T_MAX ; t++){
                tirage = rand_lx(generator);
                sens   = 2*r01int(generator)-1;
                if(EsosKaw(system,tirage,sens,Beta,0,0,lambda,rangee)){
                    hn  = system[modulo(tirage-sens,LX)];
                    hx  = system[modulo(tirage  ,LX)];
                    hxp = system[modulo(tirage+sens,LX)];
                    hnp = system[modulo(tirage+2*sens,LX)];
//                    ene += abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
//                    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
                    int tab[4] = {hn,hx,hxp,hnp};
                    for(int i=0;i<4;i++){
                        y = tab[i];
                        diffEne += sign(rangee-1-y)*sign(rangee+1-y) - sign(rangee-y) * ( sign(rangee-1-y)+sign(rangee+1-y) );
                    }
                    system[tirage]++;
                    system[modulo(tirage+sens,LX)]--;
                }
                FrEneL[cmpt] += diffEne;
            }
            FrEneL[cmpt] *= -J/2./static_cast<double>(LX*LX*T_MAX);
            cmpt++;
        }
        MPI_Gather(FrEneL,Nb_L,MPI_DOUBLE,FrTot,Nb_L,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(rang == 0 ){
            str = prefix+"/t"+to_string(Beta)+".ran"+to_string(rangee);
            ofstream flambda(str.c_str(),std::ofstream::out);
            for(int l=0;l<lambdaN;l++)
                flambda << normspace(l,0,1,lambdaN) << " " << FrTot[l] << endl;
            flambda.close();
        }
    }
    MPI_Finalize();

return 0;
}
/********** Fin main ***********/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f, double couplage, int k0){
    int hn  = array[modulo(x-s,LX)],
        hx  = array[modulo(x  ,LX)],
        hxp = array[modulo(x+s,LX)],
        hnp = array[modulo(x+2*s,LX)];
    if(hx+1 >= LY  or hxp-1 <= -LY) return false;
    int ene = abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
    double mag = Champ*(abs(hx+1)-abs(hx)+abs(hxp-1)-abs(hxp));
    double decouplage = 0;
    int tab[4] = {hn,hx,hxp,hnp};
    for(int i=0;i<4;i++){
        int y = tab[i];
        decouplage += sign(k0-1-y)*sign(k0+1-y) - sign(k0-y) * ( sign(k0-1-y)+sign(k0+1-y) );
    }

    double sos = J*ene + f*s+mag/2 - couplage*J/2*decouplage;
    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}
