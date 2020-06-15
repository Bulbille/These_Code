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
const double Tmin = 0.01, Tmax = 6;
const int Tn = 100;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 10;
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e7;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ);
double normspace(int step,double min,double max, double n){
    return (n == 1) ? max : min+(max-min)/(n-1)*1.*step;
}
double logspace(int step,double min,double max, double n){
    return (n == 1) ? pow(10,max) : pow(10,min+(max-min)/(n-1)*1.*step);
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

    int Nb_T = Tn/nbprocs;
    double* valeurs = new double[Nb_T];
    double* valTot = new double[Tn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        valeurs[cmpt] = 0;
        double B =  normspace(k,Tmin,Tmax,Tn);
        cout << k << " " << B << endl;
        Beta = 1;
        /******** Initialisation du systeme *****/
        int phi = 0;
        //Mise en place du vecteur des positions
        int* system = new int[LX]; 
        for(int x=0;x<LX;x++)
            system[x]  = LY;

        //Déclaration du PRNG 

        // Équilibrage
        int sens,tirage;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,B)){
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
        }
        //Calcul des valeurs initiales        
        phi = 0;
        for(int x=0;x<LX;x++){
            phi += abs(system[x]-LY);
        }  
       
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,B)){
                phi -= abs(system[tirage]-LY)+abs(system[modulo(tirage+sens,LX)]-LY);
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
                phi += abs(system[tirage]-LY)+abs(system[modulo(tirage+sens,LX)]-LY);
            }
            /***** Mesure du système *******/
            valeurs[cmpt] += phi;
        }
        valeurs[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_T,MPI_DOUBLE,valTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/muSOS"+to_string(LY)+'H'+to_string(Tmax);
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Tn;i++)
            fmag << normspace(i,Tmin,Tmax,Tn)<< "\t" << valTot[i] <<  "\n";
        fmag.close();
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] valTot;
return 0;
}
/********** Fin main ***********/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ){
    int hn  = array[modulo(x-s,LX)],
    hx  = array[modulo(x  ,LX)],
    hxp = array[modulo(x+s,LX)],
    hnp = array[modulo(x+2*s,LX)];
    if(hx+1 > 2*LY  or hxp-1 < 0) return false;
    int ene = abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
    double mag = -Champ*(abs(hx+1-LY)-abs(hx-LY)+abs(hxp-1-LY)-abs(hxp-LY));
    double sos = J*ene + mag;

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}
