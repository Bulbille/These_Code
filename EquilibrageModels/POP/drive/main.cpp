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
const double Tmin = 0, Tmax = 1.5*T_C;
const int Tn = 10;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 400;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e6, T_MAX = 1e6;
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

    int Nb_T = Tn/nbprocs;
    double* valeurs = new double[Nb_T];
    double* energies = new double[Nb_T];
    double* ecartTyp = new double[Nb_T];
    double* valTot = new double[Tn];
    double* eneTot = new double[Tn];
    double* ecartTot = new double[Tn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        cout << k << endl;
        valeurs[cmpt] = 0;
        energies[cmpt] = 0;
        ecartTyp[cmpt] = 0;
        double T =  normspace(k,Tmin,Tmax,Tn);
        Beta = 1./(ttc*T);
        /******** Initialisation du systeme *****/
        int phi = 0;
        //Mise en place du vecteur des positions
        int* system = new int[LX]; 
        for(int x=0;x<LX;x++)
            system[x]  = 0;
        for(int n=0;n<N;n++)
            system[rand_lx(generator)]++;

        //Déclaration du PRNG 

        // Équilibrage
        int ajout,tirage;
        int hx,hn,hp;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;
            if(EsosGlau(system,tirage,ajout,Beta,0))
                system[tirage]+=ajout;
        }
        //Calcul des valeurs initiales        
        double phi2 = 0,ene = 0;
        phi = 0;
        for(int x=0;x<LX;x++){
            phi += system[x];
            phi2 += pow(system[x],2);
            ene += abs(system[x]-system[modulo(x+1,LX)]);
        }  
       
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,0)){
                hx = system[tirage];
                hn = system[modulo(tirage-1,LX)];
                hp= system[modulo(tirage+1,LX)];
                ene -= abs(hx-hn) + abs(hx-hp);
                ene += abs(hx+ajout-hn) + abs(hx+ajout-hp);
                phi += ajout;
                phi2 += 2*hx*ajout+1;
                system[tirage]+=ajout;
            }
            /***** Mesure du système *******/
            valeurs[cmpt] += phi;
            energies[cmpt] += ene;
            ecartTyp[cmpt] += phi2;
        }
        valeurs[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        energies[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        ecartTyp[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_T,MPI_DOUBLE,valTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_T,MPI_DOUBLE,eneTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_T,MPI_DOUBLE,ecartTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/dataSOS";
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Tn;i++)
            fmag << normspace(i,Tmin,Tmax,Tn)<< "\t" << valTot[i] <<  "\t\t" <<  ecartTot[i] << "\t\t" << J*eneTot[i]<<  "\n";
        fmag.close();
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] energies;
    delete[] ecartTyp;
    delete[] valTot;
    delete[] eneTot;
    delete[] ecartTot;
return 0;
}
/********** Fin main ***********/

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
