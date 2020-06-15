/******************************************************************/
// Plan SemiInfini, dynamique de Kawasaki et de Glauber
// Modèle Chipping, hamiltonien SOS
// Parallélisé sur les cisaillements
/******************************************************************/



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
double  ttc = 3 ,Beta = 1/ttc;
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 800;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e7;
// Taux de diffusion des particules A et B
double mumax = 0.0;
const int  Dn = 100;
double dmin = 0, dmax = 6, dJ =0.1;

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f); 

double space(int step){
    if(Dn == 1 ) return dmax;
    if(step < Dn/5 ) return (2*J-dJ)*1.*step/(Dn/5);
    if(step< Dn*2/5) return (2*J-dJ)+2*dJ*1.*(step-Dn/5)/(Dn/5);
    if(step< Dn*3/5) return (2*J+dJ)+(2*J-2*dJ)*1.*(step-Dn*2/5)/(Dn/5);
    if(step< Dn*4/5) return (4*J-dJ)+2*dJ*1.*(step-Dn*3/5)/(Dn/5);
    if(step< Dn) return (4*J+dJ)+(dmax-4*J-dJ)*1.*(step-Dn*4/5)/(Dn/5);
}


/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
    ostringstream streamF,streamT,streamM;
    string str;
    streamT << fixed << setprecision(2) << ttc; 
    streamM << fixed << setprecision(2) << mumax; 
    streamF << fixed << setprecision(2) << dmax;

    int Nb_f = Dn/nbprocs;
    double* energies = new double[Nb_f];
    double* eneTot = new double[Dn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_f; k<(rang+1)*Nb_f ; k++){
        Beta = 1;
        energies[cmpt] = 0;
        double drive =  space(k); //Le -1 est pour arriver à mumax lorsque k = Nb_f au rang supérieur
        /******** Initialisation du systeme *****/
        int* system = new int[LX]; 
        double mean = 4.51;
        int N = mean*LX;
        for(int x = 0;x<LX;x++)
            system[x] = 0;
        for(int n=0;n<N;n++)
            system[rand_lx(generator)] ++;

        // Équilibrage
        int sens,tirage;
        int hx,hn,hxp,hnp;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,mumax,drive)){
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
        }
        //Calcul des valeurs initiales        
        double ene = 0;
        for(int x=0;x<LX;x++){
//            ene += pow(system[x]-system[modulo(x+1,LX)],2);
            ene += abs(system[x]-system[modulo(x+1,LX)]);
        }  
        /**** CALCUL muISTOGRAMME ****/
        int* histo = new int[LY];
        for(int y=0;y<LY;y++) histo[y] = 0;
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,mumax,drive)){
                hx  = system[tirage];
                hn  = system[modulo(tirage-sens,LX)];
                hxp  = system[modulo(tirage+sens,LX)];
                hnp = system[modulo(tirage+2*sens,LX)];
                ene += abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
                ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
//                ene -= pow(hx-hn,2)     + pow(hx-hxp,2) + pow(hxp-hnp,2);
//                ene += pow(hx+1-hn,2)   + pow(hx-hxp+2,2)+pow(hxp-1-hnp,2);
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
            energies[cmpt] += ene;
        }

        delete[] system;

        cout << k << " " << drive <<  " " << energies[cmpt] << " " << energies[cmpt]/static_cast<double>(LX*LX*T_MAX) <<endl;
        energies[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(energies,Nb_f,MPI_DOUBLE,eneTot,Nb_f,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/sos"+to_string(LX)+"_F"+streamF.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Dn;i++)
            fmag << space(i) <<  "\t\t" << J*eneTot[i]<< "\n";
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
    if(hx+1 >= LY  or hxp-1 <0) return false;
    int ene = abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
//    int ene = pow(hx+1-hn,2) + pow(hx-hxp+2,2)+pow(hxp-1-hnp,2);
//    ene -= pow(hx-hn,2) + pow(hx-hxp,2)+pow(hxp-hnp,2);
    double mag = Champ*(abs(hx+1)-abs(hx)+abs(hxp-1)-abs(hxp));
    double sos = J*ene + f*s+mag;

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

