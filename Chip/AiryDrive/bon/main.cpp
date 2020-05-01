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
double  ttc = 4 ,Beta = 1/ttc;
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 8, LX = 2<<bitLX;
const int LY = 800;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e5, T_MAX = 1e7;
// Taux de diffusion des particules A et B
double mumax = 0.01;
const int  Dn = 1;
double dmin = 0;
double dmax = 0;

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f); 
double logspace(int step,double min,double max, double n){
    return (n == 1) ? pow(min,10) : pow(10,(min-max)/(n-1)*step);
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
    ostringstream streamF,streamT,streamM;
    string str;
    streamT << fixed << setprecision(2) << ttc; 
    streamM << fixed << setprecision(2) << mumax; 
    streamF << fixed << setprecision(2) << dmax;

    int Nb_f = Dn/nbprocs;
    double* energies = new double[Nb_f];
    double* ecartTyp = new double[Nb_f];
    double* eneTot = new double[Dn];
    double* ecartTot = new double[Dn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_f; k<(rang+1)*Nb_f ; k++){
        energies[cmpt] = 0;
        ecartTyp[cmpt] = 0;
        double drive =  normspace(k,dmin,dmax,Dn); //Le -1 est pour arriver à mumax lorsque k = Nb_f au rang supérieur
        cout << k << " " << drive << endl;
        /******** Initialisation du systeme *****/
        int* system = new int[LX]; 
        for(int x=0;x<LX;x++)
            system[x]  = LY/2;

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
        double phi2 = 0,ene = 0;
        for(int x=0;x<LX;x++){
            phi2 += pow(system[x],2);
            ene += pow(system[x]-system[modulo(x+1,LX)],2);
        }  
        /**** CALCUL muISTOGRAMME ****/
        int* histo = new int[2*LY];
        for(int y=0;y<2*LY;y++) histo[y] = 0;
        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,mumax,drive)){
                hx  = system[tirage];
                hn  = system[modulo(tirage-sens,LX)];
                hxp  = system[modulo(tirage+sens,LX)];
                hnp = system[modulo(tirage+2*sens,LX)];
                ene -= pow(hx-hn,2)     + pow(hx-hxp,2) + pow(hxp-hnp,2);
                ene += pow(hx+1-hn,2)   + pow(hx-hxp+2,2)+pow(hxp-1-hnp,2);

                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
            if(t % LX == 0 ){
                for(int x=0;x<LX;x++){
                    histo[ LY+system[x] ]++;
                }
                energies[cmpt] += ene;
                ecartTyp[cmpt] += phi2;

            }
        }

        // écriture histogramme dans fichier
        str = prefix+"/histo"+to_string(LX)+"Ttc"+streamT.str()+ "_mu"+streamM.str()+"_f"+to_string(drive);
        ofstream fhisto(str.c_str(),std::ofstream::out);
        for(int y=0;y<2*LY;y++)
            fhisto << y << " " << histo[y] << " " << endl;
        fhisto.close();
        delete[] histo;
        delete[] system;

        energies[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        ecartTyp[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(energies,Nb_f,MPI_DOUBLE,eneTot,Nb_f,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_f,MPI_DOUBLE,ecartTot,Nb_f,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
//    if(rang == 0){
//        str = prefix+"/X"+to_string(LX)+"Ttc"+streamT.str()+"_mu"+streamM.str()+"_F"+streamF.str();
//        ofstream fmag(str.c_str(),std::ofstream::out);
//        for(int i=0;i<Dn;i++)
//            fmag << normspace(i,dmin,dmax,Dn) << "\t" <<   ecartTot[i] << "\t\t" << J*eneTot[i]<< "\n";
//        fmag.close();
//    }
    MPI_Finalize();
    delete[] energies;
    delete[] ecartTyp;
    delete[] eneTot;
    delete[] ecartTot;
return 0;
}
/********** Fin main ***********/

bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f){ 
    int hn  = array[modulo(x-s,LX)],
        hx  = array[modulo(x  ,LX)],
        hxp = array[modulo(x+s,LX)],
        hnp = array[modulo(x+2*s,LX)];
    if(hx+1 >= LY  or hxp-1 <= 0) return false;
//    int ene = abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
//    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
    int ene = pow(hx+1-hn,2) + pow(hx-hxp+2,2)+pow(hxp-1-hnp,2);
    ene -= pow(hx-hn,2) + pow(hx-hxp,2)+pow(hxp-hnp,2);
    double mag = Champ*(abs(hx+1)-abs(hx)+abs(hxp-1)-abs(hxp));
    double sos = J*ene + f*s+mag;

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

