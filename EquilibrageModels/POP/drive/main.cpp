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
const double Tmin = 0.5, Tmax = 10;
const int Tn = 1;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 6, LX = 2<<bitLX;
const int LY = 25;
int N; // nombre total de particules
const int T_EQ = 1e4, T_MAX = 0e6;
const double mu = -0;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosGlau(int array[], int x, int ajout,double kbeta,double Champ);
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

    int Nb_T = Tn/nbprocs;
    double* valeurs = new double[Nb_T];
    double* energies = new double[Nb_T];
    double valTot[Tn] ;
    double eneTot[Tn] ;

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        cout << k << endl;
        double T =  normspace(k,Tmin,Tmax,Tn);
        Beta = 1./(ttc*T);
        /******** Initialisation du systeme *****/
        N = round(LY*LX/2); 
        std::vector<int> syspos(N);
        int system[LX];
        for(int x=0;x<LX;x++)
            system[x]  = 0;
        for(int n=0;n<N;n++){
            syspos[n] = rand_lx(generator);
            system[syspos[n]]++;
        }
        int ene = 0;
        for(int x=0;x<LX;x++)
            ene += abs(system[x]-system[modulo(x+1,LX)]);

        // Équilibrage
        int eneq[T_EQ];
        int ajout,tirage,pos;
        int hx,hxm,hxp,hx2;
        for(long int t = 0; t<LX*T_EQ ; t++){
            ajout =2*r01int(generator)-1;
            if(ajout > 0){ //création de particule
                pos = rand_lx(generator);
                if(EsosGlau(system,pos,ajout,Beta,mu)){
                    hx  = system[pos];
                    hxp = system[modulo(pos+1,LX)];
                    hxm = system[modulo(pos-1,LX)];
                    hx2 = hx+ajout;
                    ene +=( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );
                    syspos.push_back(pos);
                    system[pos]+=ajout;
                    N+=ajout;
                }
            }
            else{ //destruction de particule
                if(N != 0){
                    tirage = std::uniform_int_distribution<int>{0,N-1}(generator);
                    pos = syspos[tirage];
                    if(EsosGlau(system,pos,ajout,Beta,mu)){
                        hx  = system[pos];
                        hxp = system[modulo(pos+1,LX)];
                        hxm = system[modulo(pos-1,LX)];
                        hx2 = hx+ajout;
                        ene +=( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );
                        syspos.erase(syspos.begin()+tirage);
                        system[pos]+=ajout;
                        N+=ajout;
                    }
                }
            }
            if(t%LX==0){
                eneq[t/LX] = ene;
//                eneq[static_cast<int>(t/LX)] = N;
//                cout << N << endl;
            }
        }

        str = prefix+"/eq"+to_string(T);
        ofstream feq(str.c_str(),std::ofstream::out);
        for(int t=0;t<T_EQ;t++)
            feq << t << " " << eneq[t] << endl;
        feq.close();


        //Calcul des valeurs initiales        
       
//        valeurs[cmpt] = 0;
//        energies[cmpt] = 0;
//        for(long int t = 0; t<LX*T_MAX ; t++){
//            tirage = rand_pos(generator);
//            ajout = 2*r01int(generator)-1;
//            if(EsosGlau(system,tirage,ajout,Beta,0)){
//                hx = system[tirage];
//                hn = system[modulo(tirage-1,LX)];
//                hp= system[modulo(tirage+1,LX)];
//                ene -= abs(hx-hn) + abs(hx-hp);
//                ene += abs(hx+ajout-hn) + abs(hx+ajout-hp);
//                //Détruire particule
//                if(ajout < 0 ){
//                    system[syspos[tirage]]--;
//                    syspos.erase(syspos.begin()+tirage);
//                    N--;
//                }
//                //Création particule
//                if(ajout > 0){
//                    syspos.push_back(tirage);
//                    N++;
//                    system[syspos[tirage]]++;
//                }
//                std::discrete_distribution<> rand_pos(0,N);
//            }
//            valeurs[cmpt] += N;
//            energies[cmpt] += ene;
//        }
//        valeurs[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
//        energies[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_T,MPI_DOUBLE,valTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_T,MPI_DOUBLE,eneTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);

//    /*** Écriture dans fichier ***/
//    if(rang == 0){
//        str = prefix+"/dataSOS";
//        ofstream fmag(str.c_str(),std::ofstream::app);
//        for(int i=0;i<Tn;i++)
//            fmag << normspace(i,Tmin,Tmax,Tn)<< "\t" << valTot[i] <<  "\t\t" <<  J*eneTot[i]<<  "\n";
//        fmag.close();
//    }
    MPI_Finalize();
//    delete[] valeurs;
//    delete[] energies;
//    delete[] valTot;
//    delete[] eneTot;
return 0;
}
/********** Fin main ***********/

bool EsosGlau(int array[], int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<0 or hx2>2*LY){
        return false;
    }
    double champ = Champ*(hx2-hx)/2;
    double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}
