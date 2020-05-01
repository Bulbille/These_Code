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
const double H = 0.1;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
int lambdaN = 8;
int rangee = 1;
const int Lmin = 10, Lmax = 20;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e5;
// Taux de diffusion des particules A et B

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
bool EsosGlau(int taille, int* array, int x, int ajout,double kbeta,double Champ, double couplage, int k0);

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
    double* IntFre = new double[Lmax-Lmin];

    /***** Boucle sur les températures ****/
    for(int LY = Lmin ; LY< Lmax ; LY++){

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
            int ajout,tirage;
            int hx,hn,hxp,y;
            /// Équilibrage
            for(long int t = 0; t<LX*T_EQ ; t++){
                tirage = rand_lx(generator);
                ajout   = 2*r01int(generator)-1;
                if(EsosGlau(LY,system,tirage,ajout,Beta,H,lambda,rangee)){
                    system[tirage]+=ajout;
                }
            }
            /// Simulation
            double diffEne = 0;
            for(int x=0;x<LX;x++){
                y = system[x];
//                ene0 += abs(y-system[modulo(x,LX)]);
                diffEne += sign(rangee-1-y)*sign(rangee+1-y) - sign(rangee-y) * ( sign(rangee-1-y)+sign(rangee+1-y) );
            }
            FrEneL[cmpt] = diffEne;

            for(long int t = 0; t<LX*T_MAX ; t++){
                tirage = rand_lx(generator);
                ajout   = 2*r01int(generator)-1;
                if(EsosGlau(LY,system,tirage,ajout,Beta,H,lambda,rangee)){
                    hn  = system[modulo(tirage-1,LX)];
                    hx  = system[tirage];
                    hxp = system[modulo(tirage+1,LX)];
                    diffEne = 0;
                    int tab[3] = {hx,hn,hxp};
                    for(int i=0;i<1;i++){
                        y = tab[i];
                        diffEne += sign(rangee-1-y)*sign(rangee+1-y) - sign(rangee-y) * ( sign(rangee-1-y)+sign(rangee+1-y) );
                    }
                    cout << diffEne << " " ;
                    system[tirage]+=ajout;
                }
                FrEneL[cmpt] += diffEne;
            }
            FrEneL[cmpt] *= -J/2./static_cast<double>(LX*LX*T_MAX);
            cmpt++;
        }
        MPI_Gather(FrEneL,Nb_L,MPI_DOUBLE,FrTot,Nb_L,MPI_DOUBLE,0,MPI_COMM_WORLD);

        if(rang == 0){
            IntFre[LY-Lmin] = 0;
            for(int i=0;i<lambdaN;i++){
                if(i == 0 or i== lambdaN-1)
                    IntFre[LY-Lmin] += FrTot[i];
                else if (i%2 == 0 )
                    IntFre[LY-Lmin] += 2*FrTot[i];
                else
                    IntFre[LY-Lmin] += 4*FrTot[i];
            }
            cout << LY << " "<< IntFre[LY-Lmin]  << endl;
            IntFre[LY-Lmin] /= (lambdaN-1)/3;
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }
    if(rang == 0 ){
        str = prefix+"/ran"+to_string(rangee);
        ofstream flambda(str.c_str(),std::ofstream::out);
        for(int LY=Lmin;LY<Lmax;LY++)
            flambda << LY << " " << IntFre[LY-Lmin] << endl;
        flambda.close();
    }
    MPI_Finalize();

return 0;
}
/********** Fin main ***********/

bool EsosGlau(int taille, int* array, int x, int ajout,double kbeta,double Champ, double couplage, int k0){
    int hx  = array[x],
    hxp = array[modulo(x+1,LX)],
    hxm = array[modulo(x-1,LX)],
    hx2 = hx+ajout;
    if(hx2<=-taille or hx2>=taille) return false;
    double champ = Champ*(abs(hx2)-abs(hx))/2;
    double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );

    double decouplage = 0;
    int tab[3] = {hxp,hx,hxm};
    for(int i=0;i<3;i++){
        int y = tab[i];
        decouplage += sign(k0-1-y)*sign(k0+1-y) - sign(k0-y) * ( sign(k0-1-y)+sign(k0+1-y) );
    }

    sos -= couplage*J/2*decouplage;

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

