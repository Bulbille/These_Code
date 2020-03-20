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
const double Tmin = 0.7, Tmax = 0.95;
const int Tn = 4;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 25;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const long int T_EQ = 1e4, T_MAX = 1e6,T_COR = 1e4;
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

    /***** Parallélisation sur les champs magnétiques ****/
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        cout << k << endl;
        double T =  normspace(k,Tmin,Tmax,Tn);
        Beta = 1./(ttc*T);
        /******** Initialisation du systeme *****/
        //Mise en place du vecteur des positions
        int* system = new int[LX]; 
        //
        //Début à T=0
        double ene = 0, phi = 0;
        for(int x=0;x<LX;x++)
            system[x]  = 0;
//        for(int n=0;n<N;n++)
//            system[rand_lx(generator)]++;
//
        // Équilibrage
        int ajout,tirage;
        int hx,hn,hp;
        double* mageq = new double[T_EQ];
        double* eneq = new double[T_EQ];
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,0)){
                hx = system[tirage];
                hn = system[modulo(tirage-1,LX)];
                hp= system[modulo(tirage+1,LX)];
                ene -= abs(hx-hn) + abs(hx-hp);
                ene += abs(hx+ajout-hn) + abs(hx+ajout-hp);
                phi += ajout;
                system[tirage]+=ajout;
            }
            /***** Mesure du système *******/
            if(t%LX == 0){
                mageq[t/LX]     = phi;
                eneq[t/LX]  = ene;
            }
        }
        str = prefix+"/eq"+to_string(T);
        ofstream feq(str.c_str(),std::ofstream::out);
        for(int t=0;t<T_EQ;t++)
            feq << t << " " << 1.*mageq[t]/LX << " " << 1.*eneq[t]/LX << endl;
        feq.close();
        double* magcor = new double[T_MAX];
        double* enecor = new double[T_MAX];
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
                system[tirage]+=ajout;
            }
            if(t%LX == 0){
                magcor[t/LX]  = phi;
                enecor[t/LX]  = ene;
            }
        }
        /** Fonction d'autocorrélation **/
        phi /=static_cast<double>(LX*T_MAX*LX);
        double* magcorfun = new double[T_MAX];
        double* enecorfun = new double[T_MAX];
        for(int t=0;t<T_COR;t++){
            double ma=0,mb=0,mc=0;
            double ea=0,eb=0,ec=0;
            for(int tp=0;tp<T_MAX-t;tp++){
                ma += magcor[tp]*magcor[tp+t] ;
                mb += magcor[tp];
                mc += magcor[tp+t];
                ea += enecor[tp]*enecor[tp+t] ;
                eb += enecor[tp];
                ec += enecor[tp+t];
            }
            magcorfun[t] = ma / static_cast<double>(T_MAX-t) - mb * mc / static_cast<double>(pow(T_MAX-t,2)) ;
            enecorfun[t] = ea / static_cast<double>(T_MAX-t) - eb * ec / static_cast<double>(pow(T_MAX-t,2)) ;
        }
        /*** Intégration pour trouver le temps d'autoccorélation ***/
        double magtau = 0, enetau = 0;
        for(int t=0;t<T_COR;t++){
            if(t == 0 or t == T_MAX-1){
                magtau += magcorfun[t]/magcorfun[0];
                enetau += enecorfun[t]/enecorfun[0];
            }
            else if(t % 2 == 0){
                magtau += 2*magcorfun[t]/magcorfun[0];
                enetau += 2*enecorfun[t]/enecorfun[0];
            }
            else{
                magtau += 4*magcorfun[t]/magcorfun[0];
                enetau += 4*enecorfun[t]/enecorfun[0];
            }
        }
        magtau /= 3;
        enetau /= 3;
        str = prefix+"/cor"+to_string(T);
        ofstream fcor(str.c_str(),std::ofstream::out);
        for(int t=0;t<T_COR;t++)
            fcor << t << " " << magcorfun[t]/magcorfun[0] << " " << enecorfun[t]/enecorfun[0] << endl;
        fcor.close();
        cout << magtau << " " << enetau << endl;
        str = prefix+"/thermo";
        ofstream fthermo(str.c_str(),std::ofstream::app);
        fthermo << "T " << magtau << " " << enetau << endl;
        fthermo.close();
    }
    MPI_Finalize();

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
