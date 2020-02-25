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
#include <Eigen/Eigenvalues> 
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;
const double T_C    = 2./log(1.+sqrt(2)), J = 1;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int LY = 200;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const int  Mun = 10;
const long int T_EQ = 1e6, T_MAX = 1e7, T_MAXCOR = 1e6 ,  T_COR = 1e4;
// Taux de diffusion des particules A et B
double mumax =-1; 
double mumin =-3; 
double drive = 0;

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
double getEquilibrium(int l,double chim);
bool EsosKaw(int* array, int x, int s, double kbeta,double Champ,double f); 
void FnCorrel(double* cor, double* val, int tmax, int tcor);
double IntCor(double* cor, int tmax);
double logspace(int step,double min,double max, double n){
    return (n == 1) ? pow(min,10) : pow(10,(min-max)/(n-1)*step);
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
    streamM << fixed << setprecision(2) << mumax;
    streamF << fixed << setprecision(2) << drive;

    int Nb_mu = Mun/nbprocs;
    double* valeurs = new double[Nb_mu];
    double* Tcor = new double[Nb_mu];
    double* energies = new double[Nb_mu];
    double* ecartTyp = new double[Nb_mu];
    double* valTot = new double[Mun];
    double* eneTot = new double[Mun];
    double* ecartTot = new double[Mun];
    double* TcorTot = new double[Mun];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_mu; k<(rang+1)*Nb_mu ; k++){
        valeurs[cmpt] = 0;
        energies[cmpt] = 0;
        ecartTyp[cmpt] = 0;
        double mu =  logspace(k,mumin,mumax,Mun); //Le -1 est pour arriver à mumax lorsque k = Nb_mu au rang supérieur
        /******** Initialisation du systeme *****/
        double mean; // = getEquilibrium(LY,mu);
        mean = 0;
        N = round(mean*LX);
        int phi = 0;
        int* system = new int[LX]; 
        for(int x=0;x<LX;x++)
            system[x]  = 0;
        for(int n=0;n<N;n++)
            system[rand_lx(generator)]++;

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
        cout << Beta << " " << mu << " " << drive << endl;
        str = prefix+"/snap"+to_string(LX)+"_mu"+to_string(mu);
        ofstream fsnap(str.c_str(),std::ofstream::out);
    //****** CALCUL TEMPS DE CORRÉLATION *********//
        //Calcul des valeurs initiales        
        double phi2 = 0,ene = 0;
        phi = 0;
        for(int x=0;x<LX;x++){
            phi += system[x];
            phi2 += pow(system[x],2);
            ene += abs(system[x]-system[modulo(x+1,LX)]);
        }  
       
        /****** CALCUL TEMPS CORRÉLATION ***/
        double* corene = new double[T_MAXCOR];
        for(long int t = 0; t<LX*T_MAXCOR ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;

            if(EsosKaw(system,tirage,sens,Beta,mu,drive)){
                hx  = system[tirage];
                hn  = system[modulo(tirage-sens,LX)];
                hxp  = system[modulo(tirage+sens,LX)];
                hnp = system[modulo(tirage+2*sens,LX)];
                ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
                ene += abs(hx+1-hn) + abs(hx-hxp+2)+abs(hxp-1-hnp);
                phi2 += 2*(hxp-hx)+2;
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;

            }
            if(t % LX == 0 )
                corene[t/LX] = ene;
            /***** Mesure du système *******/
            energies[cmpt] += ene;
            ecartTyp[cmpt] += phi2;
        }
        double* funcor = new double[T_COR];
        FnCorrel(funcor,corene,T_MAXCOR,T_COR);
        double temps = IntCor(funcor,T_COR);
        temps = 100;
        delete[] funcor;
        /**** CALCUL muISTOGRAMME ****/
        int* histo = new int[2*LY];
        for(int y=0;y<2*LY;y++) histo[y] = 0;

        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,tirage,sens,Beta,mu,drive)){
                system[tirage]++;
                system[modulo(tirage+sens,LX)]--;
            }
            if(t % LX*static_cast<int>(temps) == 0 ){
                for(int x=0;x<LX;x++){
                    histo[ LY+system[x] ]++;
                }
            }
        }

        // écriture histogramme dans fichier
        str = prefix+"/histo"+to_string(LX)+"_mu"+to_string(mu);
        ofstream fhisto(str.c_str(),std::ofstream::out);
        for(int y=0;y<2*LY;y++)
            fhisto << y << " " << histo[y] << " " << endl;
        fhisto.close();
        for(int x=0;x<LX;x++)
            fsnap << x << " " << system[x] << endl;
        fsnap.close();
        delete[] histo;
        delete[] system;

        valeurs[cmpt] = mean;
        Tcor[cmpt] = temps;
        energies[cmpt] /= static_cast<double>(LX*LX*T_MAXCOR) ;
        ecartTyp[cmpt] /= static_cast<double>(LX*LX*T_MAXCOR) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_mu,MPI_DOUBLE,valTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_mu,MPI_DOUBLE,eneTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_mu,MPI_DOUBLE,ecartTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(Tcor,Nb_mu,MPI_DOUBLE,TcorTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/X"+to_string(LX)+"_mu"+streamM.str()+"_F"+streamF.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Mun;i++)
            fmag << logspace(i,mumin,mumax,Mun) << "\t" << valTot[i] <<  "\t\t" <<  ecartTot[i] << "\t\t" << J*eneTot[i]<< " " << TcorTot[i] << "\n";
        fmag.close();
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] energies;
    delete[] ecartTyp;
    delete[] Tcor;
    delete[] valTot;
    delete[] eneTot;
    delete[] ecartTot;
    delete[] TcorTot;
return 0;
}
/********** Fin main ***********/
// Diagonalisation de la matrice
double getEquilibrium(int l,double chim){
    EigenSolver<MatrixXd> es;
    MatrixXd transfer(l,l);
    MatrixXd diag(l,l);
    for(int i=0;i<l;i++){
        diag(i,i) = i;
        for(int j=0;j<l;j++){
            transfer(i,j) = exp(-Beta*( J*abs(i-j)+ chim*(i+j)/2)) ;
            if(i!=j)
                diag(i,j) =0 ;
        }
    }
    int imax = 0;
    double max = 0;
    es.compute(transfer,true);
    VectorXcd a = es.eigenvalues();
    for(int i=0;i<l;i++){
        if(a(i).real()>max){
            imax=i;
            max=a(i).real();
        }
    }
    VectorXcd Vmax = es.eigenvectors().col(imax);
    complex<double> res = Vmax.transpose()*diag*Vmax;
    return res.real();
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

// Calcul de la fonction de corrélation
void FnCorrel(double* cor, double* val, int tmax, int tcor){
    for(int t=0;t<tcor;t++){
        double a=0,b=0,c=0;
        for(int tp=0;tp<tmax-t;tp++){
            a += val[tp]*val[modulo(tp+t,tmax)] ;
            b += val[tp];
            c += val[modulo(tp+t,tmax)];
        }
        cor[t] = a / static_cast<double>(tmax-t) - b * c / static_cast<double>(pow(tmax-t,2)) ;
    }
}

double IntCor(double* cor, int tmax){
    /*** Intégration pour trouver le temps d'autoccorélation ***/
    double tau = 0;
    for(int t=0;t<tmax;t++){
        if(t == 0 or t == tmax-1)
            tau += cor[t];
        else if(t % 2 == 0)
            tau += 2*cor[t];
        else
            tau += 4*cor[t];
    }   
    return tau/3/cor[0];
}

