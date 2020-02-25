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
#include <Eigen/Eigenvalues> 
#include <Eigen/Dense>
using namespace Eigen;

//Modules diagonalisation de matrices
using namespace std;
const double T_C    = 2./log(1.+sqrt(2)), J = 1;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX;
const int Lmat = 200, LY = 400;
int N; // nombre total de particules
//LX doit être une puissance de 2, 2<<0 = 2, 2<<n = 2^(n+1)
const int  Mun = 3;
const long int T_EQ = 1e6, T_MAX = 1e8, T_MAXCOR = 1e6 ,  T_COR = 1e4;
// Taux de diffusion des particules A et B
double mumax =0; 
double mumin=-2; 

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
double getEquilibrium(int l,double chim);
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ);
void FnCorrel(double* cor, double* val, int tmax, int tcor);
double IntCor(double* cor, int tmax);
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
    ostringstream streamT,streamM;
    string str;
    streamM << fixed << setprecision(4) << mumax+mumin;
    streamT << fixed << setprecision(2) << ttc;

    int Nb_mu = Mun/nbprocs;
    double* valeurs = new double[Nb_mu];
    double* Tcor = new double[Nb_mu];
    double* Lcor = new double[Nb_mu];
    double* energies = new double[Nb_mu];
    double* ecartTyp = new double[Nb_mu];
    double* valTot = new double[Mun];
    double* eneTot = new double[Mun];
    double* ecartTot = new double[Mun];
    double* TcorTot = new double[Mun];
    double* LcorTot = new double[Mun];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_mu; k<(rang+1)*Nb_mu ; k++){
        valeurs[cmpt] = 0;
        energies[cmpt] = 0;
        ecartTyp[cmpt] = 0;
        double mu =  logspace(k,mumin,mumax,Mun); //Le -1 est pour arriver à mumax lorsque k = Nb_mu au rang supérieur
        cout << k << " " << mu << endl;
        /******** Initialisation du systeme *****/
        double mean = getEquilibrium(Lmat,mu);
        N = round(mean*LX);
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

            if(EsosGlau(system,tirage,ajout,Beta,mu))
                system[tirage]+=ajout;
        }
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
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,mu)){
                hx = system[tirage];
                hn = system[modulo(tirage-1,LX)];
                hp= system[modulo(tirage+1,LX)];
                ene -= abs(hx-hn) + abs(hx-hp);
                ene += abs(hx+ajout-hn) + abs(hx+ajout-hp);
                phi += ajout;
                phi2 += 2*hx*ajout+1;
                system[tirage]+=ajout;
            }
            if(t % LX == 0 )
                corene[t/LX] = ene;
            /***** Mesure du système *******/
            valeurs[cmpt] += phi;
            energies[cmpt] += ene;
            ecartTyp[cmpt] += phi2;
        }
        double* funcor = new double[T_COR];
        FnCorrel(funcor,corene,T_MAXCOR,T_COR);
        double temps = IntCor(funcor,T_COR);
        delete[] funcor;
        /**** CALCUL muISTOGRAMME ****/
        double* histo = new double[2*LY];
        for(int y=0;y<2*LY;y++) histo[y] = 0;

        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,mu))
                system[tirage]+=ajout;
            if(t % LX*static_cast<int>(temps) == 0 ){
                for(int x=0;x<LX;x++){
                    histo[ LY+system[x] ]++;
                }
            }
        }
        // écriture histogramme dans fichier
        str = prefix+"/histo"+to_string(LX)+"Ttc"+streamT.str()+"_mu"+to_string(mu);
        ofstream fhisto(str.c_str(),std::ofstream::out);
        for(int y=0;y<2*LY;y++)
            fhisto << y << " " << histo[y] << " " << endl;
        fhisto.close();
        str = prefix+"/snap"+to_string(LX)+"Ttc"+streamT.str()+"_mu"+to_string(mu);
        ofstream fsnap(str.c_str(),std::ofstream::out);
        for(int x=0;x<LX;x++)
            fsnap << x << " " << system[x] << endl;
        fsnap.close();
        delete[] histo;
        delete[] system;

        Tcor[cmpt] = temps;
        valeurs[cmpt] /= static_cast<double>(LX*LX*T_MAXCOR) ;
        energies[cmpt] /= static_cast<double>(LX*LX*T_MAXCOR) ;
        ecartTyp[cmpt] /= static_cast<double>(LX*LX*T_MAXCOR) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_mu,MPI_DOUBLE,valTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_mu,MPI_DOUBLE,eneTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_mu,MPI_DOUBLE,ecartTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(Tcor,Nb_mu,MPI_DOUBLE,TcorTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(Lcor,Nb_mu,MPI_DOUBLE,LcorTot,Nb_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/X"+to_string(LX)+"Ttc"+streamT.str()+"_mu"+streamM.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Mun;i++)
            fmag << logspace(i,mumin,mumax,Mun)<< "\t" << valTot[i] <<  "\t\t" <<  ecartTot[i] << "\t\t" << J*eneTot[i]<< " " << TcorTot[i] << "\n";
        fmag.close();
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] energies;
    delete[] ecartTyp;
    delete[] Tcor;
    delete[] Lcor;
    delete[] valTot;
    delete[] eneTot;
    delete[] ecartTot;
    delete[] TcorTot;
    delete[] LcorTot;
return 0;
}
/********** Fin main ***********/
double getEquilibrium(int l,double chim){
    EigenSolver<MatrixXd> es;
    MatrixXd transfer(l,l);
    MatrixXd diag(l,l);
    for(int i=0;i<l;i++){
        diag(i,i) = i;
        for(int j=0;j<l;j++){
            transfer(i,j) = exp(-Beta*( J*abs(i-j)- chim*(abs(i-l)+abs(j-l))/2)) ;
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


bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<=-LY or hx2>=LY) return false;
    double champ = Champ*(abs(hx2)-abs(hx))/2;
    double sos    = J*( pow(hx2-hxm,2) + pow(hx2-hxp,2) - (pow(hx - hxm,2)  +  pow(hx - hxp,2) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
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
        cout << t << " " << cor[t] << endl;

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

