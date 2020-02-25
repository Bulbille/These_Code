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
//#include <gsl/gsl_fft_real.h>
#include <complex>
#define Imag dcomplex(0.0,1.0)
#include <gsl/gsl_fft_real.h>
#define TOT(z,i,l)  pow(pow(z[i],2)+pow(z[l-i],2),0.5)
typedef std::complex<double> dcomplex;

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
const int  Hn = 1;
const long int T_EQ = 1e6, T_MAX = 1e5, T_COR = 2e2;
double mumax = 8.5;
// Taux de diffusion des particules A et B
const double Da = 1, Db = 1;
//taux_dyn : pourcentage de particuls A dans le total
double taux_dyn =  0.51; 
const double drive =0; 

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
double getEquilibrium(int l,double chim);
bool EsosKaw(int* array, int x, int sens, double kbeta,double f);
void FnCorrel(double* cor, double* val, int tmax, int tcor);
double IntCor(double* cor, int tmax);
void SpCorrel(double*cor,double*val,int LX,double mbar);
double fft_maison(double* var,double q,int L){ 
    std::complex<double> sq = 0;
    for(double x=0;x<L;x++)
        sq += var[static_cast<int>(x)]*exp(2.*M_PI*Imag*q*x/static_cast<double>(L));
    return norm(sq)/L;
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
    streamT << fixed << setprecision(2) << taux_dyn;
    streamM << fixed << setprecision(2) << mumax;
    streamF << fixed << setprecision(2) << drive;

    int Nb_H = Hn/nbprocs;
    double* valeurs = new double[Nb_H];
    double* energies = new double[Nb_H];
    double* ecartTyp = new double[Nb_H];
    double* valTot = new double[Hn];
    double* eneTot = new double[Hn];
    double* ecartTot = new double[Hn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    for(int k = rang*Nb_H; k<(rang+1)*Nb_H ; k++){
        double f = (Hn==1) ? drive : static_cast<int>(k)*drive/(Hn-1); //Le -1 est pour arriver à Hmax lorsque k = Nb_H au rang supérieur
        /******** Initialisation du systeme *****/
        double mean = getEquilibrium(LY,mumax);
        N = round(mean*LX);
        int phi = 0;
        int* syspos = new int[N];
        for(int n=0;n<N;n++){
            if(phi < N*taux_dyn)
                syspos[n] = 2*rand_lx(generator);
            else
                syspos[n] = 2*rand_lx(generator)+1;  
            phi++;
        }
        //Mise en place du vecteur des positions
        int* system = new int[2*LX]; 
        double* hauteur = new double[LX]; 
        std::vector< double> weight(N);
        for(int x=0;x<2*LX;x++)
            system[x]  = 0;
        for(int n=0;n<N;n++){
            system[syspos[n]]++;
            weight[n] = Da*(syspos[n]%2  == 0)+Db*(syspos[n]%2 == 1);
        }
        //Déclaration du PRNG 
        std::discrete_distribution<> rand_pos(weight.begin(), weight.end());

        // Équilibrage
        int sens,tirage,pos;
        for(long int t = 0; t<LX*T_EQ ; t++){
            tirage = rand_pos(generator);
            pos = syspos[tirage];
            sens   = 2*r01int(generator)-1;
            if(EsosKaw(system,pos,sens,Beta,f)){
                syspos[tirage] = modulo(pos+2*sens,2*LX);
                system[pos]--;
                system[modulo(pos+2*sens,2*LX)]++;
            }
        }
        for(int x=0;x<LX;x++)
            hauteur[x] = system[2*x]+system[2*x+1];
        
    //****** CALCUL FFT  *********//
        double* fft = new double[LX/2];
        double* ffttot = new double[LX/2];
        int var = 0;
        for(int x=0;x<LX/2;x++) ffttot[x] = 0;

        for(long int t = 0; t<LX*T_MAX ; t++){
            tirage = rand_pos(generator);
            pos = syspos[tirage];
            sens   = 2*r01int(generator)-1;

            if(EsosKaw(system,pos,sens,Beta,f)){
                syspos[tirage] = modulo(pos+2*sens,2*LX);
                system[pos]--;
                system[modulo(pos+2*sens,2*LX)]++;
                hauteur[static_cast<int>(pos/2)]--;
                hauteur[modulo(static_cast<int>(pos/2)+sens,LX)]++;
            }
            if(t%LX*T_COR == 0){
                SpCorrel(fft,hauteur,LX,mean);
                gsl_fft_real_radix2_transform(fft,1,LX/2);
                for(int x=0;x<LX/2;x++) ffttot[x]+=fft[x];
                var++;
            }
        }
        str = prefix+"/X"+to_string(LX)+"_mu"+streamM.str()+"_T"+streamT.str()+"_F"+to_string(f);
        ofstream fmag(str.c_str(),std::ofstream::out);

        cout << var << " " << T_MAX/T_COR <<endl;
        for(int x=0;x<LX/2;x++)
            fmag << x << " " << TOT(ffttot,x,LX/2)*T_COR/T_MAX << endl;
        fmag.close();

        delete[] system;
        delete[] syspos;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_H,MPI_DOUBLE,valTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_H,MPI_DOUBLE,eneTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_H,MPI_DOUBLE,ecartTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/X"+to_string(LX)+"_mu"+streamM.str()+"_T"+streamT.str()+"_F"+streamF.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Hn;i++)
            fmag << ((Hn==1)? drive : i*drive/(Hn-1)) << "\t" << valTot[i] <<  "\t\t" <<  ecartTot[i] << "\t\t" << J*eneTot[i]<<"\n";
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

// Diagonalisation de la matrice
double getEquilibrium(int l,double chim){
        EigenSolver<MatrixXd> es;
        MatrixXd transfer(l,l);
        MatrixXd diag(l,l);
        for(int i=0;i<l;i++){
                diag(i,i) = i;
            for(int j=0;j<l;j++){
                transfer(i,j) = exp(-Beta*( J*abs(i-j)- chim*(i+j)/2)- (lnfact(i)+lnfact(j))/2) ;
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

// Dynamique de Kawasaki
bool EsosKaw(int* array,int x, int sens, double kbeta,double f){ 
    if(array[modulo(x,2*LX)] < 1) return false;
    int hn,hx,hxp,hnp;
    if(x%2==0){
        hn  = array[modulo(x-2*sens,2*LX)]+array[modulo(x-2*sens+1,2*LX)];
        hx  = array[modulo(x  ,2*LX)]+array[modulo(x+1  ,2*LX)];
        hxp = array[modulo(x+2*sens,2*LX)]+array[modulo(x+2*sens+1,2*LX)];
        hnp = array[modulo(x+4*sens,2*LX)]+array[modulo(x+4*sens+1,2*LX)];
    }
    else{
        hn  = array[modulo(x-2*sens,2*LX)]+array[modulo(x-2*sens-1,2*LX)];
        hx  = array[modulo(x  ,2*LX)]+array[modulo(x-1 ,2*LX)];
        hxp = array[modulo(x+2*sens,2*LX)]+array[modulo(x+2*sens-1,2*LX)];
        hnp = array[modulo(x+4*sens,2*LX)]+array[modulo(x+4*sens-1,2*LX)];
    }
    int ene = abs(hx-1-hn) + abs(hx-hxp-2)+abs(hxp+1-hnp);
    ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
    double sos = J*ene + f*sens*(x%2==0);

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

// Calcul de la fonction d'auto-corrélation
void FnCorrel(double* cor, double* val, int tmax, int tcor){
    for(int t=0;t<tcor;t++){
        double a=0,b=0,c=0;
        for(int tp=0;tp<tmax;tp++){
            a += val[tp]*val[modulo(tp+t,tcor)] ;
            b += val[tp];
            c += val[modulo(tp+t,tcor)];
        }
        cor[t] = a / static_cast<double>(tmax) - b * c / static_cast<double>(pow(tmax,2)) ;
    }
}
// Calcul de la fonction de corrélation spatiale
void SpCorrel(double*cor,double*val,int lmax,double mbar){
    double m2 = mbar*mbar;
    for(int l=0;l<lmax/2;l++){
        double a=0;
        for(int lp=0;lp<lmax;lp++)
            a += val[lp]*val[modulo(lp+l,lmax)]-m2 ;
        cor[l] = a / static_cast<double>(lmax);
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

