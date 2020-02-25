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
const int       bitLX = 7, LX = 2<<bitLX ; 
const int       LY = 200 ; // pour la diagonalisation de la matrice
//LX doit être une puissance de 2 
// 2<<0 = 2, 2<<n = 2^(n+1)
int  Hn = 20;
int  Ndiv = 20, pas = Hn/Ndiv;
double T_EQ = 3e6, T_MAX = 5e7;
double Hmax = 1e-1;
double taux_dyn = 0 ; // 0 : 100% Kawasaki, 1 : 100% Glauber
const double drive =10; 
#include "../prng.h"
#include "../utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
double getEquilibrium(int l,double h);
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ);
bool EsosKaw(int* array, int x, int ajout,int sens, double kbeta,double Champ,double f,double mean); 



/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
    ostringstream streamF,streamH;
    string str;
    streamF << fixed << setprecision(2) << taux_dyn;
    streamH << scientific << setprecision(1) << Hmax;

    int Nb_H = Hn/nbprocs;
    double* valeurs = new double[Nb_H];
    double* energies = new double[Nb_H];
    double* ecartTyp = new double[Nb_H];
    double* skew = new double[Nb_H];
    double* valTot = new double[Hn];
    double* eneTot = new double[Hn];
    double* ecartTot = new double[Hn];
    double* skewTot = new double[Hn];

    /***** Parallélisation sur les champs magnétiques ****/
    int cmpt = 0;
    valeurs[cmpt] = 0;
    energies[cmpt] = 0;
    ecartTyp[cmpt] = 0;
    skew[cmpt] = 0;
    for(int k = rang*Nb_H; k<(rang+1)*Nb_H ; k++){
        double f = static_cast<int>(k/pas)*drive/(Hn/pas-1); //Le -1 est pour arriver à Hmax lorsque k = Nb_H au rang supérieur
        int ajout,sens,tirage,tirage2;
        double phi = 0,phi2 = 0,phi3 = 0, ene = 0;
        double dynamique;
        int* system = new int[LX]; 
        /******** Initialisation du systeme *****/
        double mean = getEquilibrium(LY,Hmax);
        system[0] = getEquilibrium(LY,Hmax);
        for(int x = 0;x<LX;x++){ system[x] = system[0];         phi += system[x]; }
        while(phi < mean*LX){    system[rand_lx(generator)]++;  phi++; }
        // Équilibrage
        for(long int t = 0; t<LX*T_EQ ; t++){
            /***** Étape de Monte Carlo *******/
            tirage  = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;
            if(taux_dyn > 0){
                if(EsosGlau(system,tirage,ajout,Beta,Hmax)) system[tirage] += ajout;
            }
            else{
                sens   = 2*r01int(generator)-1;
                if(EsosKaw(system,tirage,ajout,sens,Beta,Hmax,f,mean)){
                    tirage2 = modulo(tirage+sens,LX);
                    system[tirage] += ajout;
                    system[tirage2] -= ajout;
                }
            }
        }
        phi=0;
        for(int x=0;x<LX;x++){
            phi += system[x];
            phi2 += pow(system[x],2);
            ene += abs(system[x]-system[modulo(x+1,LX)]);
        }
        mean = phi/static_cast<double>(LX);
        for(int x=0;x<LX;x++) phi3 += pow(system[x]-mean,3);

        /****** Simulation MC ***/
        for(long int t = 0; t<T_MAX ; t++){
            /***** Étape de Monte Carlo *******/
            tirage  = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;
            dynamique = rand_01(generator);

            if(dynamique < taux_dyn ) {
                if(EsosGlau(system,tirage,ajout,Beta,Hmax)){
                    ene -= abs(system[tirage]-system[modulo(tirage+1,LX)])+abs(system[tirage]-system[modulo(tirage-1,LX)]);
                    phi += ajout;
                    phi2 += 2*ajout*system[tirage]+1;// 1 = ajout^2 
                    phi3 += 3*ajout*pow(system[tirage],2)+3*system[tirage]+ajout; // ajout = ajout^3

                    system[tirage] += ajout;
                    ene += abs(system[tirage]-system[modulo(tirage+1,LX)])+abs(system[tirage]-system[modulo(tirage-1,LX)]);
                }
            }
            else{
                sens   = 2*r01int(generator)-1;
                if(EsosKaw(system,tirage,ajout,sens,Beta,Hmax,f,mean)){
                    tirage2 = modulo(tirage+sens,LX);
                    ene -= abs(system[tirage]-system[modulo(tirage-sens,LX)])
                            +abs(system[tirage2]-system[modulo(tirage2+sens,LX)])
                            +abs(system[tirage]-system[tirage2]);
                    phi2 += 2*ajout*(system[tirage]-system[tirage2])+2;// 1 = ajout^2 
                    phi3 += 3*( ajout*(pow(system[tirage],2)-pow(system[tirage2],2)) + system[tirage]-system[tirage2]);

                    system[tirage] += ajout;
                    system[tirage2] -= ajout;

                    ene += abs(system[tirage]-system[modulo(tirage-sens,LX)])
                            +abs(system[tirage2]-system[modulo(tirage2+sens,LX)])
                            +abs(system[tirage]-system[tirage2]);
                }
            }//*/
            /***** Mesure du système *******/
            valeurs[cmpt] += phi;
            energies[cmpt] += ene;
            ecartTyp[cmpt] += phi2;
            skew[cmpt] += phi3;
        }

        delete[] system;
        valeurs[cmpt] /= static_cast<double>(LX*T_MAX) ;
        energies[cmpt] /= static_cast<double>(LX*T_MAX) ;
        ecartTyp[cmpt] /= static_cast<double>(LX*T_MAX) ;
        skew[cmpt] /= static_cast<double>(LX*T_MAX) ;
//        cout << valeurs[cmpt] << " " << ecartTyp[cmpt] << " " << pow(valeurs[cmpt],2) << " "  << sqrt(ecartTyp[cmpt]-pow(valeurs[cmpt],2)) << endl;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_H,MPI_DOUBLE,valTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_H,MPI_DOUBLE,eneTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_H,MPI_DOUBLE,ecartTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
//    MPI_Gather(skew,Nb_H,MPI_DOUBLE,skewTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(rang == 0){
        str = prefix+"/X"+to_string(LX)+"_H"+streamH.str()+"_F"+streamF.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Ndiv;i++){
            double val=0,ene=0,sigma=0,troisieme=0;
            for(int p=0;p<pas;p++){
                val+= valTot[i*pas+p]/pas;
                ene+= eneTot[i*pas+p]/pas;
                sigma= ecartTot[i*pas+p]-pow(valTot[i*pas+p],2);
                troisieme += skewTot[i*pas+p]/pas;
            }
            sigma = sqrt(sigma);
//            troisieme = (troisieme-3*val*pow(sigma,2)-pow(val,3))/pow(sigma,3);
            troisieme = (troisieme-6*val)/pow(sigma,3);
            fmag << i*drive/(Hn/pas-1) << "\t" << val <<  "\t" << sigma << "\t" << troisieme << "\t" << ene<<"\n";
        }
    }//*/
    delete[] valeurs;
    delete[] energies;
    delete[] ecartTyp;
    delete[] skew;
    delete[] valTot;
    delete[] eneTot;
    delete[] ecartTot;
    delete[] skewTot;
    MPI_Finalize();
return 0;
}

// Diagonalisation de la matrice
double getEquilibrium(int l,double h){
        EigenSolver<MatrixXd> es;
        MatrixXd transfer(l,l);
        MatrixXd diag(l,l);
        for(int i=0;i<l;i++){
                diag(i,i) = i;
            for(int j=0;j<l;j++){
                transfer(i,j) = exp(-Beta*( h*(abs(i)+abs(j))/2 + J*abs(i-j)) );
            }
        }
        int imax = 0;
        double max = 0;
        es.compute(transfer,true);
        VectorXcd a = es.eigenvalues();
        for(int i=0;i<LY;i++){
            if(a[i].real()>max){
                imax=i;
                max=a[i].real();
            }
        }
    VectorXcd Vmax = es.eigenvectors().col(imax);
    complex<double> res = Vmax.transpose()*diag*Vmax;
    return res.real();
}
/****** Dynamique de MC *****/
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
    if(hx2<0) return false;
    double champ = Champ*ajout;
    double sos    = J*( abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

bool EsosKaw(int* array, int x, int ajout,int sens, double kbeta,double Champ,double f, double mean){ 
    int hn  = array[modulo(x-sens,LX)],
        hx  = array[modulo(x  ,LX)],
        hxp = array[modulo(x+sens,LX)],
        hnp = array[modulo(x+2*sens,LX)],
        hxf = hx+ajout,
        hxpf= hxp-ajout;
    if(hxf < 0  or hxpf < 0) return false;
    double ef = abs(hn-hxf)  + abs( hxf-hxpf) + abs(hxpf-hnp);
    double ei = abs(hn-hx)  + abs(hx-hxp) + abs(hxp-hnp);
    double sos = J*(ef-ei) + (hx == 0 || hxp == 0) * f*sens*abs(hx-mean);

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}


