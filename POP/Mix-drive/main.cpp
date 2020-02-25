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
const int  Hn = 50;
const long int T_EQ = 1e5, T_MAX = 5e7;
double mumax = 4;
// Taux de diffusion des particules A et B
const double Da = 1, Db = 1;
//taux_dyn : pourcentage de particules A dans le total
double taux_dyn =  0.5; 
const double drive =6; 

#include "./prng.h"
#include "./utilitaires.h"
/***********************************/
/**** Définitions des fonctions ****/
double getEquilibrium(int l,double chim);
bool EsosKaw(int* array, int x, int sens, double kbeta,double f);

/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
    ostringstream streamF,streamH,streamM;
    string str;
    streamF << fixed << setprecision(2) << taux_dyn;
    streamM << fixed << setprecision(2) << mumax;

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
    for(int k = rang*Nb_H; k<(rang+1)*Nb_H ; k++){
        cout << k << " " << rang << endl;
        valeurs[cmpt] = 0;
        energies[cmpt] = 0;
        ecartTyp[cmpt] = 0;
        skew[cmpt] = 0;
        double f = static_cast<int>(k)*drive/(Hn-1); //Le -1 est pour arriver à Hmax lorsque k = Nb_H au rang supérieur
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
        int hx,hn,hxp,hnp;
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
        //Calcul des valeurs initiales        
        double phi2 = 0,phi3 = 0, ene = 0;
        for(int x=0;x<LX;x++){
            hx  = system[2*x]+system[2*x+1];
            hxp = system[modulo(2*x+2,2*LX)]+system[modulo(2*x+3,2*LX)];
            phi2 += pow(hx,2);
            ene += abs(hx-hxp);
            phi3 += pow(hx-mean,3);
        }
        /****** Simulation MC ***/
        for(long int t = 0; t<LX*T_MAX ; t++){
            /***** Étape de Monte Carlo *******/
            tirage = rand_pos(generator);
            pos = syspos[tirage];
            sens   = 2*r01int(generator)-1;

            if(EsosKaw(system,pos,sens,Beta,f)){
                syspos[tirage] = modulo(pos+2*sens,2*LX);
                if(pos%2==0){
                    hn  = system[modulo(pos-2*sens,2*LX)]+ system[modulo(pos-2*sens+1,2*LX)];
                    hx  = system[modulo(pos  ,2*LX)]+      system[modulo(pos+1  ,2*LX)];
                    hxp = system[modulo(pos+2*sens,2*LX)]+ system[modulo(pos+2*sens+1,2*LX)];
                    hnp = system[modulo(pos+4*sens,2*LX)]+ system[modulo(pos+4*sens+1,2*LX)];
                }
                else{
                    hn  = system[modulo(pos-2*sens,2*LX)]+ system[modulo(pos-2*sens-1,2*LX)];
                    hx  = system[modulo(pos  ,2*LX)]+      system[modulo(pos-1 ,2*LX)];
                    hxp = system[modulo(pos+2*sens,2*LX)]+ system[modulo(pos+2*sens-1,2*LX)];
                    hnp = system[modulo(pos+4*sens,2*LX)]+ system[modulo(pos+4*sens-1,2*LX)];
                }
                ene -= abs(hx-hn) + abs(hx-hxp)+abs(hxp-hnp);
                ene += abs(hx-1-hn) + abs(hx-hxp-2)+abs(hxp+1-hnp);
                phi2 += 2*(hxp-hx)+2;
                phi3 += 3*( hx*(1-hx)+hxp*(1+hxp) );

                system[pos]--;
                system[modulo(pos+2*sens,2*LX)]++;
            }
            /***** Mesure du système *******/
            energies[cmpt] += ene;
            ecartTyp[cmpt] += phi2;
            skew[cmpt] += phi3;
        }
        delete[] system;
        delete[] syspos;
        valeurs[cmpt] = mean;
        energies[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        ecartTyp[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        skew[cmpt] /= static_cast<double>(LX*LX*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_H,MPI_DOUBLE,valTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_H,MPI_DOUBLE,eneTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_H,MPI_DOUBLE,ecartTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(skew,Nb_H,MPI_DOUBLE,skewTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/X"+to_string(LX)+"_mu"+streamM.str()+"_F"+streamF.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Hn;i++){
            double val=0,ene=0,sigma=0,troisieme=0;
            val= valTot[i];
            ene= J*eneTot[i];
            sigma=sqrt( ecartTot[i]-pow(valTot[i],2));
            troisieme = (skewTot[i]-6*val)/pow(sigma,3);
            fmag << i*drive/(Hn-1) << "\t" << val <<  "\t\t"  <<  ecartTot[i] << "\t\t" << troisieme << "\t\t" << ene<<"\n";
        }
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] energies;
    delete[] ecartTyp;
    delete[] skew;
    delete[] valTot;
    delete[] eneTot;
    delete[] ecartTot;
    delete[] skewTot;
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
