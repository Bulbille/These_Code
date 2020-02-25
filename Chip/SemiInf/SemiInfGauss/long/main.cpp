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
int  Hn = 1e2;
int  Ndiv = 100, pas = Hn/Ndiv;
double T_EQ = 3e6, T_MAX = 3e8;
double Hmax = 1e-5;
double taux_dyn = 0 ; // 0 : 100% Kawasaki, 1 : 100% Glauber
const double drive =10; 
// extremely fast random number generator that also produces very high quality random.
#define UNLIKELY(x) __builtin_expect((x), 0)
// original source : 
// https://martin.ankerl.com/2018/12/08/fast-random-bool/?fbclid=IwAR2r38PZymhkLk_s4JCnXBy1Gi4Fa3nPTnFsekrYnhjSdX2qpDDc_HiCWL8
class sfc64 {
    public:
        using result_type = uint64_t;

        static constexpr uint64_t(min)() { return 0; }
        static constexpr uint64_t(max)() { return UINT64_C(-1); }

        sfc64() : sfc64(std::random_device{}()) {}

        explicit sfc64(uint64_t seed) : m_a(seed), m_b(seed), m_c(seed), m_counter(1) {
            for (int i = 0; i < 12; ++i) {
                operator()();
            }
        }

        uint64_t operator()() noexcept {
            auto const tmp = m_a + m_b + m_counter++;
            m_a = m_b ^ (m_b >> right_shift);
            m_b = m_c + (m_c << left_shift);
            m_c = rotl(m_c, rotation) + tmp;
            return tmp;
        }

    private:
        template <typename T> T rotl(T const x, int k) { return (x << k) | (x >> (8 * sizeof(T) - k)); }

        static constexpr int rotation = 24; 
        static constexpr int right_shift = 11; 
        static constexpr int left_shift = 3;
        uint64_t m_a;
        uint64_t m_b;
        uint64_t m_c;
        uint64_t m_counter;
};

// Unbiased, fastest variant. Gets rid of the counter by sacrificing 1 bit of randomness.
template <typename U = uint64_t> class RandomizerWithShiftT {
    public:
        template <typename Rng> bool operator()(Rng &rng) {
            if (UNLIKELY(1 == m_rand)) {
                m_rand = std::uniform_int_distribution<U>{}(rng) | s_mask_left1;
            }
            bool const ret = m_rand & 1;
            m_rand >>= 1;
            return ret;
        }

    private:
        static constexpr const U s_mask_left1 = U(1) << (sizeof(U) * 8 - 1);
        U m_rand = 1;
};
template <typename U = uint64_t> class RandomizerInt {
    public :
        RandomizerInt(int size) : m_size(size) {}; 

        template <typename rand, typename Rng> U operator()(rand &boolrand, Rng &rng){
            int res = 0;
            for(int i = 0 ; i<m_size;i++)
                res += boolrand(rng) << i ; 
            return res;
        }
    private :
        int m_size;
};

sfc64 generator;
RandomizerWithShiftT<uint64_t> r01int;
//RandomizerInt<int> rand_lx(bitLX);
std::uniform_int_distribution<int> rand_lx(0,LX-1);
std::uniform_real_distribution<double> rand_01(0.0,1.0);

/***********************************/
/**** Définitions des fonctions ****/
int modulo(int a, int b); 
bool isOnlyDouble(const char* str);
void parametres(int argc, char* argv[]);
double getEquilibrium(int l,double h);
bool EsosGlau(int* array, int x, int ajout,double kbeta,double Champ);
bool EsosKaw(int* array, int x, int ajout,int sens, double kbeta,double Champ,double f); 



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
        system[0] = getEquilibrium(LY,Hmax);
        for(int x = 0;x<LX;x++)
            system[x] = system[0];

        /****** Simulation MC ***/
        // Équilibrage
        for(long int t = 0; t<LX*T_EQ ; t++){
            /***** Étape de Monte Carlo *******/
            tirage  = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,Hmax))
                system[tirage] += ajout;
        }
        for(int x=0;x<LX;x++){
            phi += system[x];
            phi2 += pow(system[x],2);
            phi3 += pow(system[x],3);
            ene += pow(system[x]-system[modulo(x+1,LX)],2);
        }
        // Simu
        for(long int t = 0; t<T_MAX ; t++){
            /***** Étape de Monte Carlo *******/
            tirage  = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;
            dynamique = rand_01(generator);

            if(dynamique < taux_dyn ) {
                if(EsosGlau(system,tirage,ajout,Beta,Hmax)){
                    ene -= pow(system[tirage]-system[modulo(tirage+1,LX)],2)+pow(system[tirage]-system[modulo(tirage-1,LX)],2);
                    phi += ajout;
                    phi2 += 2*ajout*system[tirage]+1;// 1 = ajout^2 
                    phi3 += 3*ajout*pow(system[tirage],2)+3*system[tirage]+ajout; // ajout = ajout^3

                    system[tirage] += ajout;
                    ene += pow(system[tirage]-system[modulo(tirage+1,LX)],2)+pow(system[tirage]-system[modulo(tirage-1,LX)],2);
                }
            }
            else{
                sens   = 2*r01int(generator)-1;
                if(EsosKaw(system,tirage,ajout,sens,Beta,Hmax,f)){
                    tirage2 = modulo(tirage+sens,LX);
                    ene -= pow(system[tirage]-system[modulo(tirage-sens,LX)],2)
                            +pow(system[tirage2]-system[modulo(tirage2+sens,LX)],2)
                            +pow(system[tirage]-system[tirage2],2);
                    phi2 += 2*ajout*(system[tirage]-system[tirage2])+2;// 1 = ajout^2 
                    phi3 += 3*ajout*pow(system[tirage],2)+3*system[tirage]+ajout; // ajout = ajout^3
                    phi3 += -3*ajout*pow(system[tirage2],2)+3*system[tirage2]-ajout; // ajout = ajout^3

                    system[tirage] += ajout;
                    system[tirage2] -= ajout;
                    ene += pow(system[tirage]-system[modulo(tirage-sens,LX)],2)
                            +pow(system[tirage2]-system[modulo(tirage2+sens,LX)],2)
                            +pow(system[tirage]-system[tirage2],2);
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
        skew[cmpt] /= static_cast<double>(LX*LX*LX*T_MAX) ;
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
            troisieme = (troisieme-3*val*pow(sigma,2)-pow(val,3))/pow(sigma,3);
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

//        stringstream autre;
//        autre << fixed << setprecision(2) << f;
//        string ssys = prefix+"/system"+autre.str();
//        ofstream fsys(ssys.c_str(),std::ofstream::out);
//        for(int x=0;x<LX;x++)
//            fsys << x << " " << system[x] << "\n";
//        fsys.close();
/***************************************************************/
/************** Fin du main ************************************/
/***************************************************************/

// Diagonalisation de la matrice
double getEquilibrium(int l,double h){
        EigenSolver<MatrixXd> es;
        MatrixXd transfer(l,l);
        MatrixXd diag(l,l);
        for(int i=0;i<l;i++){
                diag(i,i) = i;
            for(int j=0;j<l;j++){
                transfer(i,j) = exp(-Beta*( h*(i+j)/2 + J*pow(i-j,2)) );
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
    double sos    = J*( pow(hx2-hxm,2) + pow(hx2-hxp,2) - (pow(hx - hxm,2)  +  pow(hx - hxp,2) ) );

    double D_e   =  exp(-kbeta*(sos+champ));
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}

bool EsosKaw(int* array, int x, int ajout,int sens, double kbeta,double Champ,double f){ 
    int hn  = array[modulo(x-sens,LX)],
        hx  = array[modulo(x  ,LX)],
        hxp = array[modulo(x+sens,LX)],
        hnp = array[modulo(x+2*sens,LX)],
        hxf = hx+ajout,
        hxpf= hxp-ajout;
    if(hxf < 0  or hxpf < 0) return false;
    double ef = pow(hn-hxf,2)  + pow( hxf-hxpf,2) + pow(hxpf-hnp,2);
    double ei = pow(hn-hx,2)  + pow(hx-hxp,2) + pow(hxp-hnp,2);
    double sos = J*(ef-ei) + f*sens;

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}


/******* Fonctions utilitaires *****/
void parametres(int argc, char* argv[]){
    if(argc > 1){
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "-t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    Beta = 1/(ttc*T_C);
                }
            }
            else if(arg == "-h"){
                if(isOnlyDouble(argv[i+1]))
                    Hmax = atof(argv[i+1]);
            }
            else if(arg == "-prefix"){
                prefix = argv[i+1];
            }
            else if(arg == "-f"){
                if(isOnlyDouble(argv[i+1]))
                    taux_dyn = static_cast<double>(atof(argv[i+1]));
            }
            else if(arg == "-suffix"){
                suffix = argv[i+1];
            }
        }
    }
    T_MAX *= LX;
    if(taux_dyn >0)
        T_MAX /= taux_dyn;
    else
        T_EQ = 0;
    if(system(("mkdir -p " + prefix).c_str())) {;}
}
/*** Modulo **/
int modulo(int a, int b){
    while(a < 0 or a >= b){
        if(a<0) a+= b;
        else if(a >= b) a-=b;
    }
    return a;
}
/*
int modulo(int a, int b){ 
    if(a>0) return a%b;
    while(a < 0) a+= b;
    return a;
}*/

bool isOnlyDouble(const char* str){
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

    return true;
}

