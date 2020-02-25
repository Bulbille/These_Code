//Basic stuff
#include <mpi.h>
#include <sys/file.h> //Files lock
#include <random> 
#include <iostream> //cout
#include <sstream> //stringstream
#include <iomanip> // setprecision
#include <fstream> //files

#define MODEL 'A'
using namespace std;
const double T_C    = 2./log(1.+sqrt(2)), J = 1;
double  ttc = 1 ,Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
const int       bitLX = 7, LX = 2<<bitLX ; 
//LX doit être une puissance de 2 
// 2<<0 = 2, 2<<n = 2^(n+1)
int Hn = 2e0;
double T_EQ = 1e6, T_MAX = 1e4;
double Hmax = 0.005;
double taux_dyn = 0 ; // 0 : 100% Kawasaki, 1 : 100% Glauber
const double drive =2; 

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

bool EsosKaw(int* array, int x, int ajout,int sens, double kbeta,double Champ,double f){ 
    int hn  = array[modulo(x-sens,LX)],
        hx  = array[modulo(x  ,LX)],
        hxp = array[modulo(x+sens,LX)],
        hnp = array[modulo(x+2*sens,LX)],
        hxf = hx+ajout,
        hxpf= hxp-ajout;
    if(hxf < 0  or hxpf < 0) return false;
    double ef = abs(hn-hxf)  + abs( hxf-hxpf) + abs(hxpf-hnp);
    double ei = abs(hn-hx)  + abs(hx-hxp) + abs(hxp-hnp);
    double sos = J*(ef-ei) - f*sens;

    double D_e = exp(-kbeta*sos);
    return (D_e > 1) ? true : (rand_01(generator) <  D_e );
}


/***********************************/

int main(int argc,char* argv[]){
    parametres(argc,argv);

    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
    stringstream stream;
    string str;
    stream << fixed << setprecision(4) << taux_dyn;

    while(Hn%nbprocs != 0) Hn++;
    int Nb_H = Hn/nbprocs;
    double* valeurs = new double[Nb_H];
    double* energies = new double[Nb_H];
    double* correlPar = new double[Nb_H];
    double* valtot = new double[Hn];
    double* enetot = new double[Hn];
    double* ecartTyp = new double[Nb_H];
    double* ecartTot = new double[Hn];
    double* correlTot = new double[Hn];

    /***** Parallélisation sur les champs magnétiques ****/
    for(int k = 0; k<Nb_H ; k++){
        double f = 1.*(k+rang*Nb_H)*drive/(Hn-1); //Le -1 est pour arriver à Hmax lorsque k = Nb_H au rang supérieur
        int ajout,sens,tirage,tirage2;
        double phi = 0,phi_p = 0;
        double ene = 0, ene_p = 0;
        double dynamique;
        int* system = new int[LX]; 
        double* mag = new double[static_cast<int>(T_MAX)];
        for(int x = 0;x<LX;x++)
            system[x] = 60;
        ecartTyp[k] = 0;

        /****** Simulation MC ****/
        // Équilibrage
        for(long int t = 0; t<T_EQ ; t++){
            /***** Étape de Monte Carlo *******/
            //tirage  = rand_lx(r01int,generator);
            tirage  = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;

            if(EsosGlau(system,tirage,ajout,Beta,Hmax))
                system[tirage] += ajout;
        }
        for(int x=0;x<LX;x++){
            phi_p += system[x];
            ene_p += abs(system[x]-system[modulo(x+1,LX)]);
        }
        // Simu
        for(long int t = 0; t<T_MAX ; t++){
            /***** Étape de Monte Carlo *******/
            tirage  = rand_lx(generator);
            ajout   = 2*r01int(generator)-1;
            dynamique = rand_01(generator);

            if(dynamique < taux_dyn ) {
                if(EsosGlau(system,tirage,ajout,Beta,Hmax)){
                    ene_p -= abs(system[tirage]-system[modulo(tirage+1,LX)])+abs(system[tirage]-system[modulo(tirage-1,LX)]);
                    system[tirage] += ajout;
                    phi_p += ajout;
                    ene_p += abs(system[tirage]-system[modulo(tirage+1,LX)])+abs(system[tirage]-system[modulo(tirage-1,LX)]);
                }
            }//*
            else{
                sens   = 2*r01int(generator)-1;
                if(EsosKaw(system,tirage,ajout,sens,Beta,Hmax,f)){
                    tirage2 = modulo(tirage+sens,LX);
                    ene_p -= abs(system[tirage]-system[modulo(tirage-sens,LX)])
                            +abs(system[tirage2]-system[modulo(tirage2+sens,LX)])
                            +abs(system[tirage]-system[tirage2]);
                    system[tirage] += ajout;
                    system[tirage2] -= ajout;
                    ene_p += abs(system[tirage]-system[modulo(tirage-sens,LX)])
                            +abs(system[tirage2]-system[modulo(tirage2+sens,LX)])
                            +abs(system[tirage]-system[tirage2]);
                }
            }//*/
            /***** Mesure du système *******/
            phi += phi_p;
            ene += ene_p;
            ecartTyp[k] += pow(phi_p,2);
            double intpart;
            if( modf(t,&intpart) <= 1e-6)
                mag[static_cast<int>(t)] = phi_p/static_cast<double>(LX);

        }
//        stringstream autre;
//        autre << fixed << setprecision(2) << f;
//        string ssys = prefix+"/system"+autre.str();
//        ofstream fsys(ssys.c_str(),std::ofstream::out);
//        for(int x=0;x<LX;x++)
//            fsys << x << " " << system[x] << "\n";
//        fsys.close();




        double* correl = new double[static_cast<int>(T_MAX)];
        for(int t=0;t<T_MAX;t++){
            double a=0,b=0,c=0;
            for(int tp=0;tp<T_MAX-t;tp++){
                a += mag[tp]*mag[tp+t] ;
                b += mag[tp];
                c += mag[tp+t];
            }
            correl[t] = a / static_cast<double>(T_MAX-t) - b * c / static_cast<double>(pow(T_MAX-t,2)) ;
            //    correl[t] = a / static_cast<double>(T_MAX-t) - pow(phi,2) ;
        }
        /*** Intégration pour trouver le temps d'autoccorélation ***/
        double tau = 0;
        for(int t=0;t<T_MAX;t++){
            if(t == 0 or t == T_MAX-1)
                tau += correl[t]/correl[0];
            else if(t % 2 == 0)
                tau += 2*correl[t]/correl[0];
            else
                tau += 4*correl[t]/correl[0];
        }

        delete[] system;
        valeurs[k] = phi / static_cast<double>(LX*T_MAX) ;
        energies[k] = ene / static_cast<double>(LX*T_MAX) ;
        ecartTyp[k] = ecartTyp[k]/static_cast<double>(LX*LX*T_MAX) ;
        correlPar[k] = tau/3;

    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs,Nb_H,MPI_DOUBLE,valtot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_H,MPI_DOUBLE,enetot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_H,MPI_DOUBLE,ecartTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(correlPar,Nb_H,MPI_DOUBLE,correlTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(rang == 0){
        str = prefix+"/"+MODEL+"-X"+to_string(LX)+'F'+stream.str();
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Hn;i++)
            fmag << 1.*i*drive/(Hn-1) << " " << valtot[i] <<  " " << sqrt(ecartTot[i]-pow(valtot[i],2)) << " " << enetot[i] << " " << ecartTot[i] << "  " << correlTot[i] <<"\n";
    }
    MPI_Finalize();
return 0;
}
/***************************************************************/
/************** Fin du main ************************************/
/***************************************************************/

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
bool isOnlyDouble(const char* str){
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

    return true;
}

