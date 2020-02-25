//Basic stuff
#include <mpi.h>
#include <sys/file.h> //Files lock
#include <random> 
#include <iostream> //cout
#include <sstream> //stringstream
#include <iomanip> // setprecision
#include <fstream> //files
#define MODEL 'B'

using namespace std;
//Propriétés de l'Hamiltonien
double T_C    = 2./log(1.+sqrt(2)), J = 1.0, ttc = 1, Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
int       LX    = 70, LY = 30;
double Hmax = 0.0;
long int  T_EQ = 1e4,  T_MAX = 1e4;
//Générateur nombres aléatoires
default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

/***********************************/
/**** Définitions des fonctions ****/
int modulo(int a, int b); 
void ecriture(int L,double energie,double champ, string str,int* histo);
bool isOnlyDouble(const char* str);
void parametres(int argc, char* argv[]);

bool Esos(int* array, int x, int ajout,double kbeta,double Champ,int L){
    int hx  = array[x],
        hxp = array[modulo(x+1,LX)],
        hxm = array[modulo(x-1,LX)],
        hx2 = hx+ajout;
#if MODEL == 'A'
    double champ = Champ*( (L-hx2) - (L-hx) );
#elif MODEL == 'B'
    double champ = Champ*( abs(L-hx2) - abs(L-hx) );
#elif MODEL == 'C'
    double champ = -Champ*( abs(L-hx2) - abs(L-hx) );
#endif
    double sos    = abs(hx2-hxm) + abs(hx2-hxp) - (abs(hx - hxm)  +  abs(hx - hxp) );

    double D_e   = (hx2 <= 2*L and hx2 >= 0) * exp(-kbeta*(sos+champ));
    return (rand_01(generator) <  D_e );
}
/***********************************/

int main(int argc,char* argv[]){
    std::uniform_int_distribution<int> rand_lx(0,LX-1);
    parametres(argc,argv);
    generator.seed(time(NULL));

    int ajout,tirage;
    double phi_p = 0,phi=0;
    int* system = new int[LX]; 
    double* mag = new double[T_MAX];
    for(int x = 0;x<LX;x++)
        system[x] = LY;

    /*** Équilibrage du système ***/
    for(double t = 0; t<T_EQ ; t+=1./LX){
        tirage  = rand_lx(generator);
        ajout   = -1+2*round(rand_01(generator));
        if(Esos(system,tirage,ajout,Beta,Hmax,LY)){
            system[tirage] += ajout;
        }
    }
    for(int x = 0;x<LX;x++){
#if MODEL == 'A'
        phi_p += (LY-system[x]); 
#elif MODEL == 'B'
        phi_p += abs(LY-system[x]); 
#elif MODEL == 'C'
        phi_p -= abs(LY-system[x]); 
#else
        cout << "Modèle non choisi";
        exit();
#endif
    }

    for(double t = 0; t<T_MAX ; t+=1./LX){
        /***** Étape de Monte Carlo *******/
        tirage  = rand_lx(generator);
        ajout   = -1+2*round(rand_01(generator));
        if(Esos(system,tirage,ajout,Beta,Hmax,LY)){
            system[tirage] += ajout;
#if MODEL == 'A'
            phi_p += (LY-system[tirage])-(LY-(system[tirage]-ajout));
#elif MODEL == 'B'
            phi_p += abs(LY-system[tirage])-abs(LY-(system[tirage]-ajout));
#elif MODEL == 'C'
            phi_p -= abs(LY-system[tirage])-abs(LY-(system[tirage]-ajout));
#endif
        }
        phi += phi_p;
        double intpart;
        if( modf(t,&intpart) <= 1e-6)
            mag[static_cast<int>(t)] = phi_p/static_cast<double>(LX);
    }

    /** Fonction d'autocorrélation **/
    phi /=static_cast<double>(LX*T_MAX*LX);
    double* correl = new double[T_MAX];
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
    tau /= 3;
    cout << tau << endl;

    /*
    string str = prefix+"/"+MODEL+"corre"+to_string(T_MAX);
    ofstream fcorr(str.c_str(),std::ofstream::out);
    for(int t=0;t<T_MAX;t++){
        fcorr << t << " " << correl[t]/correl[0] << "\n";
    }*/

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
            else if(arg == "-lx"){
                if(isOnlyDouble(argv[i+1]))
                    LX = static_cast<int>(atof(argv[i+1]));
            }
            else if(arg == "-ly"){
                if(isOnlyDouble(argv[i+1]))
                    LY = static_cast<int>(atof(argv[i+1]));
            }
            else if(arg == "-suffix"){
                suffix = argv[i+1];
            }
        }
    }
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

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(int L, double energie,double champ, string str,int* histo){
    /*********** CREATION DOSSIER POUR RESULTATS ********/
    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    string name = prefix+'/'+str+suffix;
    FILE* fenergie;
    fenergie = fopen(name.c_str(),"a+");
    flock(fileno(fenergie),LOCK_EX);
    fprintf(fenergie,"%i %f %f \n",L,champ,energie);
    flock(fileno(fenergie),LOCK_UN);
    fclose(fenergie);

    algo = prefix +"/histo";
    {
        stringstream stream;
        stream << fixed << setprecision(2) << ttc ;
        algo += "-t"+stream.str();
    }
    {
        stringstream stream;
        stream << LY;
        algo += "-ly" + stream.str();
    }
    {
        stringstream stream;
        stream << fixed << setprecision(4) << champ;
        algo += "-h" + stream.str();
    }
    ofstream fhisto(algo.c_str());
    for(int y=0;y<2*LY+1;y++)
        fhisto << y << " " << histo[y] << "\n";

}
