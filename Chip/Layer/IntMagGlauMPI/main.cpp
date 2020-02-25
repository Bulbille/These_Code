//Basic stuff
#include <mpi.h>
#include <sys/file.h> //Files lock
#include <random> 
#include <iostream> //cout
#include <sstream> //stringstream
#include <iomanip> // setprecision
#include <fstream> //files

#define MODEL 'C'

using namespace std;
//Propriétés de l'Hamiltonien
double T_C    = 2./log(1.+sqrt(2)), J = 1.0, ttc = 1, Beta = 1/(ttc*T_C);
string prefix = "./", suffix = "",algo;
//Propriétés de la grille
int       LX    = 70, LY = 30;
int       Hn = 100; double Hmax = 0.3;
long int  T_EQ = 1e4,  T_MAX = 5e6;
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

    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);

    while(Hn%nbprocs != 0) Hn++;
    int Nb_H = Hn/nbprocs;
    double* valeurs = new double[Nb_H];
    double* valtot = new double[Hn];
    double* ecartTyp = new double[Nb_H];
    double* ecartTot = new double[Hn];

    /***** Parallélisation sur les champs magnétiques ****/
    for(int k = 0; k<Nb_H ; k++){
        cout << k << " \n";
        double H = 1.*(k+rang*Nb_H)*Hmax/(Hn-1); //Le -1 est pour arriver à Hmax lorsque k = Nb_H au rang supérieur
        int ajout,tirage;
        double phi = 0,phi_p = 0;
        int* system = new int[LX]; 
        for(int x = 0;x<LX;x++)
            system[x] = LY;
        ecartTyp[k] = 0;

        /*** Équilibrage du système ***
        for(int t = 0; t<T_EQ ; t++){
            tirage  = rand_lx(generator);
            ajout   = -1+2*round(rand_01(generator));
            if(Esos(system,tirage,ajout,Beta,H,LY)){
                system[tirage] += ajout;
            }
        }
        for(int x=0;x<LX;x++){
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
        //*/

        /****** Simulation MC ****/
        for(int t = 0; t<T_MAX ; t++){
            /***** Étape de Monte Carlo *******/
            tirage  = rand_lx(generator);
            ajout   = -1+2*round(rand_01(generator));
            if(Esos(system,tirage,ajout,Beta,H,LY)){
                system[tirage] += ajout;
#if MODEL == 'A'
                phi_p += (LY-system[tirage])-(LY-(system[tirage]-ajout)); 
#elif MODEL == 'B'
                phi_p += abs(LY-system[tirage])-abs(LY-(system[tirage]-ajout)); 
#elif MODEL == 'C'
                phi_p -= abs(LY-system[tirage])-abs(LY-(system[tirage]-ajout)); 
#endif
            }
            /***** Mesure du système *******/
            phi += phi_p;
            ecartTyp[k] += pow(phi_p,2);
        }
        delete[] system;
        valeurs[k] = phi / static_cast<double>(LX*T_MAX) ;
        ecartTyp[k] = ecartTyp[k]/static_cast<double>(LX*T_MAX) ;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs,Nb_H,MPI_DOUBLE,valtot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(ecartTyp,Nb_H,MPI_DOUBLE,ecartTot,Nb_H,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(rang == 0){
        string str = prefix+"/"+MODEL+"-X"+to_string(LX)+'Y'+to_string(LY);
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Hn;i++)
            fmag << 1.*i*Hmax/(Hn-1) << " " << valtot[i] << " " << sqrt((ecartTot[i]-pow(valtot[i],2))*2*5e1*LX/T_MAX) << "\n";
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
    T_MAX *= LX;
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
