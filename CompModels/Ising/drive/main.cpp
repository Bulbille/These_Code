#include "init.h"
#include "math.h"

double normspace(int step,double min,double max, double n){
    return (n == 1) ? max : min+(max-min)/(n-1)*1.*step;
}

void parametres(int argc, char* argv[], int grille[][TAILLE_Y]);
void generation(int grille[][TAILLE_Y],int LX,int LY);

double BETA;
double ttc = 1;
string prefix = ".";
string suffix = "";

std::default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

std::uniform_int_distribution<int> rand_lx(0,TAILLE_X-1);
std::uniform_int_distribution<int> rand_ly(1,TAILLE_Y-2);

/****************** Main ******************/
int main(int argc, char *argv[]){

	int grille[TAILLE_X][TAILLE_Y];
	generator.seed(time(NULL));
	parametres(argc,argv); /**/
    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);


    int Nb_T = Tn/nbprocs;
    double* valeurs = new double[Nb_T];
    double* energies = new double[Nb_T];
    double* valTot = new double[Tn];
    double* eneTot = new double[Tn];

    int cmpt = 0;
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        generation(grille,TAILLE_X,TAILLE_Y);
        cout << k << endl;
        valeurs[cmpt] = 0;
        energies[cmpt] = 0;
        double T =  normspace(k,Tmin,Tmax,Tn);
        Beta = 1./(ttc*T);
        /** Équilibrage **/
        for(long int t = 0; t<=T_EQ*N ; t++){
            int x = rand_lx(generator);
            int y = rand_ly(generator);
            double e_flip = grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][y-1]+grille[x][y+1];
            e_flip *= 2*J*grille[x][y];
            if(rand_01(generator) < exp(-BETA*e_flip)){
                grille[x][y] *= -1;
            }
        } 
        int phi =0, ene = 0;
        for(int x=0;x<TAILLE_X;x++){
            for(int y=1;y<TAILLE_Y-1;y++){
                phi += grille[x][y];
                ene += grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][y-1]+grille[x][y+1];
            }
        }
        /****** DYNAMIQUE DE GLAUBER ***************/
        for(long int t = 0; t<=T_MAX*N ; t++){
            int x = rand_lx(generator);
            int y = rand_ly(generator);

            double e_flip = grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][modulo(y-1,TAILLE_Y)]+grille[x][modulo(y+1,TAILLE_Y)];
            e_flip *= 2*J*grille[x][y];
            if(rand_01(generator) < exp(-BETA*e_flip)){
                grille[x][y] *= -1;
                phi += 2*grille[x][y];
                ene += e_flip;
            }
        } 
        
        valeurs[cmpt] /= static_cast<double>(N*T_MAX) ;
        energies[cmpt] /= static_cast<double>(N*T_MAX) ;
        cmpt++;
    }

    /*** Récupération des calculs parallélisés ****/
    MPI_Gather(valeurs ,Nb_T,MPI_DOUBLE,valTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(energies,Nb_T,MPI_DOUBLE,eneTot,Nb_T,MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*** Écriture dans fichier ***/
    if(rang == 0){
        str = prefix+"/dataIsing";
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int i=0;i<Tn;i++)
            fmag << normspace(i,Tmin,Tmax,Tn)<< "\t" << valTot[i] << "\t\t" << J*eneTot[i]<<  "\n";
        fmag.close();
    }
    MPI_Finalize();
    delete[] valeurs;
    delete[] energies;
    delete[] valTot;
    delete[] eneTot;
    return 0;
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[], int grille[][TAILLE_Y]){
	BETA=1./T_C;
        if(argc > 1){ 
		for(int i = 1; i<argc; ++i){
			std::string arg = argv[i];
			if(arg == "--t"){
				if(isOnlyDouble(argv[i+1])){
					ttc = atof(argv[i+1]);
					BETA = 1/(ttc*T_C);
				}
				else{ cout << "Le paramètre temperature n'est pas un double"; abort();}
			}
			else if(arg == "--prefix"){
				prefix = argv[i+1];
			}
			else if(arg == "--suffix"){
				suffix = argv[i+1];
			}
		}
	}
	generation(grille,TAILLE_X,TAILLE_Y);
}
void generation(int grille[][TAILLE_Y],int LX,int LY){
    for(int x=0;x<LX;x++){
        grille[x][0] = -1;
        grille[x][LY-1] = 1;
        for(int y=1;y<LY-1;y++)
            grille[x][y] = (rand_01(generator) < 0.5 ) ? grille[x][y] = 1 : grille[x][y] = -1 ;
    }
}
