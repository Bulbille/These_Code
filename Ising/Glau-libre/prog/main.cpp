#include "init.h"
#include "math.h"

void parametres(int argc, char* argv[], int grille[][TAILLE_Y]);
void generation(int grille[][TAILLE_Y],int LX,int LY);

double BETA;
double ttc = 1;
string prefix = ".";
string suffix = "";

std::default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

//0 → L-1 : périodique
//1 → L-2 : fixe
std::uniform_int_distribution<int> rand_lx(0,TAILLE_X-1);
std::uniform_int_distribution<int> rand_ly(0,TAILLE_Y-1);

/****************** Main ******************/
int main(int argc, char *argv[]){

	int grille[TAILLE_X][TAILLE_Y];
	generator.seed(time(NULL));
/********** OPTIONS DEPUIS LE BASH *********/
/* Le programme peut accepter deux options : */
/* --beta xxx */
/* --reprendre fichier */
/* La première option permet de faire la simulation avec une température donnée */
/* La deuxième permet de reprendre la simulation à partir d'un fichier */
/**/	parametres(argc,argv, grille); /**/
/**/    ttc = static_cast<double>(static_cast<int>(ttc*1000))/1000;
/******************************************/
    ostringstream streamT,streamTime;
    streamT << fixed << setprecision(2) << ttc;
    string str;

	/********************************************/
	/****** DYNAMIQUE DE GLAUBER ***************/
	/*******************************************/
	for(int t = 0; t<=T_MAX ; t++){
		int x = rand_lx(generator);
		int y = rand_ly(generator);

		double e_flip = grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][modulo(y-1,TAILLE_Y)]+grille[x][modulo(y+1,TAILLE_Y)];
        e_flip *= 2*J*grille[x][y];
		if(rand_01(generator) < exp(-BETA*e_flip)) grille[x][y] *= -1;

        if(t % T_PHOTO == 0){
            streamTime.str(std::string());
            cout << t << " " << endl;
            streamTime << fixed << setprecision(4) << t/T_PHOTO;
            str = prefix+"/snap"+streamT.str()+"_time"+streamTime.str();
            ofstream fnasp(str.c_str(),std::ostream::out);
            for(int x=0;x<TAILLE_X;x++){
                for(int y=0;y<TAILLE_Y;y++){
                    fnasp << x << " " << y << " " << grille[x][y] << endl;
                }
            }
            fnasp.close();
        }

	} 
        
		
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
        for(int y=0;y<LY;y++)
            grille[x][y] = (rand_01(generator) < 0.5 ) ? grille[x][y] = 1 : grille[x][y] = -1 ;
    }
}
