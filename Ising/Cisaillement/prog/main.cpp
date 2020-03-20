#include "init.h"
#include "math.h"
#include "./prng.h"

double normspace(int step,double min,double max, double n){
    return (n == 1) ? max : min+(max-min)/(n-1)*1.*step;
}

void parametres(int argc, char* argv[]);
void generation(int grille[][TAILLE_Y],int LX,int LY);

double BETA=1,ttc  = 0.9;
string prefix = ".";
string suffix = "";

/****************** Main ******************/
int main(int argc, char *argv[]){

	parametres(argc,argv); /**/
    int rang,nbprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&nbprocs);
    string str;


    int Nb_T = Cn/nbprocs;
    for(int k = rang*Nb_T; k<(rang+1)*Nb_T ; k++){
        int grille[TAILLE_X][TAILLE_Y];
        generation(grille,TAILLE_X,TAILLE_Y);
        double shear =  normspace(k,Cmin,Cmax,Cn);
        cout << k << " " << shear << endl;
        /** Équilibrage **/
        for(long int t = 0; t<=T_EQ*N ; t++){
            int x1 = rand_lx(generator);
            int y1 = rand_ly(generator);
            int dir = rand_direction(generator);
            int x2 = modulo(x1 + (dir == 0) - (dir == 1),TAILLE_X);
            int y2 = y1 + (dir == 2) - (dir == 3);
            if(y2 == 0 || y2 == TAILLE_Y -1 || grille[x1][y1] == grille[x2][y2] ) continue;
            double h_1 = grille[modulo(x1+1,TAILLE_X)][y1]+grille[modulo(x1-1,TAILLE_X)][y1]+grille[x1][y1+1]+grille[x1][y1-1]-grille[x2][y2];
            double h_2 = grille[modulo(x2+1,TAILLE_X)][y2]+grille[modulo(x2-1,TAILLE_X)][y2]+grille[x2][y2+1]+grille[x2][y2-1]-grille[x1][y1];

            double sh = shear*grille[x1][y1]*(y1 == y2)*(y1*2/(TAILLE_Y+1)-1);
            double e_flip = 2*J*grille[x1][y1]*(h_1-h_2)+shear;
           if(exp(-BETA*e_flip) >1 or rand_01(generator) < exp(-BETA*e_flip)){
                grille[x1][y1] *= -1;
                grille[x2][y2] *= -1;
            }
        } 
        cout << "Fin équilibre " <<  k << endl;
        int* mag = new int[TAILLE_Y-2];
        double* totmag = new double[TAILLE_Y-2];
        long int* tpsmag = new long int[TAILLE_Y-2];
        for(int y=1;y<TAILLE_Y-1;y++){
            mag[y] = 0;
            totmag[y] = 0;
            tpsmag[y] = 0;
            for(int x=0;x<TAILLE_X;x++){
                mag[y] += grille[x][y];
                totmag[y] += grille[x][y];
            }
        }
        /****** DYNAMIQUE DE GLAUBER ***************/
        for(long int t = 0; t<=T_MAX*N ; t++){
            int x1 = rand_lx(generator);
            int y1 = rand_ly(generator);
            int dir = rand_direction(generator);
            int x2 = modulo(x1 + (dir == 0) - (dir == 1),TAILLE_X);
            int y2 = y1 + (dir == 2) - (dir == 3);
            if(y2 == 0 || y2 == TAILLE_Y -1 || grille[x1][y1] == grille[x2][y2] ) continue;
            double h_1 = grille[modulo(x1+1,TAILLE_X)][y1]+grille[modulo(x1-1,TAILLE_X)][y1]+grille[x1][y1+1]+grille[x1][y1-1]-grille[x2][y2];
            double h_2 = grille[modulo(x2+1,TAILLE_X)][y2]+grille[modulo(x2-1,TAILLE_X)][y2]+grille[x2][y2+1]+grille[x2][y2-1]-grille[x1][y1];

            double sh = shear*grille[x1][y1]*(y1 == y2)*(y1*2/(TAILLE_Y+1)-1);
           // double sh = shear*grille[x1][y1]*(y1 == y2)*(y1 - (TAILLE_Y+1)/2.);
            double e_flip = 2*J*grille[x1][y1]*(h_1-h_2)+shear;
            if(exp(-BETA*e_flip) >1 or rand_01(generator) < exp(-BETA*e_flip)){
                grille[x1][y1] *= -1;
                grille[x2][y2] *= -1;
                if(y1 != y2){
                    mag[y1] += 2*grille[x1][y1];
                    mag[y2] += 2*grille[x1][y2];
                    totmag[y1] += mag[y1]*(t-tpsmag[y1]);
                    tpsmag[y1] = t;
                    totmag[y2] += mag[y2]*(t-tpsmag[y2]);
                    tpsmag[y2] = t;
                }
            }
        }

        for(int y=1;y<TAILLE_Y-1;y++)
            totmag[y] = (totmag[y]+mag[y]*(T_MAX*N-tpsmag[y]))/(1.*T_MAX*N*TAILLE_X);
        
        str = prefix+"/data"+to_string(shear);
        ofstream fgrille(str.c_str(),std::ofstream::out);
        for(int x=0;x<TAILLE_X;x++){
            for(int y=0;y<TAILLE_Y;y++){
                fgrille << x << " " << y << " " << grille[x][y] << endl;
            }
        }
        fgrille.close();
        str = prefix+"/mag"+to_string(shear);
        ofstream fmag(str.c_str(),std::ofstream::out);
        for(int y=1;y<TAILLE_Y-1;y++){
            fmag << y << " " << totmag[y] << endl;
        }
        fmag.close();
        delete[] mag;
    }
    MPI_Finalize();

    return 0;
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[]){
        if(argc > 1){ 
		for(int i = 1; i<argc; ++i){
			std::string arg = argv[i];
			if(arg == "--t"){
				if(isOnlyDouble(argv[i+1])){
					ttc = atof(argv[i+1]);
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
    BETA = 1/(ttc*T_C);
}
void generation(int grille[][TAILLE_Y],int LX,int LY){
    for(int x=0;x<LX;x++){
        for(int y=0;y<LY;y++)
            grille[x][y] = (y<TAILLE_Y/2) ? -1 : 1;
    }
}
