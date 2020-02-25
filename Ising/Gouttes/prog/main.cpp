#include "init.h"
#include "math.h"
#include "dynamique.h"

#include "time.h"
clock_t cputime = clock();

void parametres(int argc, char* argv[], int grille[][TAILLE_Y]);
void ecriture(int grille[][TAILLE_Y], double mag_x[], double e_x[], double e_y[], double magnetisation, double e_int, double chi, double c_v,double temps);

double BETA;
double F=F0;
double w=w0;
double energie_totale = 0;
double ttc = 1;
double H=H0;
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

	/******* VARIABLES THERMO *********/
        double magnetisation, e_int, mag2, e2;
                magnetisation=0; e_int=0;mag2=0;e2=0;
        double mag_x[TAILLE_X];
        for(int i=0;i<TAILLE_X;i++)
                mag_x[i]=0;

        double e_x[TAILLE_Y];
        double e_y[TAILLE_Y];
        for(int j=0;j<TAILLE_Y;j++){
                e_x[j]=0; e_y[j]=0;
        }


	/**** VARIABLES INTERMÉDIAIRES *****/
	int nb_liens;
        int liens[TAILLE_X * TAILLE_Y * 5];
        double acceptance_rate[TAILLE_X * TAILLE_Y * 2];
	double T_TOT = 0;
	double delta_t=0;
	int nb_photos=-1;



	/***** INITIALISATION DU SYSTÈME ***********/
	int mag_i = 0;
	for(int j=0;j<TAILLE_Y;j++){
		for(int i=0;i<TAILLE_X;i++)
			mag_i += grille[i][j];
	}
	/****************/
	/*** Glauber ****/
	/****************/
	for(double i=0;i<T_THERM;i+=1./N){
	//	if(mag_i < -0.4*N || mag_i > 0.4*N) break;
		int x = rand_lx(generator);
		int y = rand_ly(generator);
		double e_flip = grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][modulo(y-1,TAILLE_Y)]+grille[x][modulo(y+1,TAILLE_Y)];
		e_flip *= 2*J*grille[x][y];
		if(rand_01(generator) < exp(-BETA*e_flip)){
			grille[x][y] *= -1;
			mag_i += 2*grille[x][y];
		}
	} //*
	//Si les -1 sont majoritaires, on va renverser le système
	/*
	{
	int mag = 0;
	for(int j=0;j<TAILLE_Y;j++){
		for(int i=0;i<TAILLE_X;i++)
			mag += grille[i][j];
	}	
	if(mag < 0){
		for(int j=0;j<TAILLE_Y;j++){
			for(int i=0;i<TAILLE_X;i++)
				grille[i][j] *= -1;
		}	
	}
	}	//*/


	/****************/
	/**** Kawasaki **/
	/****************/
	for(double i=0;i<T_CHAMP;i+=delta_t){
		maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens);
		delta_t = 1./delta_t;
		int tirage = tirage_lien(acceptance_rate, nb_liens);
		int x1 = liens[tirage*4+0];
		int y1 = liens[tirage*4+1];
		int x2 = liens[tirage*4+2];
		int y2 = liens[tirage*4+3];
		grille[x1][y1] *= -1;
		grille[x2][y2] *= -1;
	} 
        ecriture(grille,mag_x, e_x, e_y,mag_i, 0, 0, 0,-1);
	H=0;
	/********/


	/********************************************/
	/****** DYNAMIQUE DE KAWASAKI ***************/
	/*******************************************/
	for(double t = 0; t<=T_MAX ; t+= delta_t){
		
		/****** MONTE CARLO STEP ****/
		/****************/
		/**** Kawasaki **/
		/****************/
		{
		delta_t = 0;
		maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens);
		delta_t = 1./delta_t;
		T_TOT += delta_t;

		int tirage = tirage_lien(acceptance_rate, nb_liens);
		int x1 = liens[tirage*4+0];
		int y1 = liens[tirage*4+1];
		int x2 = liens[tirage*4+2];
		int y2 = liens[tirage*4+3];
		grille[x1][y1] *= -1;
		grille[x2][y2] *= -1;
		}

		/****** CALCULS GRANDEURS THERMO **************/

		double mag_conf = 0;
		double e_conf= 0;

		for(int y=0;y<TAILLE_Y;y++){
		for(int x=0;x<TAILLE_X;x++){
			mag_conf += grille[x][y];
			mag_x[x] += 1.*grille[x][y]*delta_t;
			e_conf += hamil_loc(grille,x,y)/2;
			if(y < TAILLE_Y-1 && y > 0){ 
				e_y[y] -= 1.*grille[x][y]*grille[x][y+1]*delta_t;
				e_x[y] -= 1.*grille[x][y]*grille[modulo(x+1, TAILLE_X)][y]*delta_t; 
			}
		}
		}
		magnetisation += abs(mag_conf)*delta_t/N;
		e_int += e_conf/N*delta_t;
		mag2 += mag_conf/N*mag_conf/N*delta_t;
		e2 += e_conf/N*e_conf/N*delta_t;
	
                int modulo_photo = static_cast<int>(modulo_d(t,T_PHOTO));
                if(t > 0 && modulo_photo == 0 && static_cast<int>(t/T_PHOTO) != nb_photos){
                        double p_mag = magnetisation/T_TOT;
                        double p_e = e_int/T_TOT;
                        double p_e2 = e2/T_TOT;
                        double p_m2 = mag2/T_TOT;
                        double p_magx[TAILLE_X], p_ex[TAILLE_Y], p_ey[TAILLE_Y];
                        for(int j=0; j<TAILLE_Y;j++){
                                p_ex[j] = e_x[j] / (TAILLE_X*T_TOT);
                                p_ey[j] = e_y[j] / (TAILLE_X*T_TOT);
                        }
			for(int x=0;x<TAILLE_X;x++) p_magx[x] = mag_x[x] / (TAILLE_Y*T_TOT);
				
                        double chi = BETA*(p_m2 - pow(p_mag,2))*N;
                        double c_v = pow(BETA,2)*(p_e2 - pow(p_e,2))*N;
                        ecriture(grille,p_magx, p_ex, p_ey,p_mag, p_e, chi, c_v,static_cast<int>(t));
                        nb_photos++;
                }
	}
	/*******************************************/
	/*************** FIN  **********************/
	/*******************************************/

	/****** ÉCRITURE DANS FICHIERS ******/
        magnetisation /= T_TOT;
        e_int /= T_TOT;
        e2 /= T_TOT;
        mag2 /= T_TOT;
	for(int j=0; j<TAILLE_Y;j++){
		 e_x[j] /= (TAILLE_X*T_TOT);
		 e_y[j] /= (TAILLE_X*T_TOT);
	}
	for(int i=0; i<TAILLE_X;i++) mag_x[i] /= (TAILLE_Y*T_TOT);
        double chi = BETA*(mag2 - pow(magnetisation,2))/N;
        double c_v = pow(BETA,2)*(e2 - pow(e_int,2))/N;

	ecriture(grille, mag_x, e_x, e_y, magnetisation/N, e_int/N, chi, c_v/N,T_TOT);

return 0;
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[], int grille[][TAILLE_Y]){
	bool fichier_depart = false;
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
			else if(arg == "--w0"){
				if(isOnlyDouble(argv[i+1])){
					w = atof(argv[i+1]);
				}
				else{ cout << "Le paramètre w0 n'est pas un double"; abort();}
			}
			else if(arg == "--f0"){
				if(isOnlyDouble(argv[i+1])){
					F = atof(argv[i+1]);
				}
				else{ cout << "Le paramètre F0 n'est pas un double"; abort();}
			}
			else if(arg == "--h"){
				if(isOnlyDouble(argv[i+1])){
					H = atof(argv[i+1]);
				}
				else{ cout << "Le paramètre h n'est pas un double"; abort();}
			}
			else if(arg == "--reprendre"){
				std::ifstream depart(argv[i+1]);
				if(static_cast<bool>(depart)){
					std::string nom_fichier=argv[i+1];
					generation(grille, nom_fichier);
					fichier_depart = true;
				}
				else
					cout << "Le fichier entré en paramètre n'a pas été trouvé \n";
				depart.close();
			}
			else if(arg == "--prefix"){
				prefix = argv[i+1];
			}
			else if(arg == "--suffix"){
				suffix = argv[i+1];
			}
		}
	}
	if(!fichier_depart)
		generation(grille);
}

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(int grille[][TAILLE_Y], double mag_x[], double e_x[], double e_y[], double magnetisation, double e_int, double chi, double c_v,double temps){

        /*********** CREATION DOSSIER POUR RESULTATS ********/

        if(system(("mkdir -p " + prefix).c_str())) {;}
        string str = to_string(ttc);
	str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
        string algo = prefix + "/kaw_ttc_" + str + suffix + "-t" + to_string(temps);
	if(fabs(w) > 1e-4){
		str = to_string(w);
		str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
		algo += "-w0-" + str;
	}
	if(fabs(F) > 1e-4){
		str = to_string(F);
		str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
		algo += "-f0-" + str;
	}
	if(fabs(H) > 1e-4){
		str = to_string(H);
		str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
		algo += "-h-" + str;
	}

	/******* ÉCRITURES GRANDEURS THERMO ****/
        str = prefix + "/thermo_kaw" + suffix;
        ifstream test(str.c_str());
        ofstream thermo_res;
        if(static_cast<bool>(test)){
                thermo_res.open(str.c_str(), ios::out | ios::app);
        }
        else{
                thermo_res.open(str.c_str(), ios::out);
        }

        thermo_res << ttc  << " " << magnetisation << " " << e_int << " " << chi << " " << c_v << " " << temps << " " << float(clock()-cputime)/CLOCKS_PER_SEC << "\n";

        /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
        str = algo;
        ofstream result(str.c_str());

	for(int y=0; y<TAILLE_Y;y++){
		for(int x=0; x<TAILLE_X;x++){
			result << x << " " << y << " " << grille[x][y] << "\n";
		}//*/
	}
        result.close();

        /**** ECRITURE DANS FICHIER GRANDEURS EN FONCTION DE Y **/
        str = algo +  "-mag_x";
        ofstream res_mag_x(str.c_str());

        for(int x=0; x<TAILLE_X;x++){
                res_mag_x <<  x << " " << mag_x[x] << "\n";
        }
        res_mag_x.close();

        /**** ECRITURE DANS FICHIER ENERGIE DE BANDE **/
        str = algo + "-energie";
        ofstream energie_bande(str.c_str());

        for(int y=1; y<TAILLE_Y-1;y++){
                energie_bande <<  y << " " << e_x[y] << " " << e_y[y]<< "\n";
        }
        energie_bande.close();
}

