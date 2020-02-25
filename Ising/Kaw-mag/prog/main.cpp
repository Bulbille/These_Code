#include "init.h"
#include "math.h"
#include "dynamique.h"

#include "time.h"
#include "sys/time.h"

clock_t cputime = clock();

void parametres(int argc, char* argv[], int grille[][TAILLE_Y]);
void ecriture(int grille[][TAILLE_Y], double mag_x[], double e_x[], double e_y[], double magnetisation, double e_int, double chi, double c_v,int temps, double correl_y[][TAILLE_Y]);

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
	/***** INITIALISATION DU SYSTÈME ***********
	tnt grille_init[TAILLE_X][TAILLE_Y];
	for(int j=0;j<TAILLE_Y;j++){ 
		for( int i=0;i<TAILLE_X;i++)
			grille_init[i][j] = grille[i][j];
	}
	/*/
 	double correl_y[TAILLE_X][TAILLE_Y];
	/********/


	/******* VARIABLES THERMO *********/
        double magnetisation, e_int, mag2, e2;
                magnetisation=0; e_int=0;mag2=0;e2=0;

        double ex[TAILLE_X];
        double ey[TAILLE_X];
        double magx[TAILLE_X];
        for(int x=0;x<TAILLE_X;x++){
                ex[x]=0; ey[x]=0; magx[x]=0;
        }

	cout << L_FAISCEAU << "\n";

	/**** VARIABLES INTERMÉDIAIRES *****/
	int nb_liens;
        int liens[TAILLE_X * TAILLE_Y * 5];
        double acceptance_rate[TAILLE_X * TAILLE_Y * 2];
	double delta_t=0;
	int photo_suiv=0; int nb_photos=0;

	double H_prime = H;
	double F_prime = F;
	H=H0;	
	F=F0;
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
	H=H_prime;
        for(double i=0;i<T_THERM;i+=delta_t){
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
	F=F_prime;


	/********************************************/
	/****** DYNAMIQUE DE KAWASAKI ***************/
	/*******************************************/

	delta_t=1;
	for(double t = 0; t<=T_MAX ; t+= delta_t){
		/******************************************************/
		/*********** CALCULS GRANDEURS THERMO *****************/
		/** moyenne_n = moyenne_(n-1) * (n-1) / n + a_n / n ***/
		/******************************************************/

		int modulo_photo = static_cast<int>(modulo_d(delta_t+t,T_PHOTO));
                if(photo_suiv <= static_cast<int>(t) && modulo_photo == 0) {
                        photo_suiv += T_PHOTO; nb_photos++;
		
			////////////////////////////////////////////
			/************** Debut thermo **************/
			////////////////////////////////////////////
						/*** Grandeurs thermo scalaires ***/	

			double p_mag = 0;
			double p_e = 0;
			for(int x=0;x<TAILLE_X;x++){
			for(int y=0;y<TAILLE_Y;y++){
				p_mag += grille[x][y];
				p_e += hamil_loc(grille,x,y)/2;
			}
			}
				// Moyennes sur toutes les configurations
			magnetisation = 1.*magnetisation*(nb_photos -1)/nb_photos+abs(p_mag)/nb_photos/N;
			e_int = 1.*e_int*(nb_photos-1)/nb_photos+p_e/nb_photos/N;
			//cout << t_step << ' '  << e_int << ' ' << p_e/N << ' ' << magnetisation << ' ' << p_mag/N << '\n';
			mag2 =1.*mag2*(nb_photos-1)/nb_photos+1.*abs(p_mag)*abs(p_mag)/nb_photos/N/N;
			e2 = 1.*e2*(nb_photos-1)/nb_photos+1.*p_e*p_e/nb_photos/N/N;
			double chi = BETA*(mag2 - pow(magnetisation,2))*N;
			double c_v = pow(BETA,2)*(e2 - pow(e_int,2))*N;

			/*** Grandeurs thermo vectorielles ***/
			double p_magx[TAILLE_X], p_ex[TAILLE_X], p_ey[TAILLE_X];
			double p_correl_y[TAILLE_X][TAILLE_Y];
			int n_correl[TAILLE_X][TAILLE_Y];
			for(int x=0;x<TAILLE_X;x++){
				p_magx[x] = 0 ; p_ex[x] = 0; p_ey[x] = 0;
				// Magnetisation et énergies de bandes
				for(int y=0;y<TAILLE_Y;y++){
					p_magx[x] += 1. * grille[x][y] / TAILLE_Y;
					p_ex[x] -= 1. * grille[x][y] * grille[modulo(x+1, TAILLE_X)][y] / TAILLE_Y;
					p_ey[x] -= 1. * grille[x][y] * grille[x][modulo(y+1, TAILLE_Y)] / TAILLE_Y;

					p_correl_y[x][y] = 0;
					n_correl[x][y] = 0;
				} 
				// Fonction de corrélation 
				for(int y=0;y<TAILLE_Y;y++){
					for(int y2=y+1;y2<TAILLE_Y+y;y2++){
						p_correl_y[x][y2-y] += 1.*grille[x][y]*grille[x][modulo(y2,TAILLE_Y)];
						n_correl[x][y2-y]++;
					}
				}
			}

			for(int x=0;x<TAILLE_X;x++){
				magx[x] = 1.*magx[x]*(nb_photos-1)/nb_photos+1.*p_magx[x]/nb_photos;
				ex[x] = 1.*ex[x]*(nb_photos-1)/nb_photos+1.*p_ex[x]/nb_photos;
				ey[x] = 1.*ey[x]*(nb_photos-1)/nb_photos+1.*p_ey[x]/nb_photos;
				for(int y=0; y<TAILLE_Y;y++){
					p_correl_y[x][y] = 1.*p_correl_y[x][y]/n_correl[x][y] - p_magx[x]*p_magx[x]; 
						//p_magx est la valeur moyenne de la magnétisation par rapport aux y : OK
					correl_y[x][y] = 1.*correl_y[x][y]*(nb_photos-1)/nb_photos+1.*p_correl_y[x][y]/nb_photos;
				}
			}
			ecriture(grille,magx, ex, ey,magnetisation, e_int, chi, c_v,static_cast<int>(t), correl_y);
		}
		////////////////////////////////////////////
		/************** Fin thermo ****************/
		////////////////////////////////////////////

		/****** MONTE CARLO STEP ****/
		delta_t = 0;
		maj_liens(grille, liens, acceptance_rate,&delta_t, &nb_liens);
		delta_t = 1./delta_t;

		int tirage = tirage_lien(acceptance_rate, nb_liens);
		int x1 = liens[tirage*4+0];
		int y1 = liens[tirage*4+1];
		int x2 = liens[tirage*4+2];
		int y2 = liens[tirage*4+3];
		grille[x1][y1] *= -1;
		grille[x2][y2] *= -1;
		//*/
	}

/*************** FIN  **********************/
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
			else if(arg == "--w"){
				if(isOnlyDouble(argv[i+1])){
					w = atof(argv[i+1]);
				}
				else{ cout << "Le paramètre w0 n'est pas un double"; abort();}
			}
			else if(arg == "--f"){
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
void ecriture(int grille[][TAILLE_Y], double mag_x[], double e_x[], double e_y[], double magnetisation, double e_int, double chi, double c_v,int temps, double correl_y[][TAILLE_Y]){

        /*********** CREATION DOSSIER POUR RESULTATS ********/

        if(system(("mkdir -p " + prefix).c_str())) {;}
        string str = to_string(ttc);
	str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
        string algo = prefix + "/kaw_ttc_" + str + suffix + "-t" + to_string(0);
	if(fabs(w) > 1e-4){
		str = to_string(w);
		str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
		algo += "-w-" + str;
	}
	if(fabs(F) > 1e-4){
		cout << "Coucou";
		str = to_string(F);
		str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
		algo += "-f-" + str;
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
	thermo_res.close();

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

        for(int x=1; x<TAILLE_X-1;x++){
                energie_bande <<  x << " " << e_x[x] << " " << e_y[x]<< "\n";
        }
        energie_bande.close();

        /**** ECRITURE DANS FICHIER CORRELATION EN X **/
        str = algo + "-correl_y";
        ofstream correl(str.c_str());

        for(int x=0; x<TAILLE_X;x++){
		for(int y=0;y<TAILLE_Y;y++)
	                correl <<  x << " " << y << " " << correl_y[x][y] << "\n";
        }
        correl.close();
}

		/****************/
		/*** Glauber ****/
		/****************
		{
		delta_t=1./N;
		T_TOT += delta_t;
		int x = rand_lx(generator);
		int y = rand_ly(generator);

		double e_flip = grille[modulo(x-1,TAILLE_X)][y]+grille[modulo(x+1,TAILLE_X)][y]+grille[x][modulo(y-1,TAILLE_Y)]+grille[x][modulo(y+1,TAILLE_Y)];
                e_flip *= 2*J*grille[x][y];
		e_flip += 2*grille[x][y]*champ_mag(x,y);
		if(rand_01(generator) < exp(-BETA*e_flip)) grille[x][y] *= -1;
		} //*/
