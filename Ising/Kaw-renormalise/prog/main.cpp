#include "init.h"
#include "math.h"
#include "dynamique.h"

#include "time.h"
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
	int n_correl[TAILLE_X][TAILLE_Y];
	/********/


	/******* VARIABLES THERMO *********/
        double magnetisation, e_int, mag2, e2;
                magnetisation=0; e_int=0;mag2=0;e2=0;
        double mag_x[TAILLE_X];
        for(int i=0;i<TAILLE_X;i++)
                mag_x[i]=0;

        double e_x[TAILLE_X];
        double e_y[TAILLE_X];
        for(int x=0;x<TAILLE_X;x++){
                e_x[x]=0; e_y[x]=0;
        }


	/**** VARIABLES INTERMÉDIAIRES *****/
	int nb_liens;
        int liens[TAILLE_X * TAILLE_Y * 5];
        double acceptance_rate[TAILLE_X * TAILLE_Y * 2];
	double T_TOT = 0;
	double delta_t=0;
	int nb_photos=-1;

	double H_prime = H;
	H=H0;	
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
		 //*/

        
		/*************************************/
		/*** RENORMALISATION DE LA GRILLE ****/
		/*************************************/
		int grille_renormalisee[TAILLE_X][TAILLE_Y];
		for(int x=0;x<TAILLE_X;x++){
			for(int y=0;y<TAILLE_Y;y++){
				int a = grille[modulo(x-1,TAILLE_X)][y] + grille[modulo(x+1,TAILLE_X)][y] + grille[x][modulo(y-1,TAILLE_Y)] + grille[x][modulo(y+1,TAILLE_Y)] ;
				grille_renormalisee[x][y] = (a > 0 ) - (a < 0);
			}
		}
		
		/****** CALCULS GRANDEURS THERMO **************/

		double mag_conf = 0;
		double e_conf= 0;

		for(int y=0;y<TAILLE_Y;y++){
		for(int x=0;x<TAILLE_X;x++){
			mag_conf += grille_renormalisee[x][y];
			e_conf += hamil_loc(grille_renormalisee,x,y)/2;
			mag_x[x] += 1.*grille[x][y]*delta_t;
			if(x < TAILLE_X-1 && x > 0){/////////////////7 on s'en fout complètement.... faudrait changer à e_y[x], e_x[x] 
				e_y[x] -= 1.*grille_renormalisee[x][y]*grille_renormalisee[x][y+1]*delta_t;
				e_x[x] -= 1.*grille_renormalisee[x][y]*grille_renormalisee[modulo(x+1, TAILLE_X)][y]*delta_t; 
			}
		}
		}
		magnetisation += abs(mag_conf)*delta_t/N;
		e_int += e_conf/N*delta_t;
		mag2 += mag_conf/N*mag_conf/N*delta_t;
		e2 += e_conf/N*e_conf/N*delta_t;
	

                int modulo_photo = static_cast<int>(modulo_d(t,T_PHOTO));
                if(modulo_photo == 0 && static_cast<int>(t/T_PHOTO) != nb_photos){


			/*** Calcul grandeurs thermo ***/

                        double p_mag = magnetisation/T_TOT;
                        double p_e = e_int/T_TOT;
                        double p_e2 = e2/T_TOT;
                        double p_m2 = mag2/T_TOT;
                        double p_magx[TAILLE_X], p_ex[TAILLE_Y], p_ey[TAILLE_Y];
                        for(int x=0; x<TAILLE_X;x++){
                                p_ex[x] = e_x[x] / (TAILLE_X*T_TOT);
                                p_ey[x] = e_y[x] / (TAILLE_X*T_TOT);
                        }
			for(int x=0;x<TAILLE_X;x++){ 
				p_magx[x] = mag_x[x] / (TAILLE_Y*T_TOT);
			}
				
                        double chi = BETA*(p_m2 - pow(p_mag,2))*N;
                        double c_v = pow(BETA,2)*(p_e2 - pow(p_e,2))*N;

			double p_correl_y[TAILLE_X][TAILLE_Y];
			for(int j=0;j<TAILLE_Y;j++){
				for(int i=0;i<TAILLE_X;i++){
					p_correl_y[i][j] = 0;
					n_correl[i][j] = 0;
				}
			}

			for(int x=0;x<TAILLE_X;x++){
				for(int y=0;y<TAILLE_Y;y++){
					for(int y2=y;y2<TAILLE_Y+y;y2++){
						p_correl_y[x][y2-y] += 1.*grille_renormalisee[x][y]*grille_renormalisee[x][modulo(y2,TAILLE_Y)];
						n_correl[x][y2-y]++;
					}
				}
			}
			for(int j=0;j<TAILLE_Y;j++){
				for(int i=0;i<TAILLE_X;i++){
					p_correl_y[i][j] = 1.*p_correl_y[i][j]/n_correl[i][j] - p_magx[i]*p_magx[i];
				}
			}
			for(int x=0;x<TAILLE_X;x++){
				for(int y=0;y<TAILLE_Y;y++){
					correl_y[x][y] += p_correl_y[x][y];
			//		cout << x << " " << y << " " << correl_y[x][y] << " " << p_correl_y[x][y] << "\n";
				}
			}
			
                        ecriture(grille_renormalisee,p_magx, p_ex, p_ey,p_mag, p_e, chi, c_v,static_cast<int>(0), correl_y);
                        nb_photos++;
                }

	}

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
void ecriture(int grille[][TAILLE_Y], double mag_x[], double e_x[], double e_y[], double magnetisation, double e_int, double chi, double c_v,int temps, double correl_y[][TAILLE_Y]){

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
