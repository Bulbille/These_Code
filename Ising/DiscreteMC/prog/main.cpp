#include "init.h"
#include "math.h"
#include "dynamique.h"
#include <gsl/gsl_histogram.h>

void parametres(int argc, char* argv[], int grille[][TAILLE_Y]);
void ecriture(int grille[][TAILLE_Y], int temps, gsl_histogram histogramme,double magy[], double interface[]);

double BETA;
double ttc = 1;
string prefix = ".";
string suffix = "";
double H=H0;

std::default_random_engine generator;
std::uniform_real_distribution<double> rand_01(0.0,1.0);

//0 → L-1 : périodique
//1 → L-2 : fixe
std::uniform_int_distribution<int> rand_lx(0,TAILLE_X-1);
std::uniform_int_distribution<int> rand_ly(1,TAILLE_Y-2);

clock_t cputime = clock();


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

	/******* VARIABLES MICRO *********/
	gsl_histogram * histo = gsl_histogram_alloc (TAILLE_Y-1);
        gsl_histogram_set_ranges_uniform (histo, 0.1,TAILLE_Y-0.9);

	/**** VARIABLES INTERMÉDIAIRES *****/
	int nb_liens;
        int liens[TAILLE_X * TAILLE_Y * 5];
        double acceptance_rate[TAILLE_X * TAILLE_Y * 2];
	double delta_t=1./N;
	int photo_suiv=0; int nb_photos=1;double donnee_suiv=0;

	double magy[TAILLE_Y]; for(int y=0;y<TAILLE_Y;y++) magy[y]=0;

	/********************************************/
	/****** DYNAMIQUE DE KAWASAKI ***************/
	/*******************************************/

	for(double t = 0;t <= T_CHAMP; t+= delta_t){
	while(1){
		int x1 = rand_lx(generator);
		int y1 = rand_ly(generator);
		int x2 = rand_lx(generator);
		int y2 = rand_ly(generator);
		if(grille[x1][y1] != grille[x2][y2] && exchange(grille,x1,y1,x2,y2) <= rand_01(generator)){
			grille[x1][y1] *= -1;
			grille[x2][y2] *= -1;
			break;
		}
	}
	}

	for(double t = 0; t<=T_MAX ; t+= delta_t){
		cout << t << "\n";
		/******************************************************/
		/*********** CALCULS GRANDEURS THERMO *****************/
		/** moyenne_n = moyenne_(n-1) * (n-1) / n + a_n / n ***/
		/******************************************************/

		int modulo_photo = static_cast<int>(modulo_d(delta_t+t,T_PHOTO));
                if(photo_suiv <= static_cast<int>(t) && modulo_photo == 0) {
                        photo_suiv += T_PHOTO; nb_photos++;
			// L'interface est par construction à TAILLE_X/2
/*			double p_magy[TAILLE_Y]; for(int y=0;y<TAILLE_Y;y++) p_magy[y]=0;
			for(int y=0;y<TAILLE_Y;y++){
				for(int x=0;x<TAILLE_X;x++){
					p_magy[y] += 1.*grille[x][y] / TAILLE_X;
				}
				magy[y] = 1.*magy[y]*(nb_photos-1)/nb_photos+1.*p_magy[y]/nb_photos;
			}

			double interface[TAILLE_X]; double ancien_interface = 20;
			for(int x=0;x<TAILLE_X;x++){
				int interfacem=1;int interfacep=TAILLE_Y-2;
				bool plusmoins=true; int iterations=0; 
				while(interfacep-interfacem > 1){
					if(plusmoins){
						if(grille[x][interfacem-1] == grille[x][interfacem] 
							&& grille[x][interfacem] != grille[x][interfacem+1] 
							&& grille[x][interfacem+1] == grille[x][interfacem+2]){
							plusmoins = !plusmoins;
							iterations++;
						}
						else
							interfacem++;
					}
					else{
						if(grille[x][interfacep-1] == grille[x][interfacep] 
							&& grille[x][interfacep] != grille[x][interfacep+1] 
							&& grille[x][interfacep+1] == grille[x][interfacep+2]){
							plusmoins = !plusmoins;
							iterations++;
						}
						else
							interfacep--;
					}
					if(iterations > 1) break;
				}
				interface[x]=(interfacem+interfacep)/2.;
				cout << ancien_interface << " " << interfacep << " " << interfacem << "\n";
				if((abs(interfacep-ancien_interface) >= 7 && abs(interfacem-ancien_interface) >= 7) || abs(interfacep-interfacem) > 2){ interface[x] = -1;  continue;}
				ancien_interface = (interfacem+interfacep)/2.;
				gsl_histogram_increment(histo,(ancien_interface));
			}*/

			if(donnee_suiv <=static_cast<int>(t)){
				double interface[TAILLE_X]; double ancien_interface = 20;
				ecriture(grille,t,*histo,magy,interface);
                                donnee_suiv+=T_DONNEE;
                       } 
		}
		/****** MONTE CARLO STEP ****/
		while(1){
			int x1 = rand_lx(generator);
			int y1 = rand_ly(generator);
			double direction = rand_01(generator);
			int y2 = y1; int x2 = x1;
			if(direction < 0.25 && y1 > 1){ y2 = y2-1;}
			else if(direction < 0.5 && y1 < TAILLE_Y-1){ y2 = y2+1;}
			else if(direction < 0.75){ x2 = modulo(x1-1,TAILLE_X);}
			else x2 = modulo(x1+1,TAILLE_X);
	
			if(grille[x1][y1] != grille[x2][y2] && exchange(grille,x1,y1,x2,y2) >= rand_01(generator)){
				grille[x1][y1] *= -1;
				grille[x2][y2] *= -1;
				break;
			}
		}
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
void ecriture(int grille[][TAILLE_Y], int temps, gsl_histogram histogramme,double magy[], double interface[]){

        /*********** CREATION DOSSIER POUR RESULTATS ********/

        if(system(("mkdir -p " + prefix).c_str())) {;}
        string str = to_string(ttc);
	str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
        string algo = prefix + "/kaw_ttc_" + str + suffix + "-t" + to_string(0);
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

        thermo_res << ttc << " " << H << " " << temps << " " << float(clock()-cputime)/CLOCKS_PER_SEC << "\n";
	thermo_res.close();

        /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
        str = algo;
        ofstream result(str.c_str());
	str = algo+"-magy";
	ofstream fmagy(str.c_str());

	for(int y=0; y<TAILLE_Y;y++)
	{
		fmagy << y << " " << magy[y] << "\n";
		for(int x=0; x<TAILLE_X;x++){
			result << x << " " << y << " " << grille[x][y] << "\n";
		}//*/
	}
	fmagy.close();
	result.close();
	
        str = algo+"-histo";
        FILE* fhisto;
        fhisto = fopen(str.c_str(),"w");
        gsl_histogram_fprintf(fhisto,&histogramme, "%f","%f");
        fclose(fhisto);

	str = algo+"-interface";
	ofstream finter(str.c_str());
	for(int x=0;x<TAILLE_X;x++) {
		finter <<  interface[x] << "\n";
	}
	finter.close();
}

