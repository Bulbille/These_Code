#include "init.h"
#include "math.h"
#include "dynamique.h"
#include "base_func.h"

void parametres(int argc, char* argv[]),
     ecriture(int* grille, long int temps);

double BETA = 1/T_C,
       ttc  = 1,
       H    = H0;
string prefix = ".",
       suffix = "";

std::default_random_engine generator;
clock_t cputime = clock();


/****************** Main ******************/
int main(int argc, char *argv[]){

    int grille[N];
    std::unordered_map<int, std::unordered_map<int,int>> pointsLiens;
    std::unordered_map<int,double> liste_liens;
    int **ppv = new int*[N];
    for(int i=0; i<N ; i++)
        ppv[i] = new int[4];
    Iid liens_ppv = FactoryIid (tot_liens);
    double delta_t = 0,
           dt      = 0;
    generator.seed(time(NULL));


    long int photo_suiv  = 0,
             nb_photos   = 0,
             donnee_suiv = 0;

/********** OPTIONS DEPUIS LE BASH *********/
/* Le programme peut accepter deux options : */
/* --beta xxx */
/* --reprendre fichier */
/* La première option permet de faire la simulation avec une température donnée */
/* La deuxième permet de reprendre la simulation à partir d'un fichier */
/**/	parametres(argc,argv); /**/
/**/	generation(grille,ppv,liens_ppv,liste_liens,pointsLiens,tot_liens,&delta_t); dt=1/delta_t;
/**/    ttc = static_cast<double>(static_cast<int>(ttc*1000))/1000;
/******************************************/
    /********************************************/
    /****** DYNAMIQUE DE KAWASAKI ***************/
    /*******************************************/
    int ndt = 0;
//    cout << "++++ T=0 "<< liste_liens.size() << "\n";;
//    for (std::map<int,double>::iterator it=liste_liens.begin(); it!=liste_liens.end(); ++it){
//        cout << it->first << " " << it->second << "\n";
//    }

    for(double t = 0; t<=T_MAX ; t+= dt){		
        /******************************************************/
        /*********** CALCULS GRANDEURS THERMO *****************/
        /** moyenne_n = moyenne_(n-1) * (n-1) / n + a_n / n ***/
        /******************************************************/

        int modulo_photo = static_cast<int>(modulo_d(dt+t,T_PHOTO));
        if(photo_suiv <= static_cast<int>(t) && modulo_photo == 0) {
            photo_suiv += T_PHOTO; nb_photos++;
            if(donnee_suiv <=static_cast<int>(t)){
                cout << "T = " << t << "\n";
                ecriture(grille,t);
                donnee_suiv+=T_DONNEE;
            } 
        }
        /****** MONTE CARLO STEP ****/
    auto debut = std::chrono::high_resolution_clock::now();
        int tirage = tirage_lien(liste_liens);
    auto fin = std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
    std::cout << duree <<  std::endl;
        //        cout << tirage << " (" << liens_ppv.i[tirage] << "," << grille[liens_ppv.i[tirage]] << ") ("<< liens_ppv.j[tirage] << "," << grille[liens_ppv.j[tirage]] << ") \n";

        grille[liens_ppv.i[tirage]] *= -1;
        grille[liens_ppv.j[tirage]] *= -1;

        maj_liens(grille, ppv, liens_ppv,liste_liens, pointsLiens, &delta_t,tirage);
        
        dt = 1/delta_t;
        ndt++;

//        for (std::map<int,double>::iterator it=liste_liens.begin(); it!=liste_liens.end(); ++it){
//            cout << it->first << " " << it->second << "\n";
//        }
    }
    ecriture(grille,T_MAX);
    for ( auto it = std::next(liste_liens.begin()); it != liste_liens.end(); ++it ){
        cout << it->first << " " << it->second << "\n";
    }
/*************** FIN  **********************/
    for (int i = 0; i < N; i++) delete[] ppv[i];
        delete[] ppv;
return 0;
}

/******* LECTURE DES PARAMÈTRES DU BASH **********/
void parametres(int argc, char* argv[]){
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
			else if(arg == "--prefix"){
				prefix = argv[i+1];
			}
			else if(arg == "--suffix"){
				suffix = argv[i+1];
			}
		}
	}
}

/************* ÉCRITURE DANS FICHIER **********/
void ecriture(int* grille, long int temps){

    /*********** CREATION DOSSIER POUR RESULTATS ********/
    if(system(("mkdir -p " + prefix).c_str())) {;}
    string algo = prefix + "/";
    {
        stringstream stream;
        stream << fixed << setprecision(2) << ttc ;
        algo += "kaw_T-" + stream.str();
    }
    if(fabs(H) > 1e-6){
        stringstream stream;
        stream << fixed << setprecision(2) << H;
        algo += "-h-" + stream.str();
    }
    /**** ECRITURE DANS FICHIER DERNIER ETAT DU RESEAU **/
    string str = algo;
    ofstream result(str.c_str());

    for(int i=0; i<TAILLE_X*TAILLE_Y;i++){
        result << getCellX(i) << " " << getCellY(i) << " " << grille[i] << "\n";
    }
    result.close();
}

