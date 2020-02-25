#include "dynamique.h"
#include "math.h"
#include "init.h"
#include "base_func.h"

/************************************/
/*** EQUATIONS DE LA DYNAMIQUE ******/
/************************************/
double exchange(int* grille, int** ppv, int point1, int point2){
    int x1 = getCellX(point1);
    int y1 = getCellY(point1);
    int x2 = getCellX(point2);
    int y2 = getCellY(point2);
    double delta_h = 2.*grille[point1]*(champ_mag(x2,y2)-champ_mag(x1,y1));

    if( grille[point1] == grille[point2] || y1 <= 0 || y2 <= 0 || y1 >= TAILLE_Y-1 || y2 >= TAILLE_Y-1 ) 
        return 0;
    double h1 = -grille[point1],
           h2 = -grille[point2];

    //    for(int id=0;id<4;id++)
    //        cout << point1 << " " << id << " " <<   ppv[point1][id] << " " << "\n";

    for(int id=0; id<4; id++){
        h1 += grille[getCellPPV(point1,ppv[point1][id],ppv)]* (getCellPPV(point1,ppv[point1][id],ppv) >= 0);
        h2 += grille[getCellPPV(point2,ppv[point2][id],ppv)]* (getCellPPV(point2,ppv[point2][id],ppv) >= 0);
        //h2 += grille[ppv[point2][id]]* (ppv[point2][id] >= 0);
    }

    delta_h += 2.*J*grille[point1]*(h1-h2);
    return ((exp(-BETA*delta_h) <1 ) ? exp(-BETA*delta_h) : 1);
}

double champ_mag(int x,int y){
    (void)x ; // stoppe les warning
    double h = (y < TAILLE_Y/2) ? -H/2. : H/2. ;
    return h;
}

/*** Mise à jour des liens et tirage au sort ***/
/***********************************************/
void maj_liens(int* grille, int** ppv, Iid &liens_ppv, std::unordered_map<int,double> &liste_liens, std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens, double* delta_t, int tirage){

//    cout << "Màj tirage : " << tirage << "\n";
    for (int i = 0; i < 4; i++) {

        int voisin = ppv[liens_ppv.i[tirage]][i];
        int id_lien= couplage(liens_ppv.i[tirage],voisin,pointsLiens);
        if (id_lien >= 0){
            *delta_t -= liens_ppv.energie[id_lien];
            liens_ppv.energie[id_lien]= exchange(grille,ppv, liens_ppv.i[id_lien],liens_ppv.j[id_lien]);
//            cout << "Màj : " <<id_lien << " " << (liens_ppv.energie[id_lien] > 1e-8)  << "\n";

            if(liens_ppv.energie[id_lien] > 1e-15)
                liste_liens[id_lien] = liens_ppv.energie[id_lien];
            else
                liste_liens.erase(id_lien);

            *delta_t += liens_ppv.energie[id_lien];
        }
        voisin = ppv[liens_ppv.j[tirage]][i];
        id_lien= couplage(liens_ppv.j[tirage],voisin,pointsLiens);
        if (id_lien >= 0){
            *delta_t -= liens_ppv.energie[id_lien];
            liens_ppv.energie[id_lien]= exchange(grille,ppv, liens_ppv.i[id_lien],liens_ppv.j[id_lien]);
//            cout << "Màj : " <<id_lien << " " << (liens_ppv.energie[id_lien] > 1e-8)  << "\n";

            if(liens_ppv.energie[id_lien] > 1e-15)
                liste_liens[id_lien] = liens_ppv.energie[id_lien];
            else
                liste_liens.erase(id_lien);

            *delta_t += liens_ppv.energie[id_lien];
        }
    }
}

int tirage_lien(std::unordered_map<int,double> &liste_liens){
        auto debut = std::chrono::high_resolution_clock::now();
    //Mise dans un tableau 1D de la map
    int taille = liste_liens.size();
    double* tab_tmp = new double[taille];
    int* index   = new int[taille];
    index[0]     = liste_liens.begin()->first;
    tab_tmp[0]   = liste_liens.begin()->second;

    int i=0;
    for ( auto it = std::next(liste_liens.begin()); it != liste_liens.end(); ++it ){
        i++;
        tab_tmp[i]  = tab_tmp[i-1] + it->second;
        index[i]    = it->first;
    }

                    auto fin = std::chrono::high_resolution_clock::now();
                        auto duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
                            std::cout << "Fonction : " << duree <<  " " << taille << std::endl;
                            
//    cout << "++++ \n";
//    for(int j=0; j<=i ; j++)
//        cout << index[j] << "\n";


    //Tirage au sort et dichotomie
    std::uniform_real_distribution<double> rand_rate(0,tab_tmp[taille-1]);
    double tirage = rand_rate(generator);
    int id = 0; int ifin = taille-1;
    if(tirage < tab_tmp[0]) return index[0];
    while(ifin-id > 1){
        if(tab_tmp[ifin] > tirage && tab_tmp[id] > tirage) break;
        int m = static_cast<int>((id+ifin)/2);
        if(tirage > tab_tmp[m]) id = m;
        else  ifin = m;
    }
    return index[ifin];
}
