#include "base_func.h"
#include "math.h"
#include "init.h"

int getCellX(int i){
    return i%TAILLE_X;
}
int getCellY(int i){
    return static_cast<int>(i/TAILLE_X);
}
int getCellFromXY(int x,int y){
    return (y < 0 || y >= TAILLE_Y)? -1 : y*TAILLE_X+modulo(x,TAILLE_X);
}
int getCellPPV(int i,int lien, int** voisins){
    return (i == voisins[lien][0]) ? voisins[lien][1] : voisins[lien][0];
}

/************ Génération de la grille *********/
void etat_initial(int* grille,int M){
    //*Génération magnétisation strictement égale à 0
    for(int i=0 ; i<M ; i++){
        int x = getCellX(i);
        int y = getCellY(i);
        if(TAILLE_Y % 2 ==0){
            if(y>=TAILLE_Y/2)
                grille[i] = -1;
            else
                grille[i] = 1;
        }
        else{
            if(y>TAILLE_Y/2)
                grille[i] = -1;
            else if(y == TAILLE_Y/2 && x < TAILLE_X/2)
                grille[i] = -1;
            else
                grille[i] = 1;
        }  
     //   (y == 5 && ( x==2 || x == 3) ) ? grille[i] = 1 : grille[i] = -1;

    }
}

/******* Listage des plus proches voisins de chaque spin *******/
void def_voisins(int** voisins, int M){
    int a=0;
    for(int i=0; i<M; i++){
        int x = getCellX(i);
        int y = getCellY(i);
        int tmp = 0;
        int test_cell = getCellFromXY(x+1, y);
        if(test_cell >= 0){
            voisins[i][tmp] = test_cell;
            tmp++;
        }
        test_cell = getCellFromXY(x-1, y);
        if(test_cell >= 0){
            voisins[i][tmp] = test_cell;
            tmp++;
        }
        test_cell = getCellFromXY(x, y+1);
        if(test_cell >= 0){
            voisins[i][tmp] = test_cell;
            tmp++;
        }
        test_cell = getCellFromXY(x, y-1);
        if(test_cell >= 0){
            voisins[i][tmp] = test_cell;
            tmp++;
        }
        a+= tmp;
        while(tmp<4){
            voisins[i][tmp] = -1;
            tmp++;
        }
    }
}

// Pour éviter de compter deux fois chaque lien, on compte un point sur deux (différencier cas LX pair et LX impair) via suite()
// Utilisé dans la fonction liens()
int suite(int un, int L){
    int tmp;
    //Si on a un nombre impair, le décalage se fait tout seul
    if(modulo(L,2)!=0)            return un+2;
    if      (un == 0)             return un+2;

    tmp = un+2;
    if      (modulo(tmp-1,L) == 0 ) tmp -= 1;
    else if (modulo(tmp,L) == 0)    tmp += 1;

    return tmp;
}

int couplage(int x, int y, std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens){
   return pointsLiens[min(x,y)][max(x,y)]; 
}


void liens(int M, int** ppv,Iid &liens_ppv){
    int id_lien=0;
    for(int i=0; i<M ; i = suite(i,TAILLE_X)){
        for(int tmp = 0; tmp <4 ; tmp++){
            //étude des ppv
            if(ppv[i][tmp] > -1){
                liens_ppv.i[id_lien] = i;
                liens_ppv.j[id_lien] = ppv[i][tmp];
                id_lien++;
            }   
        }   
    }   
}

void def_couplage(std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens, Iid &liens_ppv,int nb_liens){
    for(int id=0 ; id<nb_liens ; id++){
        int x = min(liens_ppv.i[id],liens_ppv.j[id]);
        int y = max(liens_ppv.i[id],liens_ppv.j[id]);
        pointsLiens[x][y] = id; 
    }   
}

void generation(int* grille, int** ppv, Iid &liens_ppv,std::unordered_map<int,double> &liste_liens,std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens, int nb_liens, double* delta_t){
    /************ Génération de la grille *********/
    etat_initial(grille,TAILLE_X*TAILLE_Y);
    def_voisins(ppv,TAILLE_X*TAILLE_Y);
    liens(TAILLE_X*TAILLE_Y, ppv, liens_ppv);
    def_couplage(pointsLiens,liens_ppv,nb_liens);

    //Calcul des dt de chaque lien
    *delta_t = 0;
    for(int id=0; id<nb_liens; id++){
        liens_ppv.energie[id] = exchange(grille,ppv,liens_ppv.i[id],liens_ppv.j[id]);
//        cout << id << " | " <<  liens_ppv.i[id] << " " << liens_ppv.j[id] << " " << liens_ppv.energie[id] << "\n";

        //référencement des liens
        if(liens_ppv.energie[id] > 1e-8)
            liste_liens[id] = liens_ppv.energie[id];
        else
            liste_liens.erase(id);

        *delta_t+=liens_ppv.energie[id];
    }
}

