#ifndef LW
#define LW
#include "init.h"
#include "dynamique.h"
int getCellX(int i);
int getCellY(int i);
int getCellFromXY(int i);
int getCellPPV(int i,int lien, int** voisins);

void def_voisins(int** voisins, int M);
int suite(int un,int L);

void etat_initial(int* grille,int M);
int couplage(int x, int y, std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens);

void liens(int M, int** ppv,  Iid &liens_ppv);
void def_couplage(std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens, Iid &liens_ppv, int nb_liens);
void generation(int* grille, int** ppv,Iid &liens_ppv,std::unordered_map<int,double> &liste_liens, std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens, int nb_liens, double* delta_t);

#endif

