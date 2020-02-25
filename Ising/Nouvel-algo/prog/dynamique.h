#ifndef DYNAMIQUE
#define DYNAMIQUE
#include "init.h"
#include "base_func.h"

double exchange(int* grille, int** ppv,int point1, int point2);
double champ_mag(int x,int y);

//int tirage_lien(Iid &tab);
int tirage_lien(std::unordered_map<int,double> &liste_liens);
void maj_liens(int* grille, int** ppv, Iid &liens_ppv,std::unordered_map<int,double> &liste_liens, std::unordered_map<int,std::unordered_map<int,int>> &pointsLiens, double* delta_t, int tirage);

#endif
