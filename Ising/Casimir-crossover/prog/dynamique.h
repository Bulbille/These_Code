#ifndef DYNAMIQUE
#define DYNAMIQUE
#include "init.h"

double exchange(int grille[][L_Y], int x1, int y1, int x2, int y2, int K0);
double hamil_loc(int grille[][L_Y],int x,int y);
double champ_mag(int x1, int y1);
double cisaillement(int grille[][L_Y], int x1,int y1, int x2, int y2);
void generation(int grille[][L_Y]);
void maj_liens(int grille[][L_Y], int liens[], double acceptance_rate[], double* delta_t, int* nb_liens,int K0);
int tirage_lien(double tableau[], int taille);

#endif
