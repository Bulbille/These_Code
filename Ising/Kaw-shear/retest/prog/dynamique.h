#include "init.h"
#ifndef DYNAMIQUE
#define DYNAMIQUE

double exchange(int grille[][L_Y], int x1, int y1, int x2, int y2);
double champ_mag(int x1,int y1);
void generation(int grille[][L_Y],int* liens);
void maj_liens(int grille[][L_Y], int* liens, double* acceptance_rate, double* delta_t);
int tirage_lien(double* tableau);


#endif
