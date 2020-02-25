#include "init.h"
#ifndef DYNAMIQUE
#define DYNAMIQUE

double exchange(int grille[][TAILLE_Y], int x1, int y1, int x2, int y2);
void exchange_nl(int grille[][TAILLE_Y], int x1, int y1, int x2,int y2);
double hamil_loc(int grille[][TAILLE_Y],int x,int y);
double champ_mag(int x, int y);
void generation(int grille[][TAILLE_Y]);
void generation(int grille[][TAILLE_Y], string& fichier);
void maj_liens(int grille[][TAILLE_Y], int liens[], double acceptance_rate[], double* delta_t, int* nb_liens);
int tirage_lien(double tableau[], int taille);


#endif
