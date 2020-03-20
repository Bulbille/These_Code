#ifndef INIT
#define INIT

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <random>
#include <string>
#include <iomanip>
using namespace std;

// Paramètres physiques du système
/* Beta = 1/k_B T */
/* J : constante de couplage du modèle d'Ising* */
/* Intensité du champ magnétique */
/* */
const double T_C = 2./log(1+sqrt(2));
const int Tn = 8;
const double Tmin = 0, Tmax = 1;
const double J = 1.;

// Paramètres de la grille
const int bitLX = 6, TAILLE_X = 2<<bitLX;
const int TAILLE_Y = 51;
const int N = TAILLE_X*TAILLE_Y;
//Au bout de 100 étapes de Monte Carlo environ, on a la décorrélation a`peu près maximale
//On va donc prendre les données uniquement tous les 100
const double T_EQ =5e5, T_MAX = 5e5;

extern double BETA;
extern double ttc;

#endif
