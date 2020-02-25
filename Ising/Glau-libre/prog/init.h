#ifndef INIT
#define INIT

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
const double J = 1.;

// Paramètres de la grille

const int TAILLE_X =600;
const int TAILLE_Y = 600;
const int N = TAILLE_X*TAILLE_Y;
//Au bout de 100 étapes de Monte Carlo environ, on a la décorrélation a`peu près maximale
//On va donc prendre les données uniquement tous les 100
const int T_PHOTO = 1.0e1*N;
const double NB_PHOTO = 100;
const double T_MAX = 1.*T_PHOTO*NB_PHOTO;

extern double BETA;
extern double ttc;

extern std::default_random_engine generator;
extern std::uniform_real_distribution<double> rand_01;
extern std::uniform_int_distribution<int> rand_lx;
extern std::uniform_int_distribution<int> rand_ly;

#endif
