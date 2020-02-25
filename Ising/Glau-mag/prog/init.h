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
using namespace std;

// Paramètres physiques du système
/* Beta = 1/k_B T */
/* J : constante de couplage du modèle d'Ising* */
/* Intensité du champ magnétique */
/* */
const double T_C = 2./log(1+sqrt(2));
const double J = 1.;
const double H = 0;
const double w0 = 0;
const double F0 = 0;

// Paramètres de la grille

const int TAILLE_X =20;
const int TAILLE_Y = 20;
const int N = TAILLE_X*TAILLE_Y;
const int L_FAISCEAU = 0;
//Au bout de 100 étapes de Monte Carlo environ, on a la décorrélation a`peu près maximale
//On va donc prendre les données uniquement tous les 100
const int T_PHOTO = 1.0e1;
const double NB_PHOTO = 1e3;
const double T_MAX = 1.*T_PHOTO*NB_PHOTO;
const double T_THERM = 1e2;

extern double BETA;
extern double energie_totale ;
extern double w ;
extern double F ;
extern double ttc;

extern std::default_random_engine generator;
extern std::uniform_real_distribution<double> rand_01;
extern std::uniform_int_distribution<int> rand_lx;
extern std::uniform_int_distribution<int> rand_ly;

#endif
