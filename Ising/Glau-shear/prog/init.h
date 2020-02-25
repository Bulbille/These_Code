#ifndef INIT
#define INIT

#include <iostream>
#include <fstream>
#include <random>
#include <sys/file.h> //Files lock
#include <sstream>
#include <iomanip>
using namespace std;

#define ALGO 0
// 0 pour Kawasasaki, 1 pour Glauber

const double T_C = 2./log(1+sqrt(2)),J=1;
extern double BETA,H,w,F,ttc;

// Paramètres de la grille

const int TAILLE_X =200;
const int TAILLE_Y =40;
const int N = TAILLE_X*TAILLE_Y;
const int NLiensMax = (2*TAILLE_Y-1)*TAILLE_X;
//Au bout de 100 étapes de Monte Carlo environ, on a la décorrélation a`peu près maximale
//On va donc prendre les données uniquement tous les 100
extern long int T_PHOTO,T_MAX;

extern std::default_random_engine generator;
extern std::uniform_real_distribution<double> rand_01;
extern std::uniform_int_distribution<int> rand_lx;
extern std::uniform_int_distribution<int> rand_ly;

#endif
