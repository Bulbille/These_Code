#ifndef INIT
#define INIT

#include <iostream>
#include <fstream>
#include <random>
#include <sys/file.h> //Files lock
#include <sstream>
#include <iomanip>
using namespace std;

const double T_C = 2./log(1+sqrt(2)),J=1;
extern double BETA,H,w,F,ttc;

// Paramètres de la grille

const int L_X =200;
const int L_Y =40;
const int NLiensMax = 2*(L_Y-2)*L_X-L_X;
//Au bout de 100 étapes de Monte Carlo environ, on a la décorrélation a`peu près maximale
//On va donc prendre les données uniquement tous les 100
extern int T_PHOTO,T_PHO;

extern std::default_random_engine generator;

#endif
