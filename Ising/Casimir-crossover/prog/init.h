#ifndef INIT
#define INIT 

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <random>
#include <iomanip> // setprecision
#include <sstream> // stringstream


using namespace std;
const double T_C=2./log(1.+sqrt(2)), J=1, H0=0, w0=0;

// Paramètres de la grille

const int L_X = 50, L_Y = 10, N = L_X*L_Y;

const int T_PHO=1e3, T_EQ=1e0;

const double error = 1e-11;
const int Nlambda = 20, K0 = 0;

// Autres variables globales
extern int nb_liens; // compte le nombre de liens effectifs dans l'array liens
extern int nb_operations;
extern int nb_mesure;
extern double BETA;
extern double H;
extern double w;

extern default_random_engine generator;

#endif
