#ifndef INIT
#define INIT 

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <stdlib.h> //atof()
#include <sstream>
#include <gsl/gsl_histogram.h>
#include <random>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <chrono>
#include <sys/file.h> //Files lock


using namespace std;
const double T_C = 2./log(1.+sqrt(2));
const double J = 1;
const double H0 = 0;
const double w0 = 0;

// Paramètres de la grille

const int L_X = 200;
const int L_Y = 20;
const int N = L_X*L_Y;

const int T_PHOTO = 1e2;
const int T_DONNEE = 1e0;
const int T_CHAMP = 1e3;

// Autres variables globales
extern int nb_liens; // compte le nombre de liens effectifs dans l'array liens
extern int nb_operations;
extern int nb_mesure;
extern double BETA;
extern double H;
extern double w;

extern default_random_engine generator;

#endif
