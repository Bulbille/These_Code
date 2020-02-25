#ifndef INIT
#define INIT 

//Output 
#include <iostream> //standard input-output
#include <fstream>  //files
#include <iomanip>  //setprecision
#include <stdlib.h> //atof()
#include <string>
#include <sstream> //stringstream

//Benchmark
#include <chrono>

//Mathematiques
#include <cmath>
#include <random>

//Structures
#include <map>
#include <unordered_map>


using namespace std;
const double T_C = 2./log(1.+sqrt(2));
const double J   = 1;
const double H0  = 0;

// Paramètres de la grille

const int TAILLE_X = 40;
const int TAILLE_Y = 40;
const int N = TAILLE_X*TAILLE_Y;

const int T_MAX    = 1e5;
const int T_PHOTO  = 1e2;
const int T_DONNEE = T_PHOTO*1e1;
const int T_CHAMP  = 1e0;

// Autres variables globales
const int tot_liens = TAILLE_X*(2*TAILLE_Y-1); // compte le nombre de liens effectifs dans l'array liens
extern int nb_operations;
extern int nb_mesure;
extern double BETA;
extern double H;

extern default_random_engine generator;

//template <int N>
//struct Iid {
//    static const int length = N;
//    int i[N];
//    int j[N];
//    double energie[N];
//};

struct Iid {
    int length;
    int *i;
    int *j;
    double *energie;
};


#endif
