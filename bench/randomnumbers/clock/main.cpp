//Basic stuff
#include <random> 
#include <iostream> //cout
#include <chrono>

using namespace std;
//Générateur nombres aléatoires
default_random_engine defgen;
minstd_rand defminst;
mt19937 defmt;
mt19937_64 defmt64;


/***********************************/
/**** Définitions des fonctions ****/

int main(int argc,char* argv[]){
    std::uniform_int_distribution<int> rand_lx(0,1);

    int tirage;

    clock_t cputime = clock();
    for(int t=0;t<2e8;t++)
        tirage  = rand_lx(defgen);
    cout << "Default_random : " << (clock()-cputime)/CLOCKS_PER_SEC << "\n";
    cputime = clock();
    for(int t=0;t<2e8;t++)
        tirage  = rand_lx(defminst);
    cout << "minstd : "<< (clock()-cputime)/CLOCKS_PER_SEC << "\n";
    for(int t=0;t<2e8;t++)
        tirage  = rand_lx(defmt);
    cout << "mt199 : " << (clock()-cputime)/CLOCKS_PER_SEC << "\n";
    for(int t=0;t<2e8;t++)
        tirage  = rand_lx(defmt64);
    cout << "mt199 64 : "<< (clock()-cputime)/CLOCKS_PER_SEC << "\n";
return 0;
}
