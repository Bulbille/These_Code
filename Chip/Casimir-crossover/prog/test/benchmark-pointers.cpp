#include <iostream> //standard Input-Output
#include <chrono>

using namespace std;
const int TAILLE = 30;
//Propriétés de l'Hamiltonien
/***********************************/
int fonction(int* tableau,int taille){
    int somme = 0;
    for(int i=0;i<taille;i++) somme+=tableau[i];
    return somme;
}

int fonction2(int tableau[],int taille){
    int somme = 0;
    for(int i=0;i<taille;i++) somme+=tableau[i];
    return somme;
}

int main(int argc,char* argv[]){
    int tab[TAILLE];
    int* tab2 = new int[TAILLE];
    for(int i=0;i<TAILLE;i++){ tab[i] = i; tab2[i] = i;}


    auto debut = std::chrono::high_resolution_clock::now();
    for(int t=0;t<5e8;t++) fonction(&tab,TAILLE);
    auto fin = std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
    std::cout << "F1 Tab 1 : " << duree/1e9 <<  std::endl;

     debut = std::chrono::high_resolution_clock::now();
    for(int t=0;t<5e8;t++) fonction(tab2,TAILLE);
     fin = std::chrono::high_resolution_clock::now();
     duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
    std::cout << "F1 Tab 2 : " << duree/1e9 <<  std::endl;

     debut = std::chrono::high_resolution_clock::now();
    for(int t=0;t<5e8;t++) fonction2(tab,TAILLE);
     fin = std::chrono::high_resolution_clock::now();
     duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
    std::cout << "F2 Tab 1 : " << duree/1e9 <<  std::endl;

     debut = std::chrono::high_resolution_clock::now();
    for(int t=0;t<5e8;t++) fonction2(tab2,TAILLE);
     fin = std::chrono::high_resolution_clock::now();
     duree = std::chrono::duration_cast<std::chrono::nanoseconds>(fin-debut).count();
    std::cout << "F2 Tab 2 : " << duree/1e9 <<  std::endl;
return 0;
}
