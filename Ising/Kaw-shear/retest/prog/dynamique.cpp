#include "dynamique.h"
#include "math.h"
#include "init.h"

/*********** Calcul de Delta_E *******/
double exchange(int grille[][L_Y], int x1, int y1, int x2, int y2){
    if( grille[x1][y1] == grille[x2][y2] ) return 0;

    double h_1 = grille[modulo(x1+1,L_X)][y1]+grille[modulo(x1-1,L_X)][y1]+grille[x1][y1+1]+grille[x1][y1-1]-grille[x2][y2];
    double h_2 = grille[modulo(x2+1,L_X)][y2]+grille[modulo(x2-1,L_X)][y2]+grille[x2][y2+1]+grille[x2][y2-1]-grille[x1][y1];
    double delta_h = 2.*J*grille[x1][y1]*(h_1-h_2);

    double mag = 2.*grille[x1][y1]*(champ_mag(x2,y2)-champ_mag(x1,y1));
//    double shear_w = J*w*grille[x1][y1]*(y1 == y2)*(y1 - (L_Y+1)/2.);
//    double shear_f = J*F*grille[x1][y1]*(y1 == y2)*((y1==L_Y-2)-(y1==2));

    double proba = exp(-BETA*(delta_h+mag));
    return ( proba < 1) ? proba : 1 ;
}

double champ_mag(int x1,int y1){
    return (y1 < L_Y/2) ? -H/2. : H/2. ;
}

/******* Mise à jour des liens *********/
void maj_liens(int grille[][L_Y], int* liens, double* acceptance_rate, double* delta_t){
    *delta_t = 0;
    int x1,x2,y1,y2;
    for(int nblien = 0 ; nblien < NLiensMax ; nblien++){
        x1 = liens[nblien*4+0];
        y1 = liens[nblien*4+1];
        x2 = liens[nblien*4+2];
        y2 = liens[nblien*4+3];
        acceptance_rate[nblien] = exchange(grille,x1,y1,x2,y2); 
        *delta_t += acceptance_rate[nblien]; 
    }
    *delta_t = 1/ *delta_t;
}

/****** GÉNÉRATION DU RÉSEAU *******/
void generation(int grille[][L_Y],int *liens){
    /* Génération de la topologie */
    int nb_liens = 0;
    for(int y=1;y<L_Y-1;y++){
        for(int x=0;x<L_X;x++){
            //Vérification à droite
            liens[nb_liens *4] = x;
            liens[nb_liens *4+1] = y;
            liens[nb_liens *4+2] = modulo(x+1,L_X);
            liens[nb_liens *4+3] = y;
            nb_liens+= 1;
            //Vérification en haut
            if(y<L_Y-2){
                liens[nb_liens *4] = x;
                liens[nb_liens *4+1] = y;
                liens[nb_liens *4+2] = x;
                liens[nb_liens *4+3] = y+1;
                nb_liens+= 1;
            }
        }
    }
    /* Configuration initiale */
    for(int j=0;j<L_Y;j++){
        for(int i=0;i<L_X;i++){
            if(L_Y % 2 ==0){
                if(j>=L_Y/2)
                    grille[i][j] = -1;
                else
                    grille[i][j] = 1;
            }
            else{
                if(j>L_Y/2)
                    grille[i][j] = -1;
                else if(j == L_Y/2 && i < L_X/2)
                    grille[i][j] = -1;
                else
                    grille[i][j] = 1;
            }		
        }
    }
}

/********* Tirage d'un lien au hasard ************/
int tirage_lien(double* tableau){
    double tableau2[NLiensMax];
    tableau2[0] = tableau[0];
    for(int i=1;i<NLiensMax;i++) tableau2[i] = tableau[i] + tableau2[i-1];

    std::uniform_real_distribution<double> rand_rate(0,tableau2[NLiensMax-1]);     
    double tirage = rand_rate(generator);
    if(tirage < tableau2[0]) return 0;

    int id = 0, ifin = NLiensMax-1;

    while(ifin-id > 1){
        if(tableau2[ifin] > tirage && tableau2[id] > tirage) break;
        int m = static_cast<int>((id+ifin)/2);
        if(tirage > tableau2[m]) id = m;
        else  ifin = m;
    }
    return ifin;
}
