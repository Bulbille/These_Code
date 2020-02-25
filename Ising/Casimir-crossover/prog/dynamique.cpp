#include "dynamique.h"
#include "math.h"
#include "init.h"

/***********************************************/
/*** Mise à jour des liens et tirage au sort ***/
/***********************************************/

double exchange(int grille[][L_Y], int x1, int y1, int x2, int y2, int K0){
    if ( grille[x1][y1] == grille[x2][y2]) return 0;
    double bc = ( y1 <= 0 or y2 <= 0 or y1 >= L_Y-1 or y2 >= L_Y-1 );
        
    double h_1 = grille[modulo(x1+1,L_X)][y1]+grille[modulo(x1-1,L_X)][y1]+grille[x1][y1+1]+grille[x1][y1-1]-grille[x2][y2];
    double h_2 = grille[modulo(x2+1,L_X)][y2]+grille[modulo(x2-1,L_X)][y2]+grille[x2][y2+1]+grille[x2][y2-1]-grille[x1][y1];
    double delta_h = 2.*J*grille[x1][y1]*(h_1-h_2);

    int mag = 2.*grille[x1][y1]*(champ_mag(x2,y2)-champ_mag(x1,y1));
    int shear = cisaillement(grille,x1,y1,x2,y2);

    return ((exp(-BETA*(bc+delta_h)) <1 ) ? exp(-BETA*(bc+delta_h+mag+shear)) : 1);
}

double champ_mag(int x1,int y1){
    double h1 = (y1 < L_Y/2) ? -H/2. : H/2. ;
    return h1;
}

double cisaillement(int grille[][L_Y], int x1,int y1, int x2, int y2){
    //return (grille[x1][y1] == 1 && y1 == y2)*w*(y1 - (L_Y+1)/2.);
    return J*w*grille[x1][y1];
}

/********** ÉNERGIE LOCALE *******/
double hamil_loc(int grille[][L_Y],int x,int y){
    x = modulo(x,L_X);
    y = modulo(y,L_Y);
    double energie = grille[modulo(x+1,L_X)][y]+grille[modulo(x-1,L_X)][y]+grille[x][modulo(y+1,L_Y)]+grille[x][modulo(y-1,L_Y)];
    energie += champ_mag(x,y);
    return -J*energie*grille[x][y];
}

/******* CONDITIONS PÉRIODIQUES EN X ****/
/******* CONDITIONS PÉRIODIQUES EN Y ****/

void maj_liens(int grille[][L_Y], int liens[], double acceptance_rate[], double* delta_t, int* nb_liens, int K0){

    *nb_liens = 0;
    *delta_t = 0;
    for(int x=0;x<L_X;x++){
        for(int y=1;y<L_Y-1;y++){
            int x2 = modulo(x+1,L_X);
            //Vérification à droite
            if(grille[x][y] != grille[x2][y]){
                liens[*nb_liens *4] = x;
                liens[*nb_liens *4+1] = y;
                liens[*nb_liens *4+2] = x2;
                liens[*nb_liens *4+3] = y;
                acceptance_rate[*nb_liens] = exchange(grille,x,y,x2,y);
                *delta_t += acceptance_rate[*nb_liens];
                *nb_liens += 1;
            }
            //Vérification en haut
            int y2 = modulo(y+1,L_Y);
            if(y2 < L_Y-1 && grille[x][y] != grille[x][y2]){
                liens[*nb_liens *4] = x;
                liens[*nb_liens *4+1] = y;
                liens[*nb_liens *4+2] = x;
                liens[*nb_liens *4+3] = y2;
                acceptance_rate[*nb_liens] = exchange(grille,x,y,x,y2);
                *delta_t += acceptance_rate[*nb_liens];
                *nb_liens += 1;
            }
        }
    }
    *delta_t = 1/ *delta_t;
}

int tirage_lien(double tableau[], int taille){
    double tableau2[2*N];
    tableau2[0] = tableau[0];

    for(int i=1;i<taille;i++) tableau2[i] = 0;
    for(int i=1;i<taille;i++) tableau2[i] = tableau[i] + tableau2[i-1];

    std::uniform_real_distribution<double> rand_rate(0,tableau2[taille-1]);
    double tirage = rand_rate(generator);
    int id = 0; int ifin = taille-1;

    if(tirage < tableau2[0]) return 0;
    while(ifin-id > 1){
        if(tableau2[ifin] > tirage && tableau2[id] > tirage) break;
        int m = static_cast<int>((id+ifin)/2);
        if(tirage > tableau2[m]) id = m;
        else  ifin = m;
    }

    return ifin;
}




/***********************************/
/****** GÉNÉRATION DU RÉSEAU *******/
/***********************************/

void generation(int grille[][L_Y]){
    for(int i=0;i<L_X;i++){
        for(int j=0;j<L_Y;j++){
            /*
            if(L_Y % 2 ==0){
                if(j>=L_Y/2)
                    grille[i][j] = +1;
                else
                    grille[i][j] = -1;
            }
            else{
                if(j>L_Y/2)
                    grille[i][j] = +1;
                else if(j == L_Y/2 && i < L_X/2)
                    grille[i][j] = +1;
                else
                    grille[i][j] = -1;
            }
            */
            if(j >= L_Y/2-7 and j <= L_Y/2+7 )
                grille[i][j]= 1;
            else
                grille[i][j]= -1;
        }
    }
}
