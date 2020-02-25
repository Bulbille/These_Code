#include "dynamique.h"
#include "math.h"
#include "init.h"

/***********************************************/
/*** Mise à jour des liens et tirage au sort ***/
/***********************************************/
double exchange(int grille[][TAILLE_Y], int x1, int y1, int x2, int y2){
        double h_1 = grille[modulo(x1+1,TAILLE_X)][y1]+grille[modulo(x1-1,TAILLE_X)][y1]+grille[x1][y1+1]+grille[x1][y1-1];
        double h_2 = grille[modulo(x2+1,TAILLE_X)][y2]+grille[modulo(x2-1,TAILLE_X)][y2]+grille[x2][y2+1]+grille[x2][y2-1];
        double delta_h = 2.*J*grille[x1][y1]*(h_1-h_2);
	if((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) <= 1) delta_h += 4.*J;

	//// SHEAR ////
	/*Il faut savoir si l'on va dans le bon sens ou pas. */
	/*on va prendre comme référence les +1.*/
	if(y1 == y2){
		int bon_sens;
		(grille[x1][y1] == 1) ? bon_sens = 1 : bon_sens = -1;
		delta_h += bon_sens*(J*F + w*(y1 - (TAILLE_Y+1)/2.));
	}
        return ((exp(-BETA*delta_h) <1 ) ? exp(-BETA*delta_h) : 1);
}

/******* CONDITIONS PÉRIODIQUES EN X ****/
/******* CONDITIONS PÉRIODIQUES EN Y ****/
void maj_liens(int grille[][TAILLE_Y], int liens[], double acceptance_rate[], double* delta_t, int* nb_liens){
        *nb_liens = 0;
	*delta_t = 0;
        for(int y=0;y<TAILLE_Y-1;y++){
        for(int x=0;x<TAILLE_X;x++){
		int x2 = modulo(x+1,TAILLE_X);
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
                if(grille[x][y] != grille[x][y+1]){
                        liens[*nb_liens *4] = x;
                        liens[*nb_liens *4+1] = y;
                        liens[*nb_liens *4+2] = x;
                        liens[*nb_liens *4+3] = y+1;
			acceptance_rate[*nb_liens] = exchange(grille,x,y,x,y+1);
			*delta_t += acceptance_rate[*nb_liens];
                        *nb_liens += 1;
		}
        }
        }
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


/********** ÉNERGIE LOCALE *******/
double hamil_loc(int grille[][TAILLE_Y],int x,int y){
	x = modulo(x,TAILLE_X);
	y = modulo(y,TAILLE_Y);
        double energie = grille[modulo(x+1,TAILLE_X)][y]+grille[modulo(x-1,TAILLE_X)][y]+grille[x][modulo(y+1,TAILLE_Y)]+grille[x][modulo(y-1,TAILLE_Y)];
//        if(y == 0) energie -= grille[x][y-1];
 //       else if(y == TAILLE_Y -1 ) energie -= grille[x][y+1];
        return -J*energie*grille[x][y];
}

double champ_mag(int x, int y){
	if( x > (TAILLE_X/2-L_FAISCEAU) && x < (TAILLE_X/2+L_FAISCEAU)) return H;
	else return 0;
}

/***********************************/
/****** GÉNÉRATION DU RÉSEAU *******/
/***********************************/

void generation(int grille[][TAILLE_Y]){
        for(int j=0;j<TAILLE_Y;j++){
                for(int i=0;i<TAILLE_X;i++){
			(rand_01(generator) < 0.5 ) ? grille[i][j] = 1 : grille[i][j] = -1 ;
////////		if(TAILLE_Y % 2 ==0){
////////			if(j>=TAILLE_Y/2)
////////				grille[i][j] = -1;
////////			else
////////				grille[i][j] = 1;
////////		}
////////		else{
////////			if(j>TAILLE_Y/2)
////////				grille[i][j] = -1;
////////			else if(j == TAILLE_Y/2 && i < TAILLE_X/2)
////////				grille[i][j] = -1;
////////			else
////////				grille[i][j] = 1;
////////		}		
                }
        }
}

void generation(int grille[][TAILLE_Y], string& fichier){
        cout << "Generation : " << fichier << endl;
        stringstream tmp;
        tmp.str("");tmp << fichier.c_str();
        ifstream file(tmp.str());
        int x,y,spin;
        while(file >> x >> y >>  spin){
                grille[x][y]=spin;
        }
        file.close();
}
