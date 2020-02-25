#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;
#define PRECISION 10e-4
#define COL 3
void pas(double *dx, double *dy, double dl, double ex, double ey){
	double e = sqrt( ex*ex + ey*ey);
	*dx = ex*dl/e;
	*dy = ey*dl/e;

}
void champ_electrique(double x, double y, double *ex, double *ey, double tableau[][COL], int taille  ){
	int j;
	double champ_x = 0 , champ_y = 0;

	for(j=0; j<taille;j++){	
		champ_x += tableau[j][2]*(x-tableau[j][0])/pow(  ( (x-tableau[j][0])*(x-tableau[j][0]) + (y- tableau[j][1])*(y-tableau[j][1])  ), 1.5);
		champ_y += tableau[j][2]*(y-tableau[j][1])/pow(  ( (x-tableau[j][0])*(x-tableau[j][0]) + (y- tableau[j][1])*(y-tableau[j][1])  ), 1.5);
	}

	*ex = champ_x;
	*ey = champ_y;
}

int main(void){

	double x, y, dx,dy, ex, ey, theta;
	double dl = 0.001;
	int i, j, m=15;
	int taille, k;
	ofstream resultats("ligne_champ.res");


	cout << "Insérer le nombre de charges que vous voulez : " ;
	cin >> taille;
	double table[taille][COL], distance[taille];

	for(j=0; j<taille; j++){
		cout << endl << "Écrire x : ";
		cin >> table[j][0];
		cout << "Écrire y : ";
		cin >> table[j][1];
		cout << "Écrire le signe de la charge (-1 = -, 1=+) : ";
		cin >> table[j][2];
	}

	for(i=0; i<taille; i++){

		if(table[i][2] <= 0) continue;

		for(j=0; j<=m ; j++){
			theta = 2*3.141592654*(j-1)/m;
			x = table[i][0] + 0.001*cos(theta);
			y = table[i][1] + 0.001*sin(theta);
//			cout << i << " " << j << " " << theta << endl;

			for(k=0; k<= 20000;k++){
				champ_electrique( x, y, &ex, &ey, table, taille );
				pas(&dx, &dy, dl, ex, ey);
				x += dx; y+=dy;

				if(k%10 == 0) resultats << x << " " << y << endl;

			}

			resultats << endl;
		}
	}
	return 0;
}
