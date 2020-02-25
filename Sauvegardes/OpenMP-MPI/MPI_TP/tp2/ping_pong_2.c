/*
 * ping_pong_2.c   --- TP2 : Communications point à point : ping-pong
 *
 * Auteur          : Dimitri LECAS (CNRS/IDRIS - France) <dimitri.lecas@idris.fr>
 *
*/

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char *argv[]) {
  int rang,iter;
  int nb_valeurs=1000;
  int etiquette=99;
  double valeurs[nb_valeurs];
  MPI_Status statut;
  double temps_debut,temps_fin;

  /* ? */

  if (rang == 0) {
    for (iter = 0; iter<nb_valeurs; iter++)
      valeurs[iter] = rand() / (RAND_MAX + 1.);
    temps_debut=MPI_Wtime();
    /* ? */
    temps_fin=MPI_Wtime();
    printf("Moi, processus 0, j'ai envoye et recu %d valeurs"
           "(derniere = %g) du processus 1 en %f secondes.\n",
           nb_valeurs, valeurs[nb_valeurs-1], temps_fin-temps_debut);
    /* ? */
  }

  /* ? */
  return 0;
}
