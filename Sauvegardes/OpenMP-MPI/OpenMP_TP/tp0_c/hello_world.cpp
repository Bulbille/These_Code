#include <iostream>
#include <omp.h>
#include <stdio.h>

using namespace std;

int main(){

	int nb_thread=0;
	int rang;

	#pragma omp parallel private(rang) 
	{
		rang=omp_get_thread_num();
		printf("Thread au rang : %d \n", rang);
		nb_thread++;
	}
	
	cout << "Max de thread : " << nb_thread << "\n";

return 0;}
