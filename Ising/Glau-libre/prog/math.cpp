#include "math.h"
#include "init.h"


/******* Fonction modulo ********/
int modulo(int a, int b){
        while(a < 0 or a >= b){
                if(a<0) a+= b;
                else if(a >= b) a-=b;
        }
        return a;
}

double modulo_d(double a, int b){
        while(a < 0 or a >= b){
                if(a<0) a+= b;
                else if(a >= b) a-=b;
        }
        return a;
}
/**** VERIFICATION DOUBLE *****/
bool isOnlyDouble(const char* str)
{
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

	return true;
}

