bool isOnlyDouble(const char* str);
void parametres(int argc, char* argv[]);
int modulo(int a, int b); 
double fact(int n);
double lnfact(int n);

/******* Fonctions utilitaires *****/
void parametres(int argc, char* argv[]){
    if(argc > 1){ 
        for(int i = 1; i<argc; ++i){
            std::string arg = argv[i];
            if(arg == "-t"){
                if(isOnlyDouble(argv[i+1])){
                    ttc = atof(argv[i+1]);
                    Beta = 1/(ttc*T_C);
                }
            }
            else if(arg == "-m"){
                if(isOnlyDouble(argv[i+1]))
                    mumax = atof(argv[i+1]);
            }
            else if(arg == "-prefix"){
                prefix = argv[i+1];
            }
            else if(arg == "-f"){
                if(isOnlyDouble(argv[i+1]))
                    drive = static_cast<double>(atof(argv[i+1]));
            }
            else if(arg == "-suffix"){
                suffix = argv[i+1];
            }
        }
    }   
    if(system(("mkdir -p " + prefix).c_str())) {;} 
}
/*** Modulo **/
int modulo(int a, int b){ 
    if(a<0) return a+b;
    else if(a>=b) return a-b;
    else return a;
}
double fact(int n) { 
    if ((n==0)||(n==1))
        return 1; 
    else
        return n*fact(n-1);
}
double lnfact(int n){
    if(n<25) 
        return log(fact(n));
    else 
        return n*log(n)-n;
}

bool isOnlyDouble(const char* str){
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;

    return true;
}


