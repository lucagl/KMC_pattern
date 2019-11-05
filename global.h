#ifndef GLOBAL_H
#define GLOBAL_H


#include <iostream>
#include <fstream>//for file writing-reading
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctime>
#include <cstdlib>

//Shared memory parallelisation
#include "omp.h"
//-----------------





// GLOBAL VARIABLES 
extern const double  PI; //A = attachment const, F=adatom diffusion constant
extern int proc_ID;
extern int n_proc;
extern unsigned total_n_proc;
extern const int root_process;



#endif
