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





// GLOBAL VARIABLES 
extern const double  J, PI,r0,A,F; //A = attachment const, F=adatom diffusion constant
extern int compute_every;
extern int proc_ID;
extern const int n_classes;
extern int L;



#endif
