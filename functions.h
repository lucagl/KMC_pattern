#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>//for file writing-reading
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctime>
#include <cstdlib>

#include<list>
#include<tuple>
#include<vector>
#include<algorithm>

//Shared memory parallelisation
#include <omp.h>
//-----------------
//Fast Fourier Transform
#include <fftw3.h>

void read_input(int* , double* , double* , int* , bool*, double* , double *, double *, int* , int*);

inline unsigned extract (unsigned N){
    if (N>RAND_MAX){
        std :: cout << "\n Problem, cannot extract random number bigger than " << RAND_MAX << "\n";
        exit(EXIT_FAILURE);
    }
    return (rand() % N);
};

double* gauss(const int, const double);

void fftShift(double*, double*, const int , const int, int x_shift=0, int y_shift=0);


 #endif