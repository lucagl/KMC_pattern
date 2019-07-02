#ifndef	KMC_H
#define KMC_H
#include "Events.h"
#include "Island.h"
#include "Adatom.h"


const int n_classes = 5 ;


class KMC {

public:

KMC(const double,const double,const double);


void initialize (const int , const int , const double , const double);

void step(const double , int * , const bool debug_mode  = false);

void print(int frame, int flag =0 );

private: 

int L;
int radius;
double  concentration;
//int get(); get method to retrieve L, radius or density

 double J;
 double A;
 double F;
 double  current_T;

Adatom adatom;
Island island;
Events R[n_classes];

double energy(const int);
double cumulative ( double * );
int extract (const int );// extracts event
bool is_attSite(const int , const int );

};


#endif
