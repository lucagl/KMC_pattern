#ifndef	KMC_H
#define KMC_H
#include "Events.h"
#include "Island.h"
#include "Adatom.h"


const int n_classes = 5 ;


class KMC {

public:

KMC(const double,const double,const double);


void init (const int , const int , const double , const double, const bool read_old = false);

void step(const double , int * , const bool debug_mode  = false);

void print(int, int flag =0 ) const;
void print_final(int ) const;

double get_concentration() const;

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
 
double energy(const int) const;
double cumulative ( double * ) ;// calls rate evaluation, so is not a constant member function
int extract (const int ) const;// extracts event
bool is_attSite(const int , const int ) const;

};


#endif
