#ifndef	KMC_H
#define KMC_H
#include "Events.h"
#include "Island.h"
#include "Adatom.h"


const int n_classes = 6 ;


class KMC {

public:

KMC(const double,const double);


void init (const int , const int , const double , const double, const bool read_old = false);

void step(const double , const bool debug_mode  = false);

void print(int, int flag =0 ) const;
void print_final(int) const;

double get_concentration() const;
int* get_nevents() const;

private: 
int count_dnn0=0,count_dnn1=0,count_dnn2=0,count_dnn3=0,count_a=0,count_d=0;
int L;
int radius;
double  concentration;
//int get(); get method to retrieve L, radius or density

 double J;
 double A;
 double  current_T;

Adatom adatom;
Island island;
Events R[n_classes];
 
double det_rate(const int) const;
double att_rate() const;
double cumulative ( double * ) ;// calls rate evaluation, so is not a constant member function
int extract (const int ) const;// extracts event
bool is_attSite(const int , const int ) const;

};


#endif
