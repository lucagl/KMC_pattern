#ifndef	FUNCTIONS_H
#define FUNCTIONS_H
#include "Events.h"
#include "Island.h"
#include "Adatom.h"



//--

// void read_input(int* L, double* T, double* F, int* pattern_type ,int* N_add, int* N_integration, int* Number_sim, int* diffusionOnly);

double energy(int nn, double T, double J);
int extract (int N);
double cumulative (Events R []);
//int non_zero(int M [L][L]);

void KMC_step(Events R[], Island& island, Adatom& adatom, double T, int* counter);
int is_attSite(int x,int y, const Island& island);

#endif
