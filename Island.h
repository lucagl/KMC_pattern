#ifndef ISLAND_H
#define ISLAND_H 
#include "global.h"


class Island {
	//class members
private:

public:
	int ** matrix;
	int ** nn;
	int maxL;

    void print(const std::string&, int** mask1, int** mask2, int** mask3,  double temperature);

	void init_neighbours ();
	int get_neighbours( int x , int y );

	

	//=====================================
	//double KMC_step(const double,const double, int*);
	
	//class member constructor and destructor
	
	//constructor
	Island(int radius){

		maxL = L;

		matrix = new int*[L];
		nn = new int* [L];

		for (int i =0;i<L;i++){ 
			matrix[i] = new int[L] ();
			nn[i] = new int[L] ();
		}

		int x0 = int(L/2);

    	for(int i =0;i<L; i++){		
			for(int j =0;j<L;j++){
			// if (sqrt((i-x0)*(i-x0)+(j-x0)*(j-x0))<=radius) {
				if ((abs(i-x0)<=radius)&&(abs(j-x0)<=radius)) {
					matrix[i][j] = 1;
				}
			}
		}

	}
	
	// destructor
	~Island(){
		for(int i = 0; i < L; ++i) {
			delete [] matrix[i];
			delete [] nn[i];
		}
		delete [] matrix;
		delete [] nn;	
	}
};



#endif


