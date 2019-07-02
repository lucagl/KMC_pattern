#ifndef ISLAND_H
#define ISLAND_H 
#include "global.h"

typedef const int * ptr_to_const_int;


class Island {
	//class members
private:

public:
	int ** matrix;
	int ** nn;
	int L;

    void print(const std::string&, int**, int**, int**,  double temperature);

	void init_neighbours ();
	int get_neighbours( const int x ,const  int y );

	

	//=====================================
	//double KMC_step(const double,const double, int*);
	
	//class member constructor and destructor
	
	//constructor
	void init(const int radius, const int L_in){

		L = L_in;

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


