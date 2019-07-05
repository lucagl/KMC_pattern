#ifndef ISLAND_H
#define ISLAND_H 
#include "global.h"


class Island {
	//class members
private:
	int L;
public:
	bool ** matrix;
	unsigned short ** nn;
	// Island(){};
	// Island(const int L_){
	// 	std :: cout <<" Initializing island \n"<< std::flush ;
	// 	L = L_;
	// 	std :: cout << " L = " << L<< std::flush ;
	// 	matrix = new bool*[L];
	// 	nn = new unsigned short* [L];

	// 	for (int i =0;i<L;i++){ 
	// 		matrix[i] = new bool[L] ();
	// 		nn[i] = new unsigned short[L] ();
	// 	}
	// };

	void init(const int L_in, const int radius = 0){

		L = L_in;

		matrix = new bool*[L];
		nn = new unsigned short* [L];

		for (int i =0;i<L;i++){ 
			matrix[i] = new bool[L] ();
			nn[i] = new unsigned short[L] ();
		}
		
		if (radius!=0){
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
	}


    void print(const std::string&, unsigned short**, unsigned short**, unsigned short**,  const double, const double ) const;

	void init_neighbours ();
	int get_neighbours( const int  ,const  int  );

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


