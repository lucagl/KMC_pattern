#ifndef ISLAND_H
#define ISLAND_H 
#include "global.h"


class Island {
	//class members
private:
	int L;
public:
	bool ** matrix;
	unsigned short ** nn1;
	unsigned short ** nn2;

	void init(const int L_in, const bool circle =0, const int radius = 0){
	//default value in order to be able to initialise the object given the total box size when reading from previous integration file
		L = L_in;

		matrix = new bool*[L];
		nn1 = new unsigned short* [L];
		nn2 = new unsigned short* [L];

		for (int i =0;i<L;i++){ 
			matrix[i] = new bool[L] ();
			nn1[i] = new unsigned short[L] ();
			nn2[i] = new unsigned short[L] ();
		}
		
		if (radius!=0){
			int x0 = int(L/2);
			if(circle == false){
				for(int i =0;i<L; i++){		
					for(int j =0;j<L;j++){
						if ((abs(i-x0)<=radius)&&(abs(j-x0)<=radius)) {
							matrix[i][j] = 1;
						}
					}
				}
			}
			else{
				for(int i =0;i<L; i++){		
					for(int j =0;j<L;j++){
						if (sqrt((i-x0)*(i-x0)+(j-x0)*(j-x0))<=radius){
							matrix[i][j] = 1;
						}
					}
				}
			}	
		}
	}


    void print(const std::string&, const double, const double ) const;

	void init_neighbours ();//init first and second neighbours
	int get_neighbours1( const int  ,const  int  );
	int get_neighbours2( const int  ,const  int  );

	// destructor
	~Island(){
		for(int i = 0; i < L; ++i) {
			delete [] matrix[i];
			delete [] nn1[i];
			delete [] nn2[i];
		}
		delete [] matrix;
		delete [] nn1;	
		delete [] nn2
		;	
	}
};



#endif


