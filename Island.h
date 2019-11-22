#ifndef ISLAND_H
#define ISLAND_H 
#include "global.h"


class Island : public FlatLand{
	//class members

public:

	unsigned short ** nn1;
	unsigned short ** nn2;
	Island(){};
	Island(const int L_in, const bool circle =0, const int radius = 0){
	//default value in order to be able to initialise the object given the total box size when reading from previous integration file
		std :: cout << "Initialising island \n" << std :: flush;
		L = L_in;
		// std :: cout << "\n L=" << L << std :: flush;
		// std :: cout << "\n radius=" << radius << std :: flush;
		matrix = new unsigned short*[L];
		nn1 = new unsigned short* [L];
		nn2 = new unsigned short* [L];

		for (int i =0;i<L;i++){ 
			matrix[i] = new unsigned short[L] ();
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

// copy assignement operator

 Island& operator=(const Island& oldobj){
            std:: cout << "\n Island class copy assignement called"<< std :: flush; 
            L = oldobj.getBox();
            nn1 =new unsigned short*[L];
			nn2 =new unsigned short*[L];
			matrix =new unsigned short*[L];
            for (int i =0;i<L;i++){
                 nn1[i] = new unsigned short [L] ();
				 nn2[i] = new unsigned short [L] ();
				 matrix[i] = new unsigned short [L] ();
				 
                 for (int j = 0; j < L; j++)
                 {
					 matrix[i][j] = oldobj.matrix[i][j];
                     nn1[i][j] = oldobj.nn1[i][j];
					 nn2[i][j] = oldobj.nn2[i][j];
                 }
             }
            return *this;
        }

    //void print(const std::string&, const double, const double ) const;

	void init_neighbours ();//init first and second neighbours
	int get_neighbours1( const int  ,const  int  );
	int get_neighbours2( const int  ,const  int  );

	// destructor
	~Island(){
		for(int i = 0; i < L; ++i) {
			// delete [] matrix[i];
			delete [] nn1[i];
			delete [] nn2[i];
		}
		// delete [] matrix;
		delete [] nn1;	
		delete [] nn2
		;	
	}
};



#endif


