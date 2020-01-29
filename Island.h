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
	//	std :: cout << "Initialising island \n" << std :: flush;
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
		
		if (radius >= L){
			//vertical bar
			//std :: cout << "\nVertical bar\n "<< std :: endl;
			for(int i =0;i<L; i++){		
					for(int j =int(L/2-L/5);j<int(L/2+L/5);j++){
						matrix[i][j] = 1;
					}
				}
		}
		else{
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
	};



// copy assignement operator

 Island& operator=(const Island& oldobj){
            //std:: cout << "\n Island class copy assignement called"<< std :: flush; 
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

	void  init_neighbours(){
		//only of the island

		int left, right,bottom,top;
		
		for (int i = 0; i < L; i++){
			top = i+1;	
			if(top==L) top = 0;
			bottom = i-1;
			if(bottom==-1) bottom = L-1;
		
			for (int j = 0; j < L; j++)
			{
				right = j+1;
				if (right ==L) right = 0;
				left = j-1;
				if(left == -1) left = L-1; 
				nn1[i][j] = (matrix[i][right]&matrix[i][j]) + (matrix[i][left]&matrix[i][j]) 
				+ (matrix[top][j]&matrix[i][j]) + (matrix[bottom][j]&matrix[i][j]);
				nn2[i][j] = (matrix[top][right]&matrix[i][j]) + (matrix[top][left]&matrix[i][j]) 
				+ (matrix[bottom][right]&matrix[i][j]) + (matrix[bottom][left]&matrix[i][j]);
			}
				
		}
	}


	int get_neighbours1(const int x,const int y) const{


		int local_neighbour;
		int left, right,bottom,top;
		
		
		
		local_neighbour = 0;

		top = y+1;
		if(top==L) top = 0;

		bottom = y-1;
		if(bottom==-1) bottom = L-1;

		right = x+1;
		if (right ==L) right = 0;

		left = x -1;
		if(left == -1) left = L-1; 

		local_neighbour = (matrix[y][right]&matrix[y][x]) + (matrix[y][left]&matrix[y][x])
		+ (matrix[top][x]&matrix[y][x]) + (matrix[bottom][x]&matrix[y][x]);

		nn1[y][x] = local_neighbour;

	return local_neighbour;

	}

	int  get_neighbours2(const int x,const int y) const {


		int local_neighbour;
		int top, bottom,left,right;
		
		
		
		local_neighbour = 0;

		top = y+1;
		if(top==L) top = 0;

		bottom = y-1;
		if(bottom==-1) bottom = L-1;

		right = x+1;
		if (right ==L) right = 0;

		left = x -1;
		if(left == -1) left = L-1; 

		local_neighbour = (matrix[top][right]&matrix[y][x]) + (matrix[top][left]&matrix[y][x])
		+ (matrix[bottom][right]&matrix[y][x]) + (matrix[bottom][left]&matrix[y][x]);

		nn2[y][x] = local_neighbour;

	return local_neighbour;

	}

};
#endif


