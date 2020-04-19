#ifndef ISLAND_H
#define ISLAND_H 
#include "global.h"


class Island : public FlatLand{
	//class members

public:

	unsigned short ** nn1;
	unsigned short ** nn2;
	Island(){};
	Island(const unsigned L_in, const bool circle =0, const int radius = 0): FlatLand(L_in)
	{
	//default value in order to be able to initialise the object given the total box size when reading from previous integration file
	    // std :: cout << "\n Initialising island \n\n" << std :: endl;
		// std :: cout << "Island L="<< L << std:: endl;
		//matrix = new unsigned short*[L];
		
		nn1 = new unsigned short* [L];
		nn2 = new unsigned short* [L];
		for (int i =0;i<L;i++){ 
			//matrix[i] = new unsigned short[L] ();
			nn1[i] = new unsigned short[L] ();
			nn2[i] = new unsigned short[L] ();
		}
		
		if (radius >= L){
			//Multiple vertical lines
			//estimate how many bands fit..
			bool diagonal = false;	
			const int nbands = 5;
			if (proc_ID ==0)
			{
				std :: cout << "\n Straight line 0, diagonal lines (45°) 1. Insert value \n" << std:: endl;
				std:: cin >> diagonal;
				std :: cout << "\n Inserted value : " << diagonal;

			}
			
			
		
			if(!diagonal){
				std :: cout << "\n Vertical bands \n" << std:: endl;
				int thickness =int ((double)L/(2*nbands));

			for(int i =0;i<L; i++){	
					int n = 0;
					bool flag=0;
					for(int j =int(thickness/2);j<L-int(thickness/2);j++){
						if(n%thickness==0) {flag = !flag; }//std:: cout << n<<"SWITCH"<< flag <<"\n";}
						matrix[i][j] = flag;
						n+=1;	
					}
				}
			}
			else{
				std :: cout << "\n Diagonal bands \n" << std:: endl;
				int thickness = int (sqrt(2) * (double) L/(2*nbands));
				bool flag=0;
				// int switchx = (int)((double) thickness/cos(PI/4));
				int offset = ceil(((double) thickness/sin(PI/4))); 
				// int nx,ny;
				// bool flagx,flagy;
				// ny = 0;
				// nx = -1;
				// // flag = 0;
				// flagx=0;
				// flagy=0;
				// bool swap = 0;
				// //int s=0;
				// for(int i =0;i<L; i++){	
				// 	//nx= -ny + offset;
				// 	ny +=1;
				// 	if(abs(ny)%offset==0) {flagy=!flagy; }
				// 	// std:: cout << "\n"<< nx<<"\t"<<ny << std::endl;
				// 	nx =-ny;
				// 	flagx = 0;
				// 	for(int j= 0;j<L;j++){
				// 		if(nx%(offset)==0) {flagx= !flagx;}
				// 		//if(nx%(offset/2)==0) {swap= !swap;}
				// 		//matrix[((L-1)-i)*swap + i*!swap][((L-1)-j)*swap + j*!swap] = flagx^flagy;//xor and offset identical for many bands
				// 		matrix[i][j] = flagx^flagy;
				// 		nx +=1;	
				// 		//i = y, j = x
				// 	}
				// }
				//draw lines and the fill them
				int x0 = 0;
				while(x0<L){
					int x,y;
						for(int j= 0;j<L;j++){
							x = j;
							y = -x +x0 ;
							if(y>=L) {//std::cout <<"\n HERE" << L <<"\t" << y <<std::endl;
							 y=-L+y;  };
							if(y <0) {//std::cout <<"\n HERE" << L <<"\t" << y <<std::endl; 
							y= L + y; };
							//std::cout <<"\n"<< y <<std ::endl;
							matrix[y][x] = flag;
						}
						x0+=1;
						if(x0%offset==0) flag=!flag;
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
	 //Don't know why it does not call copy assignment of mother class also..
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


