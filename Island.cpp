#include "global.h"
#include "Island.h"




// void Island :: saveTxt(const std::string& file_name, const double T, const double c) const {
// //print method of the class	

// 	std :: ofstream outfile (file_name);
// 	if (outfile.is_open()){
		

// 		outfile <<"#T = " << T <<"\tc ="<<c <<"\n";;
// 		outfile << "#Island\n";
// 		for ( int i = 0;i < L;i++){
// 			for(int j =0;j<L;j++){
// 				outfile << matrix[i][j]<< "\n"; 
// 			}
// 		}
// 	}
// 	outfile.close();
// }

// double** Island :: GaussianC(double sigma){

// 	#pragma omp ecc

// };


void Island :: init_neighbours(){
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


int Island :: get_neighbours1(const int x,const int y){


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

int Island :: get_neighbours2(const int x,const int y){


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

