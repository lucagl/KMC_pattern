


#include "global.h"
#include "Island.h"



//print method of the class



void Island :: print(const std::string& file_name, int** mask1, int** mask2, int** mask3, double T){
	
	
	
	std :: ofstream outfile (file_name);
	if (outfile.is_open()){
		
		
		//outfile << "#L\tT\n";
		//outfile <<  L << "\t"<< T << "\t" << "\n";
		outfile << "#Island\tsite1\tsite2\tsite3\n";
		for ( int i = 0;i < L;i++){
			for(int j =0;j<L;j++){
				outfile << matrix[i][j]<< "\t" << mask1[i][j]<< "\t" << mask2[i][j] << "\t" << mask3[i][j] << "\n"; 
			}
		}
	}
	outfile.close();
}



// void Island :: initialize( int radius){

//     int x0 = int(L/2);

//     for(int i =0;i<L; i++){		
//         for(int j =0;j<L;j++){
//            // if (sqrt((i-x0)*(i-x0)+(j-x0)*(j-x0))<=radius) {
// 			  if ((abs(i-x0)<=radius)&&(abs(j-x0)<=radius)) {
//                 matrix[i][j] = 1;
//             }
// 		}
//     }
//   //flag = 0 ; // update all classes
// }


//------------------------

void Island :: init_neighbours(){
	//only of the island


	int local_neighbour;
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
			nn[i][j] = (matrix[i][right]&matrix[i][j]) + (matrix[i][left]&matrix[i][j]) 
			+ (matrix[top][j]&matrix[i][j]) + (matrix[bottom][j]&matrix[i][j]);
		}
			
	}
}


int Island :: get_neighbours(int x, int y){
	//only of the island


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

	//nn[i][j] = local_neighbour;

return local_neighbour;

}
























// 	else if (R[3]<d_rand<R[4]){
// 		//detach 3 (3 neighbours)
// 		attached_site=pick_event(Natt);//just int random generator up to N
// 		x,y = sort(R_att,attached_site);//1st solution pick Nth non zero element of mask matrix
// 		//2nd solution: make the class R, which contains both position, and sort method --> most elegant
// 		flag = 4;//attachment event
// 		Natt += -1;
// 		adatom[y][x] -=1;
// 		island[y][x] +=1;
// 		nn_adatom = get_adatomIslandNN(x,y);//change to function not method if i will separate adatom and island classes
// 		if (nn_adatom){
// 			Rdet1.populate(x,y);/* new element added with associated coordinate 
// 			Note that this not exclude that 2 elements have same coordinate 
// 			(2 adatom on top of each other)*/
// 		}
		
// 	}

// }



