
#include "Adatom.h"
#include "global.h"







void Adatom :: print(const std::string& file_name, int** mask, double T){
	
	
	std :: ofstream outfile (file_name);
	if (outfile.is_open()){
		
		// outfile << "#L\tT\tc \n";
		// outfile <<  L << "\t"<< T <<"\n";
		outfile <<"#T = " << T << "\n";
		outfile <<"#position \t \t attachment site\n";
		for ( int i = 0;i < L;i++){
			for(int j =0;j<L;j++){
				outfile << matrix[i][j]<<"\t" << mask[i][j]<< "\n"; //double(dx*height[int(i)]));
			}
		}
	}
	outfile.close();
}

// void Adatom :: initialize(double concentration ){

// 	int n_adatom, current_density = 0;
// 	int rand_inti,rand_intj;
//     conc = concentration;
// 	n_adatom = int(conc*L*L);

//     for(int i =0;i<L; i++){
		
//         for(int j =0;j<L;j++){
// 			if(current_density<=n_adatom){
// 				rand_inti=rand()% L;
// 				rand_intj=rand()% L;
// 				matrix[rand_inti][rand_intj] +=1; // can be on top of each other 
// 				current_density += 1;
// 		        }
// 		}
//     }
//   //flag = 0 ; // update all classes
// }
