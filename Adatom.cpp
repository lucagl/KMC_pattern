
#include "Adatom.h"
#include "global.h"




void Adatom :: print(const std::string& file_name, const double T, const double c) const{
	
	
	std :: ofstream outfile (file_name);
	if (outfile.is_open()){
		
		// outfile << "#L\tT\tc \n";
		// outfile <<  L << "\t"<< T <<"\n";
		outfile <<"#T = " << T << "\tc ="<<c <<"\n";
		outfile <<"#position\n";
		for ( int i = 0;i < L;i++){
			for(int j =0;j<L;j++){
				outfile << matrix[i][j]<< "\n"; 
			}
		}
	}
	outfile.close();
}


// void Adatom :: print(const std::string& file_name, unsigned short** mask, const double T, const double c) const{
	
	
// 	std :: ofstream outfile (file_name);
// 	if (outfile.is_open()){
		
// 		// outfile << "#L\tT\tc \n";
// 		// outfile <<  L << "\t"<< T <<"\n";
// 		outfile <<"#T = " << T << "\tc ="<<c <<"\n";
// 		outfile <<"#position \t \t attachment site\n";
// 		for ( int i = 0;i < L;i++){
// 			for(int j =0;j<L;j++){
// 				outfile << matrix[i][j]<<"\t" << mask[i][j]<< "\n"; //double(dx*height[int(i)]));
// 			}
// 		}
// 	}
// 	outfile.close();
// }

