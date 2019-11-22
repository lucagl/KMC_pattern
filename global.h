#ifndef GLOBAL_H
#define GLOBAL_H


#include <iostream>
#include <fstream>//for file writing-reading
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctime>
#include <cstdlib>

//Shared memory parallelisation
#include <omp.h>
//-----------------





// GLOBAL VARIABLES 
extern const double  PI; //A = attachment const, F=adatom diffusion constant
extern int proc_ID,n_proc;
extern unsigned seed,total_n_proc,n_threads;
extern unsigned* localseed;
extern const int root_process;



class FlatLand {
    protected :
        int L;

    public:
        unsigned short ** matrix;
        int getBox()const {return L;};
        

        FlatLand(){};

        //WRITE COPY ASSIGNEMENT OPERATOR!! THEY WILL BE DERIVED..

        FlatLand& operator=(const FlatLand& oldobj){
            std:: cout << "\n Base class copy assignement called"<< std :: flush; 
            L = oldobj.getBox();
            matrix =new unsigned short*[L];
            for (int i =0;i<L;i++){
                 matrix[i] = new unsigned short [L] ();
                 for (int j = 0; j < L; j++)
                 {
                     matrix[i][j] = oldobj.matrix[i][j];
                 }
             }
            return *this;
        }

        // FlatLand(const int L_in){
        //     L = L_in;
        //     matrix =new unsigned short*[L];
        //     for (int i =0;i<L;i++){
        //         matrix[i] = new unsigned short [L] (); 
        //     }
            
        // }

    void saveTxt(const std::string& file_name, const double T, const double c) const {
//print method of the class	

	std :: ofstream outfile (file_name);
	if (outfile.is_open()){
		outfile <<"#T = " << T <<"\tc ="<<c <<"\n";;
		outfile << "#Island\n";
		for ( int i = 0;i < L;i++){
			for(int j =0;j<L;j++){
				outfile << matrix[i][j]<< "\n"; 
			}
		}
	}
	outfile.close();
}
//Misses copy assignement and copy constructior (rule of three)
        ~FlatLand(){
        for(int i = 0; i < L; ++i) {
            delete [] matrix[i];
        }
        delete [] matrix;	
        }
};


#endif
