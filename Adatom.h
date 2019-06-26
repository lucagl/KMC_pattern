#ifndef ADATOM_H
#define ADATOM_H 
#include "global.h"

class Adatom {
    private :

    public:
        int ** matrix;
        int maxL;
        int N;
        bool is_attSite(int x,int y);
        void print(const std::string&, int** mask, double temperature);

        Adatom(double concentration){
            int n_adatom, current_density = 1;
            int rand_inti,rand_intj;

            N = int(concentration*L*L);
            maxL = L;

            matrix =new int*[L];
            for (int i =0;i<L;i++){
                matrix[i] = new int[L] (); 
            }
            
            
            for(int i =0;i<L; i++){
		
                for(int j =0;j<L;j++){
			        if(current_density<=N){
                        rand_inti=rand()% L;
                        rand_intj=rand()% L;
                        matrix[rand_inti][rand_intj] +=1; // can be on top of each other 
                        current_density += 1;
		            }
		        }
            }
        }



        ~Adatom(){
        for(int i = 0; i < L; ++i) {
            delete [] matrix[i];
        }
        delete [] matrix;	
        }
};

#endif