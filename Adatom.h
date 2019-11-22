#ifndef ADATOM_H
#define ADATOM_H 
#include "global.h"




class Adatom: public FlatLand {


    public:
        
        int N;
        // bool is_attSite(int ,int );
        //void saveTxt (const std::string&, const double, const double) const;
        Adatom() {};
        Adatom(const int L_in , const double concentration =0){
            std :: cout << "Initialising adatom \n" << std :: flush;
            int current_density = 1;
            int rand_inti,rand_intj;

            L = L_in;
            matrix =new unsigned short*[L];
            for (int i =0;i<L;i++){
                matrix[i] = new unsigned short [L] (); 
            }

            if(concentration!=0){
                N = static_cast<int>(concentration*L*L);

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
        }
// //Misses copy assignement and copy constructior (rule of three)


//         ~Adatom(){
//         for(int i = 0; i < L; ++i) {
//             delete [] matrix[i];
//         }
//         delete [] matrix;	
//         }
};

#endif