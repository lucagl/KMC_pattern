#ifndef GLOBAL_H
#define GLOBAL_H
#include "functions.h"






// GLOBAL VARIABLES 
extern const double  PI; //A = attachment const, F=adatom diffusion constant
extern int proc_ID,n_proc;
extern unsigned seed,total_n_proc,n_threads;
extern unsigned* localseed;
extern const int root_process;



class FlatLand {
    protected :
        unsigned L;
        // double * g;
        fftw_complex *gft;
        fftw_plan FT,IFT;
        bool isConvinit = false;

    public:
        unsigned short ** matrix;
        int getBox()const {return L;};
        double ** gaussianConv() const;
        

        FlatLand(){};

        FlatLand(const int L_in){
            L = L_in;
            matrix =new unsigned short*[L];
            for (int i =0;i<L;i++){
                matrix[i] = new unsigned short [L] (); 
            }
        };
        void initConv(const double);

        void resetKernel(double sigma){
            if(isConvinit) {
                double * gk = (double*) malloc (sizeof(double) * L*L);
                double * g = (double*) malloc (sizeof(double) * L*L);
                g = gauss(L,sigma);
                fftShift(g,gk,L,L);
                fftw_execute_dft_r2c(FT,gk,gft);
            }
            else {std :: cout << "Cannot reset non initialised instance. ERROR"<< std ::endl; return;}
        };

        void saveTxt(const std::string& file_name, const double T, const double c) const {
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
        };




//COPY ASSIGNEMENT OPERATOR
        FlatLand& operator=(const FlatLand& oldobj){
            //std:: cout << "\n Base class copy assignement called"<< std :: flush; 
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
        };


        ~FlatLand(){
          //  std :: cout << "\n Flatland destroyer \n" << std :: flush;
            int Lf =L/2 +1;
            for(int i = 0; i < L; ++i) {
                delete [] matrix[i];
            }
            delete [] matrix;	

        //empty arrays related to convolution method
            if(isConvinit){
               // std :: cout << "\n fftw deallocation \n" << std :: flush; 
            fftw_free(gft);
            fftw_destroy_plan(FT);
            fftw_destroy_plan(IFT);
           }
        //fftw_deallocation..
        
        
  
        }
};


#endif
