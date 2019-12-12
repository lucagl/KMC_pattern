#include "global.h"



//Define global variables
const double PI = 3.14159265358979323846;
const int root_process=0;
unsigned id,seed,n_threads;
int proc_ID,n_proc;
unsigned* localseed;

void FlatLand :: initConv(const double sigma){
    int Lf = L/2+1;
    double norm = 1./(double(L)*double(L));
    isConvinit = true;
    const int n_threadsFTT = 4;

    
    double * g = (double*) malloc (sizeof(double) * L*L);
    

    fftw_init_threads();
    fftw_plan_with_nthreads(n_threadsFTT);

    //gft = (fftw_complex*) fftw_malloc(L*Lf*sizeof(fftw_complex));
    gft = fftw_alloc_complex (L * Lf);//real to complex. The complex array has size N_x x (N_y/2 +1) 
    //gft[i][j] = gft[j+(L/2+1)*i] with j = 0..L/2 and i = 0..(L-1)
    //FFT initializations see for 2d
    //at this stage g and gft are dummy containers of the correct size

    FT = fftw_plan_dft_r2c_2d (L,L,g,gft,FFTW_ESTIMATE);//plan for real to complex
    IFT = fftw_plan_dft_c2r_2d (L,L,gft,g,FFTW_ESTIMATE);//plan for complex to real        
    
    gauss(g, L,sigma);

    //std :: cout << "\n Here 1" << std::flush;

    fftw_execute_dft_r2c(FT,g,gft); //forward fourier transform stored in gft. 

   // fftw_execute_dft_c2r(IFT,gft,g);

   //  std :: cout << "\n Here 2" << std::flush;

    // CHECK
    // std :: string file_name = "dummy/gaussianK_fromInverse.txt";
    // std :: ofstream outfile (file_name);
    // if (outfile.is_open()){
    //     for ( int i = 0;i < L;i++){
    //         for(int j =0;j<L;j++){
    //             outfile << norm * g[j+L*i]<< "\n"; 
    //         }
    //     }
    // }

    // outfile.close();

    //std :: cout << "\n Here 3" << std::flush;

            
}



double** FlatLand :: gaussianConv()const {
    int Lf = L/2+1;
    double norm = 1./(double(L)*double(L));
    double** smooth_matrix = new double*[L];
    for (int i =0;i<L;i++){
        smooth_matrix[i] = new double [L] ();
    }

    if(!isConvinit) {std :: cout <<"Convolution was not initialized. ERROR."<<std :: endl; return smooth_matrix;}
    
    
    double *direct = (double*) malloc (sizeof(double) * L*L);
    fftw_complex *fSpace = fftw_alloc_complex (L * Lf);
    fftw_complex *convolft = fftw_alloc_complex (L * Lf);


    double * g = (double*) malloc (sizeof(double) * L*L);
    gauss(g, L,5);


    fftw_execute_dft_r2c(FT,g,gft);

    #pragma omp parallel 
    {
        #pragma omp for
        for (int i = 0; i < L*L; i++)
        {//create array of contiguous memory locations
            direct[i] = matrix[i/L][i%L];//i=i%L+(i/L)*L
        }


//check-------------------
    
    // #pragma omp single
    // {
    // std :: string file_name = "dummy/final.txt";
    // std :: ofstream outfile (file_name);
    // for (int i = 0; i < L; i++)
	// {
	// 	for (int j = 0; j < L; j++)
	// 	{
	// 		outfile << direct[j+i*L]<< "\n"; 
	// 	}
		
	// }
	// outfile.close();
    // // }
// ------------------------------


        #pragma omp single 
        fftw_execute_dft_r2c(FT,direct,fSpace);

        #pragma omp for
        for (int i = 0; i < L*Lf; i++){
            // std :: complex<double> dummy1,dummy2;
            // std :: memcpy(&dummy1,gft[i],sizeof( fftw_complex ));
            // std :: memcpy(&dummy2,fSpace[i],sizeof( fftw_complex ));
            //[j+(L/2+1)*i] with j = 0..L/2 and i = 0..(L-1) i%Lf+(i/Lf)*Lf


            convolft[i][0] =  gft[i][0] * fSpace[i][0] - gft[i][1] * fSpace[i][1];  // real part
            convolft[i][1] =  gft[i][0] * fSpace[i][1] + gft[i][1] * fSpace[i][0];  // imaginary part
        }

        #pragma omp single
        fftw_execute_dft_c2r(IFT,convolft,direct);
        //fftw_execute_dft_c2r(IFT,fSpace,direct);//to check if inverse equal to direct

        #pragma omp for
        for (int i = 0; i < L*L; i++)
        {
            smooth_matrix[i/L][i%L] = norm* direct[i];//i=i%L+(i/L)*L
        }
    } 


// 	std :: string file_name1 = "dummy/final_smooth.txt";
//     std :: ofstream outfile1 (file_name1);

// // check -----------------------------------
// 	for (int i = 0; i < L; i++)
// 	{
// 		for (int j = 0; j < L; j++)
// 		{
// 			outfile1 << norm*direct[j+i*L]<< "\n"; 
// 		}
		
// 	}
// 	outfile1.close();


double * a = (double*) malloc (sizeof(double) * L*L);
fftw_execute_dft_c2r(IFT,gft,a); // is a call of the inverse modyfying also first argument?


std :: string file_name = "dummy/gaussianK_fromInverse.txt";
    std :: ofstream outfile (file_name);
    if (outfile.is_open()){
        for ( int i = 0;i < L;i++){
            for(int j =0;j<L;j++){
                outfile << norm * a[j+L*i]<< "\n"; 
            }
        }
    }

    outfile.close();


    //---------------------------

    return smooth_matrix;

}
