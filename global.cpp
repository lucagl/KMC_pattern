#include "global.h"



//Define global variables
const double PI = 3.14159265358979323846;
const int root_process=0;
unsigned id,seed,n_threads, max_threads,n_threadsFTT;
int proc_ID,n_proc;
unsigned* localseed;
//_______________________________



void FlatLand :: initConv(const double sigma){
    int Lf = L/2+1;
    double norm = 1./(double(L)*double(L));

    if(isConvinit) {resetKernel(sigma); return;}

    isConvinit = true;
    
   // std :: cout <<"Initialization of convolution routines"<< std :: endl;
    double * gk = (double*) malloc (sizeof(double) * L*L);
    double * g = (double*) malloc (sizeof(double) * L*L);
    //g contains gaussian centered in the center (L/2)
    g=gauss(L,sigma);
    fftShift(g,gk,L,L);
    gft = fftw_alloc_complex (L * Lf);//real to complex. The complex array has size N_x x (N_y/2 +1) 

    fftw_init_threads();
    n_threadsFTT = max_threads;
    fftw_plan_with_nthreads(n_threadsFTT);
    

    FT = fftw_plan_dft_r2c_2d (L,L,gk,gft,FFTW_ESTIMATE);//plan for real to complex
    IFT = fftw_plan_dft_c2r_2d (L,L,gft,gk,FFTW_ESTIMATE);//plan for complex to real        

    fftw_execute_dft_r2c(FT,gk,gft); //forward fourier transform stored in gft. 
    /*CAREFUL : do not back transofrm gft giving as argument to c2r gft. 
    The function affects the first argument..*/
    free(g);
    free(gk);



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

    
    #pragma omp parallel num_threads(n_threadsFTT)
    {
        #pragma omp for 
        for (int i = 0; i < L*L; i++)
        {//create array of contiguous memory locations
            direct[i] =  matrix[i/L][i%L];
        }

         #pragma omp single 
        fftw_execute_dft_r2c(FT,direct,fSpace);

        #pragma omp  for
        for (int i = 0; i < L*Lf; i++){
            convolft[i][0] =  gft[i][0] * fSpace[i][0] - gft[i][1] * fSpace[i][1];  // real part
            convolft[i][1] =  gft[i][0] * fSpace[i][1] + gft[i][1] * fSpace[i][0];  // imaginary part
        }

        #pragma omp single
        fftw_execute_dft_c2r(IFT,convolft,direct);//convolf overwritten careful
        
        //fftw_execute_dft_c2r(IFT,fSpace,direct);//to check if inverse equal to direct

        #pragma omp for
        for (int i = 0; i < L*L; i++)
        {
            smooth_matrix[i/L][i%L] = norm* direct[i];//i=i%L+(i/L)*L
        }
    } 

    fftw_free(convolft);
    fftw_free(fSpace);
    free (direct);


    return smooth_matrix;

}

// Routines to dispose fft routines 