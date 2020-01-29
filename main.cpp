/*######################### 2D SINGLE LAYER CLUSTER GROWTH COUPLED TO A CONCENTRATION FIELD ######################

 Author:								Luca Gagliardi IIT
 Copiright:							© 2019 Gagliardi IIT, ILM-CNRS
 
 
###############################################################################################

------------ COMPILING INSTRUCTIONS ----------------

mpic++ *.cpp -o kmc.ex -lfftw3_omp -lfftw3 -fopenmp

-Build libraries using..


-Make a python scritp--> Input file read by the script (make a user interface ? :) )
					  "Live plots" when calling render command
					  MPI must be displaced on python side? Not sure, I can still return observations which are results of averages. 
					  However I think is a better idea to parallelize from python to have more explicit control.



------------------- CONVENTIONS -----------
System of coordinates is x: left-right-wise and y: top-bottom wise. The indexes of arrays start with 0.

-------------------- TODO -----------------
- [minor]displace kmc.init() within constructor?
- Redefine all 2d arrays in contiguous memory to avoid copying them in temporal variable for fourier transform
	(However if the convolution is done just a few times should not be too heavy)

-Check coherency variable types everywhere (es. int and unsigned..)



/////
- MAJOR :Time dependent Temperature

- Program architecture can be improved.. 
	for instance L in both master and dependent classes is redundant.
	Creating a parent class for island and adatom could also be interesting
- Special class for diffusion that since is simultaneous has some specificities..

-------------------- COMMENTS ------------------
- Specific implementations for special scenario of NN1 =1, NN2 =0 and NN1 =0, NN2 =1.
  This because the fact that I have only on neighbour can be exploited for more efficient update.
- The isotropic case is BR=sqrt(2)/2 ~ 0.70711 (zeta = J/J')
-------------------- QUESTIONS ------------------

- Change attachment site criteria based on being on the diagonal ? --> No but tempting..

################################################################################
*/


// #############################################
//					MAIN	
// #############################################

#include "kmc.h"

// --- PARALLELISATION
#include <thread>
#include <mpi.h>

//------

int ierr;

/*  Need:
	methods to initialise and init
		env.init() ok
		env.reset() ok
	methods to print
		env.render()
	methods for extracting informations:
		env.action_space.sample() ..
	step:
	observation, reward, done, info = env.step(action)
So what I called kmc could become environment.. However I need to set number of episodes..
*/
int main(int argc, char **argv){

// int argc; //MPI stuff
// char **argv; //MPI stuff

const double J =0.2;

double T0,conc0, A, BR, E_shift;
int L,radius, n_steps;
int frame, print_every;

bool read_old,is_circle;

std::time_t start, end,curr_time;
clock_t t1,t2;

float seconds;

long elapsed_time;

t1=clock();
start = std::time(NULL);
	

ierr = MPI_Init(&argc, &argv);

MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

MPI_Comm_rank(MPI_COMM_WORLD, &proc_ID);


// ----------------- READ INPUT AND PRINT INITIAL INFO -------------
if(proc_ID == root_process){
	// read from input file and broadcast
	time(&curr_time);
	std :: cout << "\n Start time " << ctime(&curr_time) << "\n" << std :: flush;	
	std :: cout << "\n Number of processors= "<< n_proc << std :: endl ;

	read_input(&L, &T0,& conc0, &radius,&is_circle, &A, &BR, &E_shift, &n_steps, &print_every);//, &read_old);

	std :: cout << "\n J= " << J << "  |  L= "<< L<< "  |  T=" << T0 <<"  |  concentration= "<< conc0 <<
		"  |  initial island radius= "<< radius <<  "  |  attachment parameter= " << A << "	|	Energy shift= " << E_shift << "	|	Bond energy ratio= "<< BR <<"  |  kmc steps= " << n_steps<<
	"  |  print each= " << print_every <<"\n";
  
	double c_eq = exp((-2*J*(1+BR) + E_shift)/T0);
	
	std :: cout << "\n Equilibrium concentration at T=0, is  " << c_eq << "\n";
}

seed = time(NULL)*(proc_ID+1);

int syst_cores = std :: stoi(exec("grep -c ^processor /proc/cpuinfo"));

//std :: cout << syst_cores << std :: endl;
#pragma omp parallel num_threads(syst_cores)
{
	#pragma omp single
	{
	max_threads = omp_get_num_threads();
	// max_threads = omp_get_max_threads(); Use this if threads are set from environment variable
	if(proc_ID == root_process){
		std :: cout << "\n Max number of threads per process used = " << max_threads << "\n \n";
	}
	localseed = new unsigned[max_threads];//one seed per potential thread
	}
	
	unsigned id = omp_get_thread_num();	

	localseed[id] = seed *(id + 1);
	
}
	

//RootID broadcast data to other processors
MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&T0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&conc0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&radius, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&is_circle, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&BR, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&E_shift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&n_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&print_every, 1, MPI_INT, 0, MPI_COMM_WORLD);
//MPI_Bcast(&read_old, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD); skip for now, I wnat to handle this directly in python

// make directories for each thread

auto aString = "mkdir plots" + (std::to_string(proc_ID));
const char* makeDir = aString.c_str();
auto path = "plots" + (std::to_string(proc_ID));
auto bString ="rm " + path + "/island* " + path + "/adatom* ";
const char* remove_old = bString.c_str();
system(makeDir);
system(remove_old);

// -------------------------------------------------------


// 	___________________INITIALIZATION _____________________________
	
	srand (seed);// initialise random generator differently for each MPI thread
	
	KMC kmc(J,BR,A,E_shift,L,is_circle,radius,conc0,T0);
	kmc.init();
	kmc.initConv_adatom(double(L)/20);
	kmc.initConv_island(2);

//_____________________________ EPISODES ______________________________

// int e=0;
// while(){
// e+=1;

//kmc.reset();

// _________________________RUN KMC ___________________________

if(proc_ID == root_process){
	std :: cout << "\n *Starting integration* \n" << std :: endl;
}

kmc.saveTxt(path,0,true,true);
//kmc.saveTxt(path,0);
frame = 0;
double t =0;
for (int k = 1; k <= n_steps; k++){

	/* Temperature can be changed here..
	T = T0 + ..
	A = A0 + ..
	*/
	// EVOLUTION STEP 
	
	t+=kmc.step(T0);
	
	if ((k%print_every)== 0){
		frame+=1;
<<<<<<< HEAD
		kmc.print(frame);

		//update number threads
		//Parallelization seems conveninet only for sizes about 200X200 onwards.
		//Make it non parallelized if using sistemaically smaller sizes..
		n_threads= ceil(float(kmc.get_classN()[25])/3000);
=======
		kmc.saveTxt(path,frame,true,true);//save convolved images
		//kmc.saveTxt(path,frame);//same usual ones
		n_threads= ceil(float(kmc.get_classN()[25])/3500);//update number threads based on number of diffusing adatoms (very empirical..)
		if(n_threads>max_threads) n_threads = max_threads;		
>>>>>>> devel
	}
	if(k%(1+n_steps/10)==0 && proc_ID == root_process){
		//std :: cout << floor(float(k)/n_steps)*100 << "% (threads per process= "<< n_threads << " )"<< std :: flush;
		std :: cout  << " || "<< std :: flush;
		
		}
}

//__________________ FINAL MESSAGES ____________________________
if(proc_ID==0){
	


	int * counter = kmc.get_nevents();
	//int * N_class = kmc.get_classN();
	
 	std :: cout << "\n\nDetachment # nn1=0, nn2=0 " << counter[0] << "\tDetachment # nn1= 1,nn2=0 " << counter[1] <<"\tDetachment # nn1= 2,nn2=0 " 
	 << counter[2]<<"\tDetachment # nn1= 3,nn2=0 " << counter[3] << "\tDetachment # nn1= 4,nn2=0 " << counter[4] 
	 << "\nDetachment # nn1=0,nn2=1 " << counter[5] << "\tDetachment # nn1= 1,nn2=1 " << counter[6] <<"\tDetachment # nn1= 2,nn2=1 " 
	 << counter[7]<<"\tDetachment # nn1= 3,nn2=1 " << counter[8] << "\tDetachment # nn1= 4,nn2=1 " << counter[9] 
	 << "\nDetachment # nn1=0,nn2=2 " << counter[10] << "\tDetachment # nn1= 1,nn2=2 " << counter[11] <<"\tDetachment # nn1= 2,nn2=2 " 
	 << counter[12]<<"\tDetachment # nn1= 3,nn2=2 " << counter[13] << "\tDetachment # nn1= 4,nn2=2 " << counter[14]
	 << "\nDetachment # nn1=0,nn2=3 " << counter[15] << "\tDetachment # nn1= 1,nn2=3 " << counter[16] <<"\tDetachment # nn1= 2,nn2=3 " 
	 << counter[17]<<"\tDetachment # nn1= 3,nn2=3 " << counter[18] << "\tDetachment # nn1= 4,nn2=3 " << counter[19] 
	 << "\nDetachment # nn1=0,nn2=4 " << counter[20] << "\tDetachment # nn1= 1,nn2=4 " << counter[21] <<"\tDetachment # nn1= 2,nn2=4 " 
	 << counter[22]<<"\tDetachment # nn1= 3,nn2=4 " << counter[23]
	 <<"\nAttachment # = " << counter[24] << "\t Diffusion # = " << counter[25] ;

	
	end = std :: time(NULL);
	elapsed_time = end-start;
	t2 = clock();
	seconds = ((float)t2-(float)t1)/ CLOCKS_PER_SEC;
	time(&curr_time);
	std :: cout << "\n Elapsed CPU time = " << seconds <<"s \n"<< std:: endl;
	std :: cout << "\n Elapsed time = " << elapsed_time <<"s \n"<< std:: endl;
	std :: cout << "\n End time " << ctime(&curr_time) << "\n";
	std :: cout << "-----------"<<"\n"<< std:: endl;

	//Other final messages like physical time
}

// _______________ FINAL PRINTS ___________________________




// kmc.reset();
// auto path2 = "dummy";
// kmc.saveTxt(path2,0);

kmc.print_final(n_steps/print_every,1);

// double sigma = 1;
// double** convResult_isl, ** convResult_adt; 
// for (int i = 0; i < 20; i++)
// {
	
// 	kmc.initConv_adatom(sigma*3);
// 	kmc.initConv_island(sigma);
// 	convResult_isl = kmc.getIslandConv();
// 	convResult_adt = kmc.getAdatomConv();
	
// 	std:: string file1 = "dummy/island" + std::to_string(i) + ".txt";
// 	std:: string file2 = "dummy/adatom" + std::to_string(i) + ".txt";
// 	printFile(convResult_isl,L,file1,std::to_string(sigma));
// 	printFile(convResult_adt,L,file2,std::to_string(sigma*3));
// 	sigma+=0.2;

// }


MPI_Finalize();




return 0;
 
}



