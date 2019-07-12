/*######################### 2D SINGLE LAYER CLUSTER GROWTH COUPLED TO A CONCENTRATION FIELD ######################

 Author:								Luca Gagliardi IIT
 Copiright:							© 2019 Gagliardi IIT, ILM-CNRS
 
 
###############################################################################################

 ------------ COMPILING INSTRUCTIONS ----------------
place in working folder *.cpp and *.h
compile: mpic++ *.cpp -o "executable".ex 
run:  mpirun -np "cores number" executable.out
input file: file Input.txt containing input informations must be present 

------------------- CONVENTIONS -----------
System of coordinates is x: left-right-wise and y: top-bottom wise. The indexe start with 0.


--------------------- TODO -----------------
- Time computation 
-second nn implementation and new parameter given by ratio J/J'

- program architecture can be improved.. for instance L in both master and dependent classes is redundant.


- Spostare bordello di KMC_step in una funzione dedicata?



- Would be faster if I already individuate contiguos elements between detachment 2 and 3
- Sort the event list for better organization of data?
--------------------


---------- QUESTIONS ------------------

- How to derive c_eq ?



################################################################################
*/



// Useful stuff
/*
INLINE: command make speed up for short function. Indeed this tells to the 
		compiler to expand the function inline when a call is made with significant speed-up
TEMPLATE: to make function or classes versatile on different types
 "" string
 '' character
*/

// GLOBAL VARIABLES

// use stat keyword in case this istance must not be accessible by linked files of the same project


// #############################################
//					MAIN	
// #############################################
#include "global.h"
#include "kmc.h"
#include "functions.h"

// --- PARALLELISATION
#include <mpi.h>
//------


//Define global variables
const double PI = 3.14159265358979323846;
const int root_process=0;
int proc_ID;
int ierr,n_proc;
	
int main(int argc, char **argv){
	

//	const double  J=1, A= 0.1, F = 0.025;//Maybe some from input file?


	

	//int x,y; //just dummy variables fr clarity
	

	

	clock_t t1,t2;
	t1=clock();
	float seconds;
	

	
	ierr = MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	MPI_Comm_rank(MPI_COMM_WORLD, &proc_ID);
	


	const double J =0.2;

	double T0,conc0, A;
	int L,radius, n_steps;
	int frame, print_every;
	
	bool read_old;
	
	if(proc_ID == root_process){
		// read from input file and broadcast	
		std :: cout << "\n Number of processors= "<< n_proc << std :: endl ;
		read_input(&L, &T0,& conc0, &radius, &A, &n_steps, &print_every, &read_old);

		std :: cout << "\n J= " << J << "  L= "<< L<< " T =" << T0 <<" concentration = "<< conc0 <<
		 " radius = "<< radius <<  " 'attachment parameter ' =" << A<<" kmc steps =" << n_steps<<
		" print each =" << print_every << " read old file? =" << read_old<<"\n";
	}

	//RootID broadcast data to other processors
	MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&T0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&conc0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&radius, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&print_every, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&read_old, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

	// make directories for each thread

	auto aString = "mkdir plots" + (std::to_string(proc_ID));
	const char* makeDir = aString.c_str();
	auto path = "plots" + (std::to_string(proc_ID));
	auto bString ="rm " + path + "/island* " + path + "/adatom* ";
	const char* remove_old = bString.c_str();
	system(makeDir);
	system(remove_old);


// 	___________________INITIALIZATION _____________________________

	srand (time(NULL)*(proc_ID+1));// initialise random generator differently for each thread
	KMC kmc(J,A);
	kmc.init(L,radius,conc0,T0, read_old);
	kmc.print(0);


// _________________________RUN KMC ___________________________


frame = 0;

for (int k = 0; k < n_steps; k++){
	// Temperature function..

	kmc.step(T0);

	if ((k%print_every)== 0){
		frame+=1;
		kmc.print(frame);
	}
	if(k%(n_steps/10)==0&&proc_ID == root_process){
		std :: cout  << " | "<< std :: flush;
		}
}

//__________________ FINAL MESSAGES ____________________________
if(proc_ID==0){
	int * counter;
	counter = kmc.get_nevents();
 	std :: cout << "\n\nDetachment # nn0= " << counter[0] << "\tDetachment # nn1= " << counter[1] <<"\tDetachment # nn2= " 
	 << counter[2]<<"\tDetachment # nn3= " << counter[3] << "\tAttachment # = " << counter[4] << "\t Diffusion # = " << counter[5] ;


	std :: cout << "\n Numeber of parallel KMCs ="<< n_proc<<"\n";
	
		
	t2 = clock();
	seconds = ((float)t2-(float)t1)/ CLOCKS_PER_SEC;
	std :: cout << "\n Elapsed time = " << seconds <<"s \n"<< std:: endl;
	std :: cout << "-----------"<<"\n"<< std:: endl;

	//Other final messages like physical time
}

// _______________ FINAL PRINTS ___________________________

	kmc.print_final(n_steps/print_every);

MPI_Finalize();


 return 0;
 
}


