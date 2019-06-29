/*############ 2D SURFACE EVOLUTION UNDER THE EFFECT OF A LOAD USING KMC ######################

 Author:								Luca Gagliardi IIT
 Copiright:				© 2019 Gagliardi IIT, ILM-CNRS
 
 
###############################################################################################

 ------------ COMPILING INSTRUCTIONS ----------------
place in working folder *.cpp and *.h
compile: mpic++ *.cpp -o "executable".ex 
run:  mpirun -np "cores number" executable.out
input file: file Input.txt containing input informations must be present 

------------------- CONVENTIONS -----------
System of coordinates is x: left-right-wise and y: top-bottom wise. The indexe start with 0.


--------------------- TODO -----------------
- give directly new value of nn without useless obvous computation ?
-change outmut format (directly make a 2d matrix or different columns for each L*L matrix)
- A bit unatisfactory to not have to separate classes for adatom and island.
However if I do so the coupling between the two must happen in a separate routine, but it can be done 
- Spostare bordello di KMC_step in una funzione dedicata?
- Add possibility to read from file
- Would be faster if I already individuate contiguos elements between detachment 2 and 3
- Sort the event list for better organization of data?
--------------------



------------- CONTENT -------------------
* parallelization over different replicas.

---------- QUESTIONS ------------------
+ Interesting to do boolean operations with sort of masks?
+ Shall I always update physical position ? (yes)
+ What happens if 2 or more adatoms one on the other when this is close to an attachment site?
	for now the chances of getting attached are highr (more adatom in the event list of attachment)


-------------------- NEWS--------------------


################################################################################
*/


#include "global.h"
#include "functions.h"
#include "Island.h"
#include "Adatom.h"
#include "Events.h"
#include <vector>

// --- PARALLELISATION
#include <mpi.h>
//------


//Define global variables
const double  J=1, A= 0.5, F = 0.1;//Maybe some from input file?
const double PI = 3.14159265358979323846;
const int n_classes = 5;
int L;

int n_proc,ierr;
int proc_ID;




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


	
int main(int argc, char **argv){
	
	int print_every,n_steps;
	int root_process=0;
	int radius;


	//int x,y; //just dummy variables fr clarity
	double T,conc0;

	int int_steps,frame;

	clock_t t1,t2;
	t1=clock();
	float seconds;
	int counter [n_classes];
	const bool debug_mode = false;
	
// MPI stuff ------ 
	/* Now replicate this process to create parallel processes.
	 * From this point on, every process executes a seperate copy
	 * of this program */
	
	ierr = MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	MPI_Comm_rank(MPI_COMM_WORLD, &proc_ID);
	
	if(proc_ID == root_process){
		std :: cout << "\n Number of processors= "<< n_proc << std :: endl ;
	}
// ----------------
	

	
    
		
		//read_input(&L,&T, &int_steps);
		
		std::string aString = "mkdir plots" + (std::to_string(proc_ID));
		const char* makeDir = aString.c_str();
		std::string path = "plots" + (std::to_string(proc_ID));
		std::string bString ="rm " + path + "/island* " + path + "/adatom* ";
		const char* remove_old = bString.c_str();
		system(makeDir);
		system(remove_old);
		
		

	
	
	
	//RootID broadcast data to other processors
	
	//MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&F, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&int_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&compute_every, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&N_sim, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&pattern_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	
	// ++++++++++++ Class declarations ++++++++++++++++
	//std :: cout <<"M="<< M << "\n"<<std :: endl;


	//M=int_steps/compute_every;
	
	
	srand (time(NULL)*(proc_ID+1));// initialise random generator differently for each thread
	

int i_rand =i_rand = rand() % 4 +1;
std:: cout << i_rand;

	L = 40;
	n_steps = 10000;
	print_every = 100;
	conc0 = 0.6;
	radius = 4;
	T =0.1;


// INITIALIZATION -----------------
	Island island(radius);
	Adatom adatom(conc0);

	Events R[5] = {Events(L),Events(L),Events(L),Events(L),Events(L)};
	R[0].D = energy(1,T,J);
	R[1].D = energy(2,T,J);
	R[2].D = energy(3,T,J);
	R[3].D = A;
	R[4].D = 4.0*F; //4 possible movements
	island.init_neighbours(); 

// Fill classes
	for (int y = 0; y < L; y++){
		for (int x = 0; x < L; x++){		
			if (island.nn[y][x]==1){
				R[0].populate(x,y);
			}
			if (island.nn[y][x]==2){
				R[1].populate(x,y);
			}
			if (island.nn[y][x]==3){
				R[2].populate(x,y);
			}
			
			for(int k=0; k<adatom.matrix[y][x]; k++){
				// I consider possibility of multiple adatoms on top of each other 
				R[4].populate(x,y);
				if(is_attSite(x,y,island)) R[3].populate(x,y);
			}
			}
		}
	
// std :: cout << "\nTRYING EXIST\n" ;
// int x =rand() %L;
// int y =rand() %L;

// std :: cout << "\n coordinate " << x <<" ,  "<< y << ":\t"<<R[4].exist(x,y) ;


// -----------------------------------------


	frame = 0;

	//save in parallel images
	
	std :: ofstream outfile (path+"/configuration.txt");
	if (outfile.is_open()){
		outfile << L << "\t" << (n_steps/print_every);
	}
	outfile.close();

	island.print(path + "/island" + std::to_string(frame) + ".txt" ,R[0].mask, R[1].mask, R[2].mask,T);
	adatom.print(path + "/adatom" + std::to_string(frame)+ ".txt",R[3].mask,T);

	//print configuration file..




// CHECKS ------------------------------
if(proc_ID==0){
	// std :: cout<< "\n island \n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< island.matrix[i][j]; 
	// 	}
	// 	std :: cout <<"\n";
	// }
	// std :: cout<< "\n adatoms \n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< adatom.matrix[i][j]; 
	// 	}
	// 	std :: cout <<"\n";
	// }
	std :: cout << "\n total number adatoms"<< (adatom.N) ;
	std :: cout << "\n Elements in det class 1:  " << R[0].N << "\t elements in det class 2:  " << R[1].N << "\t elements in det class 3:  " << R[2].N  << "\n";
 	std :: cout << "\n Elements in att class:" << R[3].N << "\t elements in diffusion class" << R[4].N << "\n";

	// for (int i = 0; i < R[3].N; i++)
 	//  {
	// 	std :: cout << i <<"\t(" << R[4].where(i)[0]<< ","<< R[4].where(i)[1] << ")\n";
 	// }
	// std :: cout<< "\n diffusion mask \n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< R[4].mask[i][j]; 
	// 	}
	// std :: cout <<"\n";
	// }

	// for (int i = 0; i < R[3].N; i++)
 	//  {
	// 	std :: cout << i <<"\t(" << R[3].where(i)[0]<< ","<< R[3].where(i)[1] << ")\n";
 	// }
	// std :: cout<< "\n attachment mask\n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< R[3].mask[i][j]; 
	// 	}
	// std :: cout <<"\n";
	// }


}




// 	std :: cout<< "\n neighbours \n";
// 	for ( int i = 0;i <L;i++){
// 		for(int j =0;j<L;j++){
// 			std::cout << "\t"<< island.nn[i][j]; 
// 			}
// 		std :: cout <<"\n";
// 		}

	
// std :: cout << "\n Elements in the classes:" << R[0].N << ", \t" << R[1].N << ", \t" << R[2].N << "\n";
// for (int i = 0; i < R[1].N; i++)
//  {
// 	std :: cout << i <<"\t(" << R[1].where(i)[0]<< ","<< R[1].where(i)[1] << ")\n";
//  }
//  std :: cout << "\n \n";
// for (int i = 0; i < R[2].N; i++)
//  {
// 	std :: cout << i <<"\t(" << R[2].where(i)[0]<< ","<< R[2].where(i)[1] << ")\n";
//  }

//---------------------------


// RUN KMC -----------------------------

T =1 ;

for (int k = 0; k < n_steps; k++)
{
	KMC_step(R, island, adatom, T, counter);

	if ((k%print_every)== 0){
		frame+=1;
		island.print(path + "/island" + std::to_string(frame) + ".txt" ,R[0].mask, R[1].mask, R[2].mask,T);
		adatom.print(path + "/adatom" + std::to_string(frame)+ ".txt",R[3].mask,T);
	}
	

if(proc_ID==0 &&debug_mode){
	//CHECKS*********************
	std :: cout<< "\n  " << k+1 <<"  KMC step \n";

		std :: cout<< "\n island \n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< island.matrix[i][j]; 
		}
		std :: cout <<"\n";
	}

	std :: cout<< "\n neighbours \n";

	for ( int i = 0;i <L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< island.nn[i][j]; 
			}
		std :: cout <<"\n";
		}

	std :: cout<< "\n adatoms \n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< adatom.matrix[i][j]; 
		}
		std :: cout <<"\n";
	}
 	std :: cout << "\n total number adatoms"<< (adatom.N) ;
	std :: cout << "\n Elements in det class 1:  " << R[0].N << "\t elements in det class 2:  " << R[1].N << "\t elements in det class 3:  " << R[2].N  << "\n";
	std :: cout << "\n Elements in att class:" << R[3].N << "\t  elements in diffusion class" << R[4].N << "\n";

// 	std :: cout << "\n detachment 1nn list \n";
// 	 for (int i = 0; i < R[0].N; i++)
//  	  {
// 	 	std :: cout << i <<"\t(" << R[0].where(i)[0]<< ","<< R[0].where(i)[1] << ")\n";
//  	 }

// 	std :: cout << "\n detachment 2nn list \n";
// 	 for (int i = 0; i < R[1].N; i++)
//  	  {
// 	 	std :: cout << i <<"\t(" << R[1].where(i)[0]<< ","<< R[1].where(i)[1] << ")\n";
//  	 }

// 	std :: cout << "\n detachment 3nn list \n";
// 	 for (int i = 0; i < R[2].N; i++)
//  	  {
// 	 	std :: cout << i <<"\t(" << R[2].where(i)[0]<< ","<< R[2].where(i)[1] << ")\n";
//  	 }

	// std :: cout << "\n attachment list \n";
	//  for (int i = 0; i < R[3].N; i++)
 	//   {
	//  	std :: cout << i <<"\t(" << R[3].where(i)[0]<< ","<< R[3].where(i)[1] << ")\n";
 	//  }

	// std :: cout << "\n diffusion list \n";
	// for (int i = 0; i < R[4].N; i++)
 	//   {
	//  	std :: cout << i <<"\t(" << R[4].where(i)[0]<< ","<< R[4].where(i)[1] << ")\n";
 	//  }


	std :: cout<< "\n det 1 mask \n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< R[0].mask[i][j]; 
		}
	std :: cout <<"\n";
	}

std :: cout<< "\n det 2 mask \n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< R[1].mask[i][j]; 
		}
	std :: cout <<"\n";
	}

	std :: cout<< "\n det 3 mask \n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< R[2].mask[i][j]; 
		}
	std :: cout <<"\n";
	}


	std :: cout<< "\n attachment mask\n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< R[3].mask[i][j]; 
		}
	std :: cout <<"\n";
	}
}
}
if(proc_ID==0){
 	std :: cout << "\tDetachment # nn1= " << counter[0] << "\tDetachment # nn2= " << counter[1] <<"\tDetachment # nn3= " << counter[2] << "\tAttachment # = " << counter[3] << "\t Diffusion # = " << counter[4] ;




	std :: cout << "\n Numeber of parallel KMCs ="<< n_proc<<"\n";
		//output(av_obs,n_proc,"observables_F"+ std::to_string(F)+"_T"+std::to_string(T)+".txt",N_sim);
		
		//std :: cout << "\n  Averaged final physical time = "<< av_obs.time[M-1]/double(n_proc)<< std:: endl;
		
		//elapsed time for master (which is the last to finish since it does more stuff)
		
	t2 = clock();
	seconds = ((float)t2-(float)t1)/ CLOCKS_PER_SEC;
	std :: cout << "\n Elapsed time = " << seconds <<"s \n"<< std:: endl;
	std :: cout << "-----------"<<"\n"<< std:: endl;
}

//Consistency CHECK
// int counter =0;
// 	for (int i = 0; i < L; i++){
//         for (int j = 0; j < L; j++){
//             counter += adatom.matrix[i][j];
//         }
//     }
// 	if(counter != adatom.N){
// 		std :: cout << "PROBLEM! Inconsistency 1";
// 	}
// MORE..?


MPI_Finalize();


 return 0;
 
}



