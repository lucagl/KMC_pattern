#include "global.h"
#include "functions.h"

void read_input(int* L, double* T, double* c, int* radius, double* F, double* A, int* n_steps, int* print_every, bool* read_old){

	std :: ifstream finput("Input.txt");
    std :: string line;
	

	if (finput.is_open()){
		std :: getline(finput,line); 
		
		*L= std::stoi(line); //parse to integer
		

		std :: getline(finput,line); 

		*T= std::stod(line);//initial temperature
        
        std :: getline(finput,line); 
		
		*c= std::stod(line);//initial concentration

        std :: getline(finput,line); 
        
		*radius= std::stoi(line);//initial concentration
		
		std :: getline(finput,line);
		
		*F= std::stod(line); //diffusion rate (constant)
		
		std :: getline(finput,line); 
		
		*A= std::stod(line);//attachment rate
		
		std :: getline(finput,line); 
		
		*n_steps= std::stoi(line); // kmc steps
		
		std :: getline(finput,line); 
		
		*print_every= std::stoi(line);
		
		std :: getline(finput,line, ' '); 
		
		auto answ = line;
    
        *read_old = false;
        if (answ == "yes"){
            std:: cout << "\n Reading from last configuration. L, and temperature (and concentration) will be overwritten \n";
            *read_old = true;
        }
		finput.close();
	}
	else {
		std :: cout << "input file not existing. Abort \n";
		exit(EXIT_FAILURE);
	}
	
	
	// std :: cout << "Size of the box is "<< *L<< " X "<< *L << "\n";
	// std :: cout << "Temperature " << *T << "\n";
	// std :: cout << "Density of adatoms " << << "\n";
    // std :: cout << "Attachment rate " <<  << "\n";
    // std :: cout << "Diffusion rate (per direction) " <<  << "\n";
	
	

	// std :: cout << "Number of KMC steps " << *N_integration << "\n";
	// std :: cout << "Frequency plot " << print_every << "\n";


	if(((*L)*(*L))>RAND_MAX){
		std :: cout << "L too big . Abort\n";
		exit(EXIT_FAILURE);
	}
	

}
