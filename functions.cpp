#include "global.h"
#include "functions.h"



void read_input(int* L, double* T, double* c, int* radius, bool* is_circle, 
double* A, double* BR, double* E_shift, int* n_steps, int* print_every, bool* read_old){

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
        
		*radius= std::stoi(line);//initial size or radius 
		
		std :: getline(finput,line);

		//std::string delimiter = " ";

		//line = line.substr(0, line.find(delimiter)); // token is "scott"

		//std :: cout << "\n " << line[0] << std:: flush;
 		//int answ = std::stoi(line);
		 

		if(line[0] == 'y'){	
			*is_circle= true;// Square or circle? 
		}
		else if (line[0] == 'n'){
			*is_circle = false;
		}
		else{
			std :: cout << "\n\n ABORT, not legal shape inserted\n" <<std:: flush; 
			exit(EXIT_FAILURE);
		}
		
		
		std :: getline(finput,line);
		
		*A= std::stod(line);//attachment parameter
		
		std :: getline(finput,line); 

		*BR= std::stod(line);//bond energy ratio parameter

		std :: getline(finput,line); 

		*E_shift= std::stod(line);//"Excess concentration term"
		
		std :: getline(finput,line); 
		
		*n_steps= std::stoi(line); // kmc steps
		
		std :: getline(finput,line); 
		
		*print_every= std::stoi(line);
		
		// std :: getline(finput,line, ' '); 
		
		// //auto answ = line;
    
        // *read_old = false;
        // if (line == "yes"){
        //     std:: cout << "\n Reading from last configuration. L, and temperature (and concentration) will be overwritten \n";
        //     *read_old = true;
        // }
		// else if (line=="no")
		// {
		// 	//do nothing
		// }
		// else
		// {
		// 	std:: cout << "\n ABORT. Not legal answer";
		// 	exit(EXIT_FAILURE);
		// }
		
		finput.close();
	}
	else {
		std :: cout << "input file not existing. Abort \n";
		exit(EXIT_FAILURE);
	}

	if(((*L)*(*L))>RAND_MAX){
		std :: cout << "L too big . Abort\n";
		exit(EXIT_FAILURE);
	}
	

}

