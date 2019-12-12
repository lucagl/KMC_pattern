#include "global.h"
#include "functions.h"



void read_input(int* L, double* T, double* c, int* radius, bool* is_circle, 
double* A, double* BR, double* E_shift, int* n_steps, int* print_every){

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
//subroutine
void gauss(double *g,const int L, const double sigma ){
	
	double mu = double(L)/2;
	double x,y;
 
	for (int i = 0; i < L; i++)
	{	
		for (int j = 0; j < L; j++)
		{	x = double(i)-mu;
			y = double(j) - mu;
			g[j+L*i]= 1/(sigma*2*PI) * exp(-(x*x+y*y)/(2*sigma*sigma));
		}
		
	}

// CHECK
	std :: string file_name = "dummy/gaussianK.txt";
	std :: ofstream outfile (file_name);
        if (outfile.is_open()){
            for ( int i = 0;i < L;i++){
                for(int j =0;j<L;j++){
                    outfile << g[j+L*i]<< "\n"; 
                }
            }
        }
        outfile.close();


	return;
	

}