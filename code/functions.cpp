#include "global.h"
#include "functions.h"



#include <memory>


void read_input(std :: string inName,int* nthreads, int* L, double* T, double* c, int* radius, bool* is_circle, 
double* A, double* BR, double* E_shift, unsigned long* n_steps, int* print_every){

	std :: ifstream finput;

	if(inName!=""){
		finput.open(inName.c_str());
		//std :: cout << "Reading " << inName.c_str()<< "\n";
	}
	else {finput.open("Input.txt");}


    std :: string line;
	

	if (finput.is_open()){

		//NEW get n threads from user
		std :: getline(finput,line);
		*nthreads = std::stoi(line);


		std :: getline(finput,line); 
		
		*L= std::stoi(line); //parse to integer
		
		std :: getline(finput,line); 

		*T= std::stod(line);//initial temperature
        
        std :: getline(finput,line); 
		
		*c= std::stod(line);//initial concentration

        std :: getline(finput,line); 
        
		*radius= std::stoi(line);//initial size or radius 
		
		std :: getline(finput,line);

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
		
		*n_steps= std::stol(line); // kmc steps
		
		std :: getline(finput,line); // frames
		
		*print_every= std::stoi(line);


		// std :: getline(finput,line); //sigma
		// *sigma = std::stod(line);
		
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
		std :: cout << "input file not provided. Abort \n";
		exit(EXIT_FAILURE);
	}

	if(((*L)*(*L))>RAND_MAX){
		std :: cout << "L too big . Abort\n";
		exit(EXIT_FAILURE);
	}
	

}

double* gauss(const int L, const double sigma ){
	
	double mu =double(L)/2;
	double x,y;
	double * g = (double*) malloc (sizeof(double) * L*L);
 
	for (int i = 0; i < L; i++)
	{	
		for (int j = 0; j < L; j++)
		{	x = double(i)-mu;
			y = double(j) - mu;
			g[j+L*i]= 1/(sigma*sigma*2*PI) * exp(-(x*x+y*y)/(2*sigma*sigma));
		}
		
	}

// CHECK
	// std :: string file_name = "dummy/gaussianK.txt";
	// std :: ofstream outfile (file_name);
    //     if (outfile.is_open()){
    //         for ( int i = 0;i < L;i++){
    //             for(int j =0;j<L;j++){
    //                 outfile << g[j+L*i]<< "\n"; 
    //             }
    //         }
    //     }
    //     outfile.close();

//need to shift for fourier space to be correctly centered in frequency.. don't really understand why..

return g;

}

void fftShift(double* in, double* out,const int xdim, const int ydim, int x_shift, int y_shift){

	if(!x_shift && !y_shift){
		x_shift = xdim/2;
		y_shift = ydim/2;
	}

	
  for (int i = 0; i < xdim; i++) {
    int ii = (i + x_shift) % xdim;
    for (int j = 0; j < ydim; j++) {
      int jj = (j + y_shift) % ydim;
      out[ii * ydim + jj] = in[i * ydim + j];
    }
  }

}

void printFile(double** in, const int L, const std::string& file_name,const std::string&flag )  {
	std :: ofstream outfile (file_name);
	if (outfile.is_open()){
		if(flag!="Null") outfile <<flag << "\n";
		for ( int i = 0;i < L;i++){
			for(int j =0;j<L;j++){
				outfile << in[i][j]<< "\n"; 
			}
		}
	}
	outfile.close();
}


std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}