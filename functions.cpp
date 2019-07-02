void read_input(int* L, double* T, double* F, int* pattern_type,int* N_add, int* N_integration, int* Number_sim, int* diffusionOnly){

	double sigma,E;
	std :: string line,dummy;
	std :: ifstream finput("Input.txt");
	
	if (finput.is_open()){
		std :: getline(finput,line, '\t'); //extracts until character delimiter '\t'
		std :: getline(finput,dummy, '\n');
		*L= std::stoi(line); //parse to integer
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		*T= std::stod(line);//parse to double
		
		std :: getline(finput,line, '\t');
		std :: getline(finput,dummy, '\n');
		
		*F= std::stod(line); //deposition rate (constant)
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		load= std::stod(line);
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		*pattern_type= std::stoi(line);
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		*N_add= std::stoi(line);
		
		if (*N_add>=((*L)*(*L))) {
			std :: cout << "Number of additives must be lower than total number. Abort\n";
			exit(EXIT_FAILURE);
		}
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		*N_integration= std::stoi(line);
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		compute_every= std::stoi(line);
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		print_surface= std::stoi(line);//used only by root_ID
		
		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		*Number_sim = std::stoi(line);

		std :: getline(finput,line, '\t'); 
		std :: getline(finput,dummy, '\n');
		
		*diffusionOnly = std::stoi(line);
	}
	else {
		std :: cout << "input file not existing. Abort\n";
		exit(EXIT_FAILURE);
	}
	
	sigma = 1;
	E=1;	
	A =(1+sigma)*(2*sigma-1)*sigma/(4*PI*PI*E);	
	
	std :: cout << "Size of the box is "<< *L<< " X "<< *L << "\n";
	std :: cout << "Temperature " << *T << "\n";
	std :: cout << "Density of impurities " << double(*N_add)/double((*L)*(*L))<< "\n";
	//std :: cout << "Interaction strenght " << J << "\n";
	if (!(*pattern_type==0)) {
	std :: cout << "Pattern type " << *pattern_type << "\n";
	std :: cout << "Load " << load << "\n";
	std :: cout << "Elastic strenght " << A << "\n";
	}
	else {std :: cout << "No indent, Load =0 " << "\n"; }
	//std :: cout << "Young modulus " << E << "\n";
	//std :: cout << "Poisson ratio " << sigma << "\n";
	if (*diffusionOnly ==1)
	{
	std :: cout << "Only surface diffusion \n";
	}
	else{
	std :: cout << "Deposition rate " << *F << "\n";
	}

	std :: cout << "Number of KMC steps " << *N_integration << "\n";
	std :: cout << "Frequency observable update " << compute_every << "\n";
	std :: cout << "Frequency surface plot " << print_surface << "\n";
	std :: cout << "Simulation per core " << *Number_sim << "\n";
	if(((*L)*(*L))>RAND_MAX){
		std :: cout << "L too big . Abort\n";
		exit(EXIT_FAILURE);
	}
	

}
