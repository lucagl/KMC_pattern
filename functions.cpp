
#include "global.h"
#include "functions.h"
#include "Adatom.h"
#include "Island.h"
#include "Events.h"


// Global variables have to be redeclared somewhere


// void read_input(int* L, double* T, double* F, int* pattern_type,int* N_add, int* N_integration, int* Number_sim, int* diffusionOnly){

// 	double sigma,E;
// 	std :: string line,dummy;
// 	std :: ifstream finput("Input.txt");
	
// 	if (finput.is_open()){
// 		std :: getline(finput,line, '\t'); //extracts until character delimiter '\t'
// 		std :: getline(finput,dummy, '\n');
// 		*L= std::stoi(line); //parse to integer
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		*T= std::stod(line);//parse to double
		
// 		std :: getline(finput,line, '\t');
// 		std :: getline(finput,dummy, '\n');
		
// 		*F= std::stod(line); //deposition rate (constant)
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		load= std::stod(line);
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		*pattern_type= std::stoi(line);
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		*N_add= std::stoi(line);
		
// 		if (*N_add>=((*L)*(*L))) {
// 			std :: cout << "Number of additives must be lower than total number. Abort\n";
// 			exit(EXIT_FAILURE);
// 		}
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		*N_integration= std::stoi(line);
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		compute_every= std::stoi(line);
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		print_surface= std::stoi(line);//used only by root_ID
		
// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		*Number_sim = std::stoi(line);

// 		std :: getline(finput,line, '\t'); 
// 		std :: getline(finput,dummy, '\n');
		
// 		*diffusionOnly = std::stoi(line);
// 	}
// 	else {
// 		std :: cout << "input file not existing. Abort\n";
// 		exit(EXIT_FAILURE);
// 	}
	
// 	sigma = 1;
// 	E=1;	
// 	A =(1+sigma)*(2*sigma-1)*sigma/(4*PI*PI*E);	
	
// 	std :: cout << "Size of the box is "<< *L<< " X "<< *L << "\n";
// 	std :: cout << "Temperature " << *T << "\n";
// 	std :: cout << "Density of impurities " << double(*N_add)/double((*L)*(*L))<< "\n";
// 	//std :: cout << "Interaction strenght " << J << "\n";
// 	if (!(*pattern_type==0)) {
// 	std :: cout << "Pattern type " << *pattern_type << "\n";
// 	std :: cout << "Load " << load << "\n";
// 	std :: cout << "Elastic strenght " << A << "\n";
// 	}
// 	else {std :: cout << "No indent, Load =0 " << "\n"; }
// 	//std :: cout << "Young modulus " << E << "\n";
// 	//std :: cout << "Poisson ratio " << sigma << "\n";
// 	if (*diffusionOnly ==1)
// 	{
// 	std :: cout << "Only surface diffusion \n";
// 	}
// 	else{
// 	std :: cout << "Deposition rate " << *F << "\n";
// 	}

// 	std :: cout << "Number of KMC steps " << *N_integration << "\n";
// 	std :: cout << "Frequency observable update " << compute_every << "\n";
// 	std :: cout << "Frequency surface plot " << print_surface << "\n";
// 	std :: cout << "Simulation per core " << *Number_sim << "\n";
// 	if(((*L)*(*L))>RAND_MAX){
// 		std :: cout << "L too big . Abort\n";
// 		exit(EXIT_FAILURE);
// 	}
	

// }



// pass by reference the object because potentially big (if I pass by value it will generate a potentially big copy)

double energy(int nn, double T, double J){

    double r;
    r = exp(-J*nn/T);

    return r;
}

double cumulative (Events R[], double* r){
    
    double R_sum =0;
    r[0] = 0;
    for (int i = 1; i <=n_classes; i++)
    {
        r[i] = r[i-1] + R[i-1].rate();
    }
    R_sum = r[n_classes];
    return R_sum;
}

int extract (int N){

    return (rand() % N);
}









int is_attSite(int x,int y,const Island &island){
	
	int contact = 0 ;
	int top,bottom,left,right;


	top = y+1;
	if(top==L) top = 0;

	bottom = y-1;
	if(bottom==-1) bottom = L-1;

	right = x+1;
	if (right ==L) right = 0;
	
	left = x -1;
	if(left == -1) left = L-1; 

	contact = (island.matrix[bottom][x] | island.matrix[top][x] | island.matrix[y][left] | island.matrix[y][right])^island.matrix[y][x];
	
	return contact;
}

// Counts elements in LxL matrix (not necessarly columns only 1)

// int non_zero(const int M[L][L]){
//     int counter =0;

//     for (int i = 0; i < L; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             counter += M[i][j];
//         }
        
//     }
//     return counter;
// }



 void KMC_step(Events R[], Island &island, Adatom &adatom , double T, int* event_counter){

    static int count_dnn1=0,count_dnn2=0,count_dnn3=0,count_a=0,count_d=0;
   
	int index,i_rand;
	double R_sum,d_rand;
	int left,right,top,bottom,x,y;
    double r [n_classes+1];

    event_counter[0] = count_dnn1;
    event_counter[1] = count_dnn2;
    event_counter[2] = count_dnn3;
    event_counter[3] = count_a;
    event_counter[4] = count_d;

    
   
    


 	R_sum = cumulative(R,r); 

    std :: cout << "\n \n \n \n "<< r[0] <<  "\t" << r[1] << "\t"<< r[2] << "\t"<< r[3] << "\t"<< r[4] << "\t"<< r[5] << "\n";


 	d_rand =  ((double) rand() / (RAND_MAX)) * R_sum;

    std :: cout << "\n \n  " <<d_rand;

// DETACHMENT AT 1 NN SITE
 	if (r[0]<d_rand && d_rand<r[1]){

        count_dnn1+=1;
		index = extract(R[0].N);//simple uniform random generator
		// locate event
	     std :: cout << "\n detachment event 1 , index="<< index;
        // std:: cout << "\n counter 1    " << count_dnn1 << "\n";
 		x = R[0].where(index)[0];
		y = R[0].where(index)[1];

		//update with flag 1

		//update(R,1);

		R[0].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;

		//update R3 and  neighbours
		island.nn[y][x] = 0;
		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 
		

		if(island.matrix[y][right]){
			 island.nn[y][right] -= 1;
            if(island.nn[y][right] ==1){
			 	R[0].populate(right,y); 
				R[1].destroy_byPosition(right,y); 
			 }
			 else if(island.nn[y][right] ==2){
			 	R[1].populate(right,y); 
				R[2].destroy_byPosition(right,y); 
			 }
			 else if(island.nn[y][right] ==3){
				 R[2].populate(right,y); 
			 }
			 else {
				 std :: cout << "Error" ;
				 exit(EXIT_FAILURE);
			 }
		}
		else if(island.matrix[y][left]) {
			island.nn[y][left] -= 1; 
            if(island.nn[y][left] ==1){
			 	R[0].populate(left,y); 
				R[1].destroy_byPosition(left,y); 
			 }
			else if(island.nn[y][left] ==2){
			 	R[1].populate(left,y); 
				R[2].destroy_byPosition(left,y); 
			 }
			else if(island.nn[y][left] ==3){
				 R[2].populate(left,y); 
			 }
			else {
				std :: cout  << "Error" ;
				 exit(EXIT_FAILURE);
			 }
			}
		else if(island.matrix[top][x]) {
			island.nn[top][x] -= 1; 
            if(island.nn[top][x] ==1){
			 	R[0].populate(x,top); 
				R[1].destroy_byPosition(x,top); 
			 }
			else if(island.nn[top][x] ==2){
			 	R[1].populate(x,top); 
				R[2].destroy_byPosition(x,top); 
			 }
			else if(island.nn[top][x] ==3){
				 R[2].populate(x,top); 
			 }
			else {
				std :: cout  << "Error" ;
				 exit(EXIT_FAILURE);
			 }
			}
		else if(island.matrix[bottom][x]){
			island.nn[bottom][x] -= 1;
            if(island.nn[bottom][x] ==1){
			 	R[0].populate(x,bottom); 
				R[1].destroy_byPosition(x,bottom); 
			 } 
			else if(island.nn[bottom][x] ==2){
			 	R[1].populate(x,bottom); 
				R[2].destroy_byPosition(x,bottom); 
			 }
			else if(island.nn[bottom][x] ==3){
				 R[2].populate(x,bottom); 
			 }
			else {
				std :: cout  << "Error" ;
				 exit(EXIT_FAILURE);
			 }
			}

		//Attachment list update
        R[3].populate(x,y);
        //Diffusion list upadate
        R[4].populate(x,y);

	}


// DETACHMENT AT 2 NN SITE
 	else if (r[1]<d_rand && d_rand<r[2]){
		count_dnn2+=1;

		index = extract(R[1].N);
		std :: cout << "\n detachment event 2 , index="<< index;
		// std:: cout << "\n counter 2    " << count_dnn2 << "\n";


		x = R[1].where(index)[0];
		y = R[1].where(index)[1];
		R[1].destroy(index);
		//update R3, R2 and R1 and nn
	

		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;

		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 

		

		if (island.nn[y][left] == 4)
		{
			 R[2].populate(left,y);
		}
		
		else if(island.nn[y][left] ==3){
			 R[1].populate(left,y);
			 R[2].destroy_byPosition(left,y);
		}
		else if(island.nn[y][left] ==2){
			R[0].populate(left,y);
			R[1].destroy_byPosition(left,y);
		}

		if (island.nn[y][right] == 4)
		{
			 R[2].populate(right,y);
		}
		else if(island.nn[y][right] ==3){
			R[1].populate(right,y);
			R[2].destroy_byPosition(right,y);
		}
		else if(island.nn[y][right] ==2){
			R[0].populate(right,y);
			R[1].destroy_byPosition(right,y);
		}

		if(island.nn[top][x] ==4){
			R[2].populate(x,top);
		}
		else if(island.nn[top][x] ==3){
			R[1].populate(x,top);
			R[2].destroy_byPosition(x,top);
		}
		else if(island.nn[top][x] ==2){
				R[0].populate(x,top);
				R[1].destroy_byPosition(x,top);
		}

		if(island.nn[bottom][x] ==4){
			R[2].populate(x,bottom);
		}
		else if(island.nn[bottom][x] ==3){
				R[1].populate(x,bottom);
				R[2].destroy_byPosition(x,bottom);
		}
		else if(island.nn[bottom][x] ==2){
			R[0].populate(x,bottom);
			R[1].destroy_byPosition(x,bottom);
		}
		
		
		
	// Attachment list update 

        R[3].populate(x,y);
     //Diffusion list update
        R[4].populate(x,y);


	//Many possible combinations, this is why no specific update above*/

		island.nn[y][x] = 0;
		island.nn[top][x] = island.get_neighbours(x,top); 
		island.nn[bottom][x] = island.get_neighbours(x,bottom); 

		island.nn[y][left] =island.get_neighbours(left,y); 
		island.nn[y][right] = island.get_neighbours(right,y); 

		island.nn[bottom][left] = island.get_neighbours(left,bottom);
		island.nn[bottom][right] =  island.get_neighbours(right,bottom);
		island.nn[top][left] = island.get_neighbours(left,top);
		island.nn[top][right] =  island.get_neighbours(right,top);


	
	}



    // DETACHMENT AT 3 NN SITE
	else if (r[2]<d_rand && d_rand<r[3]){
		index = extract(R[2].N);
		std :: cout << "\n detachment event 3 , index ="<< index;
         count_dnn3+=1;

		x = R[2].where(index)[0];
		y = R[2].where(index)[1];
		R[2].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;


		// update R2 and R3 and nn

		island.nn[y][x] =0;

		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 

		

		if (!(island.matrix[top][x])){
			R[1].populate(left,y);
			R[1].populate(right,y);
			island.nn[y][left] = 2;
			island.nn[y][right] = 2;
			island.nn[bottom][x] =3;//a 3 neighbour site appears below detached element
			R[2].populate(x,bottom);
		}
		else if (!(island.matrix[bottom][x])){
			R[1].populate(left,y);
			R[1].populate(right,y);
			island.nn[y][left] = 2;
			island.nn[y][right] = 2;
			island.nn[top][x] =3;
			R[2].populate(x,top);
		}
		else if (!(island.matrix[y][right])){
			R[1].populate(x,top);
			R[1].populate(x,bottom);
			island.nn[top][x] = 2;
			island.nn[bottom][x] = 2;
			island.nn[y][left] = 3;
			R[2].populate(left,y);
		}
		else if (!(island.matrix[y][left])){
			R[1].populate(x,top);
			R[1].populate(x,bottom);
			island.nn[top][x] = 2;
			island.nn[bottom][x] = 2;
			island.nn[y][right] = 3;
			R[2].populate(right,y);
		}

 		//ATTACHMENT list update

        R[3].populate(x,y);

        //Diffusion list update

        R[4].populate(x,y);
	}

    // GET ATTACHED
	else if (r[3]<d_rand && d_rand<r[4]){

        count_a+=1;
        index = extract(R[3].N);
        std :: cout << "\n attachment event, index="<< index;
        // std:: cout << "\n counter att    " << count_a << "\n";
        

        x = R[3].where(index)[0];
		y = R[3].where(index)[1];
		R[3].destroy(index);
        R[4].destroy_byPosition(x,y);

        adatom.matrix[y][x]-=1;
        adatom.N -=1.0;
        island.matrix[y][x] = 1;


        // Update island classes involed and nn

        

        top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 

        //R[3]: new island site might imply new attachment site for neighbouring adatoms

        if(adatom.matrix[y][right]>=0 && R[3].exist(right,y)) R[3].populate(right,y);//Problem: I have to establish if that site was already filled..
        if(adatom.matrix[y][left]>=0 && R[3].exist(left,y)) R[3].populate(left,y);
        if(adatom.matrix[top][x]>=0 && R[3].exist(x,top)) R[3].populate(x,top);
        if(adatom.matrix[bottom][x]>=0 && R[3].exist(x,bottom)) R[3].populate(x,bottom);


        //----

        //on island:
        int  s = 0;
        if(island.matrix[y][right]){
            s+=1;
            if(island.nn[y][right]==3){
                R[2].destroy_byPosition(right,y);
            }
            else if(island.nn[y][right]==2){
                R[2].populate(right,y);
                R[1].destroy_byPosition(right,y);
            }
            else if(island.nn[y][right]==1){
                R[1].populate(right,y);
                R[0].destroy_byPosition(right,y);
            }
            island.nn[y][x] +=1;
            island.nn[y][right] +=1;
        }
        if(island.matrix[y][left]){
            s+=1;
            if(island.nn[y][left]==3){
                R[2].destroy_byPosition(left,y);
            }
            else if(island.nn[y][left]==2){
                R[2].populate(left,y);
                R[1].destroy_byPosition(left,y);
            }
            else if(island.nn[y][left]==1){
                R[1].populate(left,y);
                R[0].destroy_byPosition(left,y);
            }
            island.nn[y][x] +=1;
            island.nn[y][left] +=1;
        }
        if(island.matrix[top][x]){
            s+=1;
            if(island.nn[top][x]==3){
                R[2].destroy_byPosition(x,top);
            }
            else if(island.nn[top][x]==2){
                R[2].populate(x,top);
                R[1].destroy_byPosition(x,top);
            }
            else if(island.nn[top][x]==1){
                R[1].populate(x,top);
                R[0].destroy_byPosition(x,top);
            }
            island.nn[y][x] +=1;
            island.nn[top][x] +=1;
        }
        if(island.matrix[bottom][x]){
            s+=1;
            if(island.nn[bottom][x]==3){
                R[2].destroy_byPosition(x,bottom);
            }
            else if(island.nn[bottom][x]==2){
                R[2].populate(x,bottom);
                R[1].destroy_byPosition(x,bottom);
            }
            else if(island.nn[bottom][x]==1){
                R[1].populate(x,bottom);
                R[0].destroy_byPosition(x,bottom);
            }
            island.nn[bottom][x] +=1;
            island.nn[y][x] +=1;
        }
        R[s-1].populate(x,y);
        if(s>=3){
            std :: cout << "Error in attachment. Badly recognized environment";
            exit (EXIT_FAILURE);
        }
    
     }

// Simple diffusion event
	else if (r[4]<d_rand && d_rand<r[5]){

        count_d+=1;
        index = extract(R[4].N);

        std :: cout << "\n *HERE*\n";
        std :: cout << "\n diffusion event , index="<< index;
        // std:: cout << "\n counter diff    " << count_d << "\n";
        


        x = R[4].where(index)[0];
		y = R[4].where(index)[1];

        R[4].destroy(index);
        adatom.matrix[y][x] -=1;
        i_rand = rand() % 4 +1;

        if(is_attSite(x,y,island)){
        //remove if it was (in the previous position) on an attached site
            R[3].destroy_byPosition(x,y);
        }
        
        if(i_rand ==1){

            top = y+1;
            if(top==L) top = 0;
            R[4].populate(x,top);
            adatom.matrix[top][x] += 1;

            if(is_attSite(x,top,island)){
                R[3].populate(x,top);
            }
        }
        else if(i_rand ==2){

		    bottom = y-1;
		    if(bottom==-1) bottom = L-1;
            R[4].populate(x,bottom);
            adatom.matrix[bottom][x] += 1;
            
            if(is_attSite(x,bottom,island)){
                R[3].populate(x,bottom);
            }
       }
        else if(i_rand ==3){

		    right = x+1;
		    if (right ==L) right = 0;

            R[4].populate(right,y);
            adatom.matrix[y][right] += 1;

            if(is_attSite(right,y,island)){
                R[3].populate(right,y);
            }
       }
        else if(i_rand ==4){

		    left = x -1;	
		    if(left == -1) left = L-1; 

            R[4].populate(left,y);
            adatom.matrix[y][left] += 1;

            if(is_attSite(left,y,island)){
                R[3].populate(left,y);
            }
       }

    }


}