
#include "global.h"
#include "kmc.h"
#include "Adatom.h"
#include "Island.h"
#include "Events.h"




KMC:: KMC(const double J_read, const double A_read){

        if(proc_ID == root_process){
            std :: cout << "\n Starting KMC with " << n_classes << " classes of events \n";
        }
        J = J_read; // link strenght
        A =A_read; // attachment over diffusion parameter >0 =few attachement, <0 many attachement (diffusion dominated)
       // F =F_read; // diffusion constant

// not supposed to change if a read from previous conf file

}

void KMC :: init (const int L_read, const int radius_read, const double conc_read, const double T0, const bool old_conf){


    if (old_conf) {
        std :: cout << "\n Insert full path to last configuration file \n";
        std :: string filename,line;
        std :: cin >> filename;
       // std :: cout << filename;

	    std :: ifstream finput(filename);
        
	
	    if (finput.is_open()){

            std :: getline(finput,line, '\t'); //extracts until character delimiter '\t'
            L= std::stoi(line); //parse to integer
            
            island.init(L);
            adatom.init(L);

            std :: getline(finput,line, '\t'); //skip
            std :: getline(finput,line, '\t'); 
            current_T = std::stod(line);
            //std :: cout << current_T<< "\n";

            std :: getline(finput,line,'\t'); 
            concentration = std::stod(line);
            //std :: cout << concentration<< "\n";

            std :: cout << "L = " << L << " T =" << current_T << "initial concentration = " << concentration << "\n"<< std::flush ;
            std :: getline(finput,line);//skip one line
        // read island and adatom from file
            int dummy;
            
            for (int i = 0; i < L; i++){
                
                std :: getline(finput,line);
                //std :: cout << line <<std::flush;
                std::istringstream ss(line);
                for (int j = 0; j < L; j++){              
                    ss >> dummy;
                    //std:: cout << dummy<<" " <<std::flush;
                    island.matrix[i][j] =(dummy ? true : false);//a bit involved way to convert int to bool
                    //island.matrix[i][j] =0;
                 //  std:: cout << island.matrix[i][j]<<" " <<std::flush;
                }
                //std:: cout <<"\n";
            }


            std :: getline(finput,line);//empty line
            int counter =0;
            for (int i = 0; i < L; i++){
                std :: getline(finput,line);
                std::istringstream ss(line);
                for (int j = 0; j < L; j++){
                    ss>> adatom.matrix[i][j];
                   // std:: cout << adatom.matrix[i][j]<<" " <<std::flush;
                    counter += adatom.matrix[i][j];
                }
                adatom.N = counter;
                //std:: cout <<"\n";
            }
        finput.close();
        

        }
        else{
        std :: cout << "input file not found. Abort \n";
		exit(EXIT_FAILURE);
        }
    }

    else{
        L = L_read;
        radius = radius_read;
        current_T = T0;// initial temperature
        concentration = conc_read; //initial concentration
        island.init(L,radius);
        adatom.init(L,concentration);
    }

    
    

    for (int i = 0; i < n_classes; i++)
    {
        R[i].init(L);
    }
    R[0].D = det_rate(0); //nn=0
    R[1].D = det_rate(1);//I can use constant expression (finezza..)
    R[2].D = det_rate(2);
    R[3].D = det_rate(3);
    R[4].D = att_rate();
    R[5].D = 4.0 *1; //4 possible movements, diffusion constant =1 (normalization)

    island.init_neighbours(); 

// Fill classes
    for (int y = 0; y < L; y++){
        for (int x = 0; x < L; x++){
            if (island.nn1[y][x]==0&&island.matrix[y][x]){
                R[0].populate(x,y);
            }		
            if (island.nn1[y][x]==1){
                R[1].populate(x,y);
            }
            if (island.nn1[y][x]==2){
                R[2].populate(x,y);
            }
            if (island.nn1[y][x]==3){
                R[3].populate(x,y);
            }
            
            for(int k=0; k<adatom.matrix[y][x]; k++){
                // I consider possibility of multiple adatoms on top of each other 
                R[5].populate(x,y);
                if(is_attSite(x,y)) R[4].populate(x,y);
            }
        }
    }

    
}
 

double KMC ::  det_rate(const int nn) const{

    double rate;
    rate = exp((-J*nn - A)/current_T);//" attachment strenght"A=E_A-E_D with E_D: diff energy, E_A: attachment energy

    return rate;
}

double KMC ::  att_rate() const{

    double rate;
    rate = exp(-A/current_T);

    return rate;
}


double KMC :: cumulative (double* r){
    
    double R_sum =0;
    r[0] = 0;
    for (int i = 1; i <=n_classes; i++)
    {
        r[i] = r[i-1] + R[i-1].rate();
    }
    R_sum = r[n_classes];
    return R_sum;
}

int KMC ::  extract (int N) const{

    return (rand() % N);
}



bool KMC :: is_attSite(const int x,const int y) const{
	
	int contact = false ;
	int top,bottom,left,right;


	top = y+1;
	if(top==L) top = 0;

	bottom = y-1;
	if(bottom==-1) bottom = L-1;

	right = x+1;
	if (right ==L) right = 0;
	
	left = x -1;
	if(left == -1) left = L-1; 

	contact = ~island.matrix[y][x] &//not island 
    ((island.matrix[bottom][x] | island.matrix[top][x] | island.matrix[y][left] | island.matrix[y][right]));
	
	return contact;
}

double KMC :: get_concentration() const {
    return concentration;
}

int* KMC :: get_nevents() const {
    static int event_counter [n_classes];
    event_counter[0] = count_dnn0;
    event_counter[1] = count_dnn1;
    event_counter[2] = count_dnn2;
    event_counter[3] = count_dnn3;
    event_counter[4] = count_a;
    event_counter[5] = count_d;
    return event_counter;
}


void KMC :: print (int frame, int flag) const{

    auto path = "plots" + (std::to_string(proc_ID));
    auto name_a = path + "/adatom" + std::to_string(frame) + ".txt";
    auto name_b = path + "/island" + std::to_string(frame) + ".txt";

    if (flag==0){
        adatom.print(name_a, R[4].mask,current_T,concentration);
        island.print(name_b, R[1].mask,R[2].mask,R[3].mask,current_T,concentration);
    }
    else if (flag ==1){
        island.print(name_b, R[1].mask,R[2].mask,R[3].mask,current_T, concentration);
    }

}

void KMC :: print_final (const int n_frames) const{

    auto path = "plots" + (std::to_string(proc_ID));
	std :: ofstream outfile (path+"/configuration.txt");
 

	if (outfile.is_open()){

        outfile << L << "\t" << n_frames << "\t" << current_T << "\t" << concentration << "\t" << "\n";

        for ( int i = 0;i < L;i++){
            for(int j =0;j<L;j++){
                outfile << "\t"<< island.matrix[i][j]; 
            }
            outfile <<"\n";
	    }
        outfile <<"\n";
        for ( int i = 0;i < L;i++){
            for(int j =0;j<L;j++){
                outfile << "\t"<< adatom.matrix[i][j]; 
            }
            outfile <<"\n";
	    }
    }
    
    
    outfile.close();
}




 void KMC:: step(const double T, const bool debug_mode){

    
    static int step=0;
    int who;
   
	int index,i_rand;
	double R_sum,d_rand;
	int left,right,top,bottom,x,y;
    double r [n_classes+1];
    bool error;
    

    // update_rate(T)
    current_T = T;

 	R_sum = cumulative(r); 

    //std :: cout << "\n \n \n \n "<< r[0] <<  "\t" << r[1] << "\t"<< r[2] << "\t"<< r[3] << "\t"<< r[4] << "\t"<< r[5] << "\n";


 	d_rand =  ((double) rand() / (RAND_MAX)) * R_sum;

   // std :: cout << "\n \n  " <<d_rand;


/*===================================
DETACHEMENT EVENT AT A 0 NN SITE
===================================
 */


if (r[0]<d_rand && d_rand<r[1]){

    who = 0;
    count_dnn0 ++;
    index = extract(R[0].N);//simple uniform random generator
    x = R[0].where(index)[0];
    y = R[0].where(index)[1];

    R[0].destroy(index);//delete from list
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;
//no island neighbours to update : isolated element

// update attachment list
		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 

        if(R[4].exist(right,y)&& !is_attSite(right,y)){
             error=R[4].destroy_coordinates(right,y);
        }
        if(R[4].exist(left,y)&& !is_attSite(left,y)){
                error=R[4].destroy_coordinates(left,y);
        }


        if(R[4].exist(x,top)&& !is_attSite(x,top)){
               error= R[4].destroy_coordinates(x,top);
        }

        if(R[4].exist(x,bottom)&& !is_attSite(x,bottom)){
               error= R[4].destroy_coordinates(x,bottom);
        }

        R[5].populate(x,y);
}
 

/*===================================
DETACHEMENT EVENT AT A 1 NN1 and 0 NN2 SITE
===================================
 */
else if (r[1]<d_rand && d_rand<r[2]){

         who = 1;

        count_dnn1+=1;
		index = extract(R[1].N);//simple uniform random generator
		// locate event
	     
        // std:: cout << "\n counter 1    " << count_dnn1 << "\n";

 		x = R[1].where(index)[0];
		y = R[1].where(index)[1];

     // std :: cout << "\n DETACHMENT event 1 , coordinate="<< x << " , " << y<< "\n \n";

		R[1].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;


    // ****************************************

		//update R3 and  neighbours
		island.nn1[y][x] = 0;
		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 
		

		if(island.matrix[y][right]){
			 island.nn1[y][right] -= 1;
            if(island.nn1[y][right] ==1){
			 	R[1].populate(right,y); 
				error=R[2].destroy_coordinates(right,y); 
			 }
			 else if(island.nn1[y][right] ==2){
			 	R[2].populate(right,y); 
				error=R[3].destroy_coordinates(right,y); 
			 }
			 else if(island.nn1[y][right] ==3){
				 R[3].populate(right,y); 
			 }
			 else {
                 if(island.nn1[y][right] ==0){
                error= R[1].destroy_coordinates(right,y);
                 R[0].populate(right,y);
              //   exit(EXIT_FAILURE);
                 }
                //  else{
                //  exit(EXIT_FAILURE);
                //  }
                }
		}
		else if(island.matrix[y][left]) {
			island.nn1[y][left] -= 1; 
            if(island.nn1[y][left] ==1){
			 	R[1].populate(left,y); 
				error=R[2].destroy_coordinates(left,y); 
			 }
			else if(island.nn1[y][left] ==2){
			 	R[2].populate(left,y); 
				error=R[3].destroy_coordinates(left,y); 
			 }
			else if(island.nn1[y][left] ==3){
				 R[3].populate(left,y); 
			 }
			else {
                 if(island.nn1[y][left] ==0){
              //   {std :: cout << "COMPLETE DISSOLUTION" ;
                error= R[1].destroy_coordinates(left,y);
                 R[0].populate(left,y);
              //   exit(EXIT_FAILURE);
                 }
                //  else{
                //  exit(EXIT_FAILURE);
                //  }
                }
			}
		else if(island.matrix[top][x]) {
			island.nn1[top][x] -= 1; 
            if(island.nn1[top][x] ==1){
			 	R[1].populate(x,top); 
				error=R[2].destroy_coordinates(x,top); 
			 }
			else if(island.nn1[top][x] ==2){
			 	R[2].populate(x,top); 
				error=R[3].destroy_coordinates(x,top); 
			 }
			else if(island.nn1[top][x] ==3){
				 R[3].populate(x,top); 
			 }
			else {
                if(island.nn1[top][x] ==0){
               // {std :: cout << "COMPLETE DISSOLUTION" ;
               error= R[1].destroy_coordinates(x,top);
                R[0].populate(x,top);
               // exit(EXIT_FAILURE);
                }
                //  else{
                //  exit(EXIT_FAILURE);
                //  }
                }
			}
		else if(island.matrix[bottom][x]){
			island.nn1[bottom][x] -= 1;
            if(island.nn1[bottom][x] ==1){
			 	R[1].populate(x,bottom); 
				error=R[2].destroy_coordinates(x,bottom); 
			 } 
			else if(island.nn1[bottom][x] ==2){
			 	R[2].populate(x,bottom); 
				error=R[3].destroy_coordinates(x,bottom); 
			 }
			else if(island.nn1[bottom][x] ==3){
				 R[3].populate(x,bottom); 
			 }
			else {
                if(island.nn1[bottom][x] ==0){
              
               error= R[1].destroy_coordinates(x,bottom);
                R[0].populate(x,bottom);
                }
                //  else{
                //  exit(EXIT_FAILURE);
                //  }
                }
			}

        //provvisional, just to update correctly 2nn
        island.nn2[y][x] = 0;
        island.nn2[top][right] = island.get_neighbours2(right,top);
        island.nn2[top][left] = island.get_neighbours2(left,top);
        island.nn2[bottom][right] = island.get_neighbours2(right,bottom);
        island.nn2[bottom][left] = island.get_neighbours2(left,bottom);


//********************************


		//Attachment list update: dissolved island atom becomes an adatom on same site


        for(int i =0; i< adatom.matrix[y][x];i++){
            //for loop is useful if another adatom above
            R[4].populate(x,y); //I already updated # adatoms
       }

        //Take care of deleting elements from R[3]
        if(R[4].exist(right,y)&& !is_attSite(right,y)){
           error=  R[4].destroy_coordinates(right,y);
        }
        if(R[4].exist(left,y)&& !is_attSite(left,y)){
              error=  R[4].destroy_coordinates(left,y);
        }


        if(R[4].exist(x,top)&& !is_attSite(x,top)){
              error=  R[4].destroy_coordinates(x,top);
        }

        if(R[4].exist(x,bottom)&& !is_attSite(x,bottom)){
            error=    R[4].destroy_coordinates(x,bottom);
        }

// ************************************
        //Diffusion list upadate
        R[5].populate(x,y);
        


	}
/*===================================
DETACHEMENT EVENT AT A 2 NN1 SITE AND 0 NN2 SITE
===================================
 */

 	else if (r[2]<d_rand && d_rand<r[3]){

        who = 2;


	    count_dnn2+=1;

		index = extract(R[2].N);
	
		// std:: cout << "\n counter 2    " << count_dnn2 << "\n";


		x = R[2].where(index)[0];
		y = R[2].where(index)[1];

       //std :: cout << "\n DETACHMENT event 2 , coordinate="<< x << " , " << y<< "\n \n";

		R[2].destroy(index);
        adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;

     //*******************************
		//update R3, R2 and R1 and nn
	



		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 

		
        island.nn1[y][x] = 0;

		if (island.nn1[y][right] == 4){
            island.nn1[y][right] -=1;
			R[3].populate(right,y);
		}
		else if(island.nn1[y][right] ==3){
            island.nn1[y][right] -=1;
			R[2].populate(right,y);
		    error=	R[3].destroy_coordinates(right,y);
		}
		else if(island.nn1[y][right] ==2){
            island.nn1[y][right] -=1;
			R[1].populate(right,y);
		    error=	R[2].destroy_coordinates(right,y);
		}

        // strange case, maybe existing only for nearest neighbour. //note I cannot make it directly an adatom
		else if(island.nn1[y][right] ==1){
			island.nn1[y][right] -=1;
		    error=	R[1].destroy_coordinates(right,y);
          //island.nn[y][right] =0; not needed, updates of nn later
            R[0].populate(right,y);
		}
        //---------

        if (island.nn1[y][left] == 4){
            island.nn1[y][left] -=1;
			R[3].populate(left,y);
		}
		else if(island.nn1[y][left] ==3){
            island.nn1[y][left] -=1;
			R[2].populate(left,y);
		    error=R[3].destroy_coordinates(left,y);
		}
		else if(island.nn1[y][left] ==2){
            island.nn1[y][left] -=1;
			R[1].populate(left,y);
		    error=R[2].destroy_coordinates(left,y);
		}

        // strange case, maybe existing only for nearest neighbour. Furthermore clearly pose problems with det bal if I make of t an adatom.
        else if(island.nn1[y][left] ==1){
			island.nn1[y][left] -=1;
		    error=R[1].destroy_coordinates(left,y);
            R[0].populate(left,y);

		}

		if(island.nn1[top][x] ==4){
            island.nn1[top][x] -=1;
			R[3].populate(x,top);
		}
		else if(island.nn1[top][x] ==3){
            island.nn1[top][x] -=1;
			R[2].populate(x,top);
		    error=R[3].destroy_coordinates(x,top);
		}
		else if(island.nn1[top][x] ==2){
            island.nn1[top][x] -=1;
			R[1].populate(x,top);
			error=R[2].destroy_coordinates(x,top);
		}

        // strange case, maybe existing only for nearest neighbour 

		else if(island.nn1[top][x] ==1){
			island.nn1[top][x] -=1;
		    error=R[1].destroy_coordinates(x,top);
            R[0].populate(x,top);
		}
        //---------

		if(island.nn1[bottom][x] ==4){
            island.nn1[bottom][x] -=1;
			R[3].populate(x,bottom);
		}
		else if(island.nn1[bottom][x] ==3){
                island.nn1[bottom][x] -=1;
				R[2].populate(x,bottom);
		        error=R[3].destroy_coordinates(x,bottom);
		}
		else if(island.nn1[bottom][x] ==2){
            island.nn1[bottom][x] -=1;
			R[1].populate(x,bottom);
			error=R[2].destroy_coordinates(x,bottom);
		}
        // strange case, maybe existing only for nearest neighbour 

		else if(island.nn1[bottom][x] ==1){
			island.nn1[bottom][x] -=1;
			error=R[1].destroy_coordinates(x,bottom);
            R[0].populate(x,bottom);

		}
        //---------
		
		
	//**********************************************
	// Attachment list update 
    for(int i =0; i< adatom.matrix[y][x];i++){
            R[4].populate(x,y);
    }


        //Take care of deleting elements from R[3]
        if(R[4].exist(right,y)&& !is_attSite(right,y)){
             error=R[4].destroy_coordinates(right,y);
        }
        if(R[4].exist(left,y)&& !is_attSite(left,y)){
                error=R[4].destroy_coordinates(left,y);
        }


        if(R[4].exist(x,top)&& !is_attSite(x,top)){
                error=R[4].destroy_coordinates(x,top);
        }

        if(R[4].exist(x,bottom)&& !is_attSite(x,bottom)){
                error=R[4].destroy_coordinates(x,bottom);
        }

    // *******************************
     //Diffusion list update

        R[5].populate(x,y);

    //*********************************

	//neighbours update, displace within ifs case per case 

		// island.nn1[y][x] = 0;
		// island.nn1[top][x] = island.get_neighbours1(x,top); 
		// island.nn1[bottom][x] = island.get_neighbours1(x,bottom); 

		// island.nn1[y][left] =island.get_neighbours1(left,y); 
		// island.nn1[y][right] = island.get_neighbours1(right,y); 

        //PROVVISORIO
        island.nn2[y][x] = 0;
        island.nn2[top][right] = island.get_neighbours2(right,top);
        island.nn2[top][left] = island.get_neighbours2(left,top);
        island.nn2[bottom][right] = island.get_neighbours2(right,bottom);
        island.nn2[bottom][left] = island.get_neighbours2(left,bottom);

		// island.nn[bottom][left] = island.get_neighbours(left,bottom);
		// island.nn[bottom][right] =  island.get_neighbours(right,bottom);
		// island.nn[top][left] = island.get_neighbours(left,top);
		// island.nn[top][right] =  island.get_neighbours(right,top);


	
	}

/*===================================
DETACHMENT EVENT AT A 3 NN1 SITE AND 0 NN2 site
===================================
 */

    // DETACHMENT AT 3 NN SITE
	else if (r[3]<d_rand && d_rand<r[4]){

        who = 3;

	    index = extract(R[3].N);
		
        count_dnn3+=1;

		x = R[3].where(index)[0];
		y = R[3].where(index)[1];

     // std :: cout << "\n DETACHMENT event 3 , coordinate="<< x << " , " << y<< "\n \n";

		R[3].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;
        

//***************************************
		// update R2 and R3 and nn

		island.nn1[y][x] =0;

		top = y+1;
		if(top==L) top = 0;
		bottom = y-1;
		if(bottom==-1) bottom = L-1;
		right = x+1;
		if (right ==L) right = 0;
		left = x -1;	
		if(left == -1) left = L-1; 

		
        if (island.nn1[y][left] == 4){
            island.nn1[y][left]-=1;
            R[3].populate(left,y);
            }
            
        else if(island.nn1[y][left] ==3){
            island.nn1[y][left]-=1;
            R[2].populate(left,y);
            error=R[3].destroy_coordinates(left,y);
        }
        else if(island.nn1[y][left] ==2){
            island.nn1[y][left]-=1;
            R[1].populate(left,y);
            error=R[2].destroy_coordinates(left,y);
        }
        // strange case, maybe existing only for nearest neighbour 
        else if(island.nn1[y][left] ==1){
            island.nn1[y][left]-=1;		
            error=R[1].destroy_coordinates(left,y);
            R[0].populate(left,y);
        }


        //---------
        //Update other detch classes

		if (island.nn1[y][right] == 4){
            island.nn1[y][right]-=1;
			 R[3].populate(right,y);
		}
		else if(island.nn1[y][right] ==3){
            island.nn1[y][right]-=1;
			R[2].populate(right,y);
			error=R[3].destroy_coordinates(right,y);
		}
		else if(island.nn1[y][right] ==2){
            island.nn1[y][right]-=1;
			R[1].populate(right,y);
			error=R[2].destroy_coordinates(right,y);
		}

        // strange case, maybe existing only for nearest neighbour 

		else if(island.nn1[y][right] ==1){
            island.nn1[y][right]-=1;
			error=R[1].destroy_coordinates(right,y);
            R[0].populate(right,y);
		}
        //---------

		if(island.nn1[top][x] ==4){
            island.nn1[top][x]-=1;
			R[3].populate(x,top);
		}
		else if(island.nn1[top][x] ==3){
            island.nn1[top][x]-=1;
			R[2].populate(x,top);
			error=R[3].destroy_coordinates(x,top);
		}
		else if(island.nn1[top][x] ==2){
            island.nn1[top][x]-=1;
			R[1].populate(x,top);
			error=R[2].destroy_coordinates(x,top);
		}
		else if(island.nn1[top][x] ==1){
            island.nn1[top][x]-=1;
			error=R[1].destroy_coordinates(x,top);
            R[0].populate(x,top);
		}
        //---------

		if(island.nn1[bottom][x] ==4){
            island.nn1[bottom][x]-=1;
			R[3].populate(x,bottom);
		}
		else if(island.nn1[bottom][x] ==3){
            island.nn1[bottom][x]-=1;
			R[2].populate(x,bottom);
			error=R[3].destroy_coordinates(x,bottom);
		}
		else if(island.nn1[bottom][x] ==2){
            island.nn1[bottom][x]-=1;
			R[1].populate(x,bottom);
			error=R[2].destroy_coordinates(x,bottom);
		}
		else if(island.nn1[bottom][x] ==1){
            island.nn1[bottom][x]-=1;
			error=R[1].destroy_coordinates(x,bottom);
            R[0].populate(x,bottom);
		}
        //---------


		// island.nn1[top][x] = island.get_neighbours1(x,top); 
		// island.nn1[bottom][x] = island.get_neighbours1(x,bottom); 

		// island.nn1[y][left] =island.get_neighbours1(left,y); 
		// island.nn1[y][right] = island.get_neighbours1(right,y); 

        //provvisional, just to update correctly 2nn
        island.nn2[y][x] = 0;
        island.nn2[top][right] = island.get_neighbours2(right,top);
        island.nn2[top][left] = island.get_neighbours2(left,top);
        island.nn2[bottom][right] = island.get_neighbours2(right,bottom);
        island.nn2[bottom][left] = island.get_neighbours2(left,bottom);




//When second nearest..
		// island.nn[bottom][left] = island.get_neighbours(left,bottom);
		// island.nn[bottom][right] =  island.get_neighbours(right,bottom);
		// island.nn[top][left] = island.get_neighbours(left,top);
		// island.nn[top][right] =  island.get_neighbours(right,top);
// *********************************
 		//Attachemnt list update
        for(int i =0; i< adatom.matrix[y][x];i++){
            R[4].populate(x,y);
        }


        //Take care of deleting elements from R[3]
        if(R[4].exist(right,y)&& !is_attSite(right,y)){
             error=R[4].destroy_coordinates(right,y);
        }
        if(R[4].exist(left,y)&& !is_attSite(left,y)){
               error= R[4].destroy_coordinates(left,y);
        }


        if(R[4].exist(x,top)&& !is_attSite(x,top)){
                error=R[4].destroy_coordinates(x,top);
        }

        if(R[4].exist(x,bottom)&& !is_attSite(x,bottom)){
                error=R[4].destroy_coordinates(x,bottom);
        }

// ******************************
        //Diffusion list update

        R[5].populate(x,y);
	}


/*===================================
ATTACHMENT EVENT
===================================
 */

	else if (r[4]<d_rand && d_rand<r[5]){

        who = 4;

        count_a+=1;
        index = extract(R[4].N);
        
        // std:: cout << "\n counter att    " << count_a << "\n";
        

        x = R[4].where(index)[0];
		y = R[4].where(index)[1];

     //   std :: cout << "\n ATTACHMENT event, coordinate="<< x << " , " << y << "\n \n";


		//R[3].destroy(index); NO!! because all multiples adatoms are not on att site anymore!
        error=R[4].destroy_coordinates(x,y);//I have to remove in this case all multiples adatoms on the site
        error=R[5].destroy_singleCoordinate(x,y);

        adatom.matrix[y][x]-=1;
        adatom.N -=1;
        island.matrix[y][x] = 1;

// ****************************

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
        // if R[3] on site does not exist, populate if there are adatoms there.
        // If it already existed (it contains one or mutiple time) the coordinate: doesn't need to be added 
        if(! R[4].exist(right,y)&& !island.matrix[y][right]){ //or use is_attSite
            for(int i =0; i < adatom.matrix[y][right];i++){
             R[4].populate(right,y);
             }
        }
        if(!R[4].exist(left,y)&& !island.matrix[y][left]){
            for(int i =0; i < adatom.matrix[y][left];i++){
                R[4].populate(left,y);
            }
        }


        if(!R[4].exist(x,top)&& !island.matrix[top][x]){
            for(int i =0; i < adatom.matrix[top][x];i++){
                R[4].populate(x,top);
            }
        }

        if(!R[4].exist(x,bottom)&& !island.matrix[bottom][x]){
            for(int i =0; i < adatom.matrix[bottom][x];i++){
                R[4].populate(x,bottom);
            }
        }

        

    //*******************************

        // Island classes :
        int  s = 0;
        if(island.matrix[y][right]){
            s+=1;
            if(island.nn1[y][right]==3){
                error=R[3].destroy_coordinates(right,y);
            }
            else if(island.nn1[y][right]==2){
                R[3].populate(right,y);
                error=R[2].destroy_coordinates(right,y);
            }
            else if(island.nn1[y][right]==1){
                R[2].populate(right,y);
                error=R[1].destroy_coordinates(right,y);
            }
            //other exotic case..
            //island is there and no neighbours = isolated element
             else if(island.nn1[y][right]==0 && island.matrix[y][right]){
                R[1].populate(right,y);
                error=R[0].destroy_coordinates(right,y);
            }
            //island.nn1[y][x] +=1;
            island.nn1[y][right] +=1;
            

        }
        if(island.matrix[y][left]){
            s+=1;
            if(island.nn1[y][left]==3){
               error= R[3].destroy_coordinates(left,y);
            }
            else if(island.nn1[y][left]==2){
                R[3].populate(left,y);
                error=R[2].destroy_coordinates(left,y);
            }
            else if(island.nn1[y][left]==1){
                R[2].populate(left,y);
               error= R[1].destroy_coordinates(left,y);
            }
            //other exotic case..
             else if(island.nn1[y][left]==0 && island.matrix[y][left]){
                R[1].populate(left,y);
                error=R[0].destroy_coordinates(left,y);
            }
            //island.nn1[y][x] +=1;
            island.nn1[y][left] +=1;
        }
        if(island.matrix[top][x]){
            s+=1;
            if(island.nn1[top][x]==3){
                error=R[3].destroy_coordinates(x,top);
            }
            else if(island.nn1[top][x]==2){
                R[3].populate(x,top);
                error=R[2].destroy_coordinates(x,top);
            }
            else if(island.nn1[top][x]==1){
                R[2].populate(x,top);
                error=R[1].destroy_coordinates(x,top);
            }
            //other exotic case..
            else if(island.nn1[top][x]==0 && island.matrix[top][x]){
                R[1].populate(x,top);
                error=R[0].destroy_coordinates(x,top);
            }
            //island.nn1[y][x] +=1;
            island.nn1[top][x] +=1;
        }
        if(island.matrix[bottom][x]){
            s+=1;
            if(island.nn1[bottom][x]==3){
                error=R[3].destroy_coordinates(x,bottom);
            }
            else if(island.nn1[bottom][x]==2){
                R[3].populate(x,bottom);
                error=R[2].destroy_coordinates(x,bottom);
            }
            else if(island.nn1[bottom][x]==1){
                R[2].populate(x,bottom);
                error=R[1].destroy_coordinates(x,bottom);
            }
            //other exotic case..
            else if(island.nn1[bottom][x]==0 && island.matrix[bottom][x]){
                R[1].populate(x,bottom);
                error=R[0].destroy_coordinates(x,bottom);
            }
            island.nn1[bottom][x] +=1;
            //island.nn1[y][x] +=1;
        }
        
        if(s<=3){
            R[s].populate(x,y);
        }
        if(s==0){
            std :: cout << "Problem, tried to attach adatom not on att site";
            exit(EXIT_FAILURE);}

        // new ------
            island.nn1[y][x] = island.get_neighbours1(x,y);
            island.nn2[y][x] = island.get_neighbours2(x,y);
            // PROVVISORIO
            island.nn2[bottom][right] = island.get_neighbours2(right,bottom);
            island.nn2[bottom][left] = island.get_neighbours2(left,bottom);
            island.nn2[top][right] = island.get_neighbours2(right,top);
            island.nn2[top][left] = island.get_neighbours2(left,top);
            // ------
    
     }


/*===================================
DIFFUSION EVENT
===================================
 */

	else if (r[5]<d_rand && d_rand<r[6]){

        who = 5;

        count_d+=1;
        index = extract(R[5].N);

        x = R[5].where(index)[0];
		y = R[5].where(index)[1];

     //  std :: cout << "\n\n DIFFUSION event , coordinate="<< x << " , " << y << "\n \n";

        R[5].destroy(index);
        adatom.matrix[y][x] -=1;
        i_rand = rand() % 4 +1;

        if(is_attSite(x,y)){
        //remove if it was (in the previous position) on an attachment site
            error=R[4].destroy_singleCoordinate(x,y);
        }
        
        if(i_rand ==1){

            top = y+1;
            if(top==L) top = 0;
            R[5].populate(x,top);
            adatom.matrix[top][x] += 1;

            if(is_attSite(x,top)){
                R[4].populate(x,top);
            }
        }
        else if(i_rand ==2){

		    bottom = y-1;
		    if(bottom==-1) bottom = L-1;
            R[5].populate(x,bottom);
            adatom.matrix[bottom][x] += 1;
            
            if(is_attSite(x,bottom)){
                R[4].populate(x,bottom);
            }
       }
        else if(i_rand ==3){

		    right = x+1;
		    if (right ==L) right = 0;

            R[5].populate(right,y);
            adatom.matrix[y][right] += 1;

            if(is_attSite(right,y)){
                R[4].populate(right,y);
            }
       }
        else if(i_rand ==4){

		    left = x -1;	
		    if(left == -1) left = L-1; 

            R[5].populate(left,y);
            adatom.matrix[y][left] += 1;

            if(is_attSite(left,y)){

                R[4].populate(left,y);
            }
       }

    }

    concentration = static_cast<double>(adatom.N)/(L*L);//update average concentration of adatoms

    


    if(proc_ID==root_process && debug_mode){

//CHECKS*********************
    step++;
	std :: cout<< "\n  " << step <<"  KMC step \n";
    if(who==0){
        std :: cout << "\n DETACHMENT event 0, coordinate="<< x << " , " << y << "\n \n";
    }

    else if(who==1){
        std :: cout << "\n DETACHMENT event 1, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==2){
        std :: cout << "\n DETACHMENT event 2, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==3){
        std :: cout << "\n DETACHMENT event 3, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==4){
        std :: cout << "\n ATTACHMENT event, coordinate="<< x << " , " << y << "\n \n";
    }    
    else if(who==5){
        std :: cout << "\n DIFFUSION event, coordinate="<< x << " , " << y << "\n \n";
    }

    // if(error==0){
    //     std :: cout << "\n\n ERROR: Unsucesful attempt to remove"<<" x = " << x <<" y = " << y << "unexisting MULTIPLE OR SINGLE element in class list." << std :: flush;
    // }
    // else if (error==-1){
    //     std :: cout << "\n\n ERROR: Unsucesful attempt to remove"<<" x = " << x <<" y = " << y << "unexisting SINGLE element in class list."<< std :: flush;
    // }



	std :: cout<< "\n island \n";
	for ( int i = 0;i < L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< island.matrix[i][j]; 
		}
		std :: cout <<"\n";
	}

	std :: cout<< "\n First neighbours \n";

	for ( int i = 0;i <L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< island.nn1[i][j]; 
			}
		std :: cout <<"\n";
		}
std :: cout<< "\n Second neighbours \n";
    for ( int i = 0;i <L;i++){
		for(int j =0;j<L;j++){
			std::cout << "\t"<< island.nn2[i][j]; 
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
	std :: cout << "\n Elements in det class0 (isolated)"<< R[0].N <<"\telements in det class 1:  " << R[1].N << 
    "\t elements in det class 2:  " << R[2].N << "\t elements in det class 3:  " << R[3].N  << "\n";
	std :: cout << "\n Elements in att class:" << R[4].N << "\t  elements in diffusion class" << R[5].N << "\n";

    std :: cout << "\n detachment 0nn list \n";
	 for (int i = 0; i < R[0].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[0].where(i)[0]<< ","<< R[0].where(i)[1] << ")\n";
 	 }
	std :: cout << "\n detachment 1nn list \n";
	 for (int i = 0; i < R[1].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[1].where(i)[0]<< ","<< R[1].where(i)[1] << ")\n";
 	 }

	std :: cout << "\n detachment 2nn list \n";
	 for (int i = 0; i < R[2].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[2].where(i)[0]<< ","<< R[2].where(i)[1] << ")\n";
 	 }

	std :: cout << "\n detachment 3nn list \n";
	 for (int i = 0; i < R[3].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[3].where(i)[0]<< ","<< R[3].where(i)[1] << ")\n";
 	 }

	std :: cout << "\n attachment list \n";
	 for (int i = 0; i < R[4].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[4].where(i)[0]<< ","<< R[4].where(i)[1] << ")\n";
 	 }

	// std :: cout << "\n diffusion list \n";
	// for (int i = 0; i < R[5].N; i++)
 	//   {
	//  	std :: cout << i <<"\t(" << R[5].where(i)[0]<< ","<< R[5].where(i)[1] << ")\n";
 	//  }
    // std :: cout<< "\n det 1 mask \n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< R[1].mask[i][j]; 
	// 	}
	// std :: cout <<"\n";
	// }

	// std :: cout<< "\n det 2 mask \n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< R[2].mask[i][j]; 
	// 	}
	// std :: cout <<"\n";
	// }

	// std :: cout<< "\n det 3 mask \n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< R[3].mask[i][j]; 
	// 	}
	// std :: cout <<"\n";
	// }
	// std :: cout<< "\n attachment mask\n";
	// for ( int i = 0;i < L;i++){
	// 	for(int j =0;j<L;j++){
	// 		std::cout << "\t"<< R[4].mask[i][j]; 
	// 	}
	// std :: cout <<"\n";
	// }


	}

if (error){
    exit(EXIT_FAILURE);
}

}