Skip to content
Search or jump to…
Pull requests
Issues
Marketplace
Explore
 
@lucagl 
lucagl
/
KMC_pattern
Private
0
0
0
Code
Issues
4
Pull requests
Actions
Projects
1
Security
Insights
Settings
KMC_pattern/kmc.cpp
@lucagl
lucagl Stable version with second neighbours.
Latest commit dd107e4 on Sep 5, 2019
 History
 1 contributor
4061 lines (3451 sloc)  132 KB
 

#include "global.h"
#include "kmc.h"
#include "Adatom.h"
#include "Island.h"
#include "Events.h"


enum det_classes {
    // #nn1 X #nn2
    // N = nn1 + 5 * nn2
    _0x0 = 0, _1x0=1, _2x0=2, _3x0=3, _4x0=4,
    _0x1 = 5,_1x1 = 6, _2x1 = 7, _3x1 = 8, _4x1 = 9,
    _0x2 = 10, _1x2 = 11, _2x2 = 12, _3x2 = 13, _4x2 = 14,
    _0x3=15, _1x3 = 16, _2x3 = 17, _3x3 = 18, _4x3 = 19,
    _0x4 = 20, _1x4 = 21, _2x4 = 22, _3x4 = 23,  attachment=24, diffusion =25
};


KMC:: KMC(const double J_read, const double BR_read, const double A_read){

        if(proc_ID == root_process){
            std :: cout << "\n Starting KMC with " << n_classes << " classes of events \n";
        }
        J = J_read; // link strenght
        BR = BR_read;
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
    R[_0x0].D = det_rate(0,0); //nn=0
    R[_1x0].D = det_rate(1,0);//I can use constant expression (finezza..)
    R[_2x0].D = det_rate(2,0);
    R[_3x0].D = det_rate(3,0);
    R[_4x0].D = det_rate(4,0);
    R[_0x1].D = det_rate(0,1);
    R[_1x1].D = det_rate(1,1);
    R[_2x1].D = det_rate(2,1);
    R[_3x1].D = det_rate(3,1);
    R[_4x1].D = det_rate(4,1);
    R[_0x2].D = det_rate(0,2);
    R[_1x2].D = det_rate(1,2);
    R[_2x2].D = det_rate(2,2);
    R[_3x2].D = det_rate(3,2);
    R[_4x2].D = det_rate(4,2);
    R[_0x3].D = det_rate(0,3);
    R[_1x3].D = det_rate(1,3);
    R[_2x3].D = det_rate(2,3);
    R[_3x3].D = det_rate(3,3);
    R[_4x3].D = det_rate(4,3);
    R[_0x4].D = det_rate(0,4);
    R[_1x4].D = det_rate(1,4);
    R[_2x4].D = det_rate(2,4);
    R[_3x4].D = det_rate(3,4);
    
    R[attachment].D = att_rate();
    R[diffusion].D = 4.0 *1; //4 possible movements, diffusion constant =1 (normalization)

    island.init_neighbours(); 

// Fill classes
    for (int y = 0; y < L; y++){
        for (int x = 0; x < L; x++){
            if (island.nn1[y][x]==0 && island.nn2[y][x]==0 && island.matrix[y][x]){
                R[_0x0].populate(x,y);
            }		
            if (island.nn1[y][x]==1 && island.nn2[y][x]==0){
                R[_1x0].populate(x,y);
            }
            if (island.nn1[y][x]==2 && island.nn2[y][x]==0){
                R[_2x0].populate(x,y);
            }
            if (island.nn1[y][x]==3 && island.nn2[y][x]==0){
                R[_3x0].populate(x,y);
            }
            //-------------------new
            //R[4]
            if (island.nn1[y][x]==4 && island.nn2[y][x]==0){
                R[_4x0].populate(x,y);
            }
            //R[5]
            if (island.nn1[y][x]==0 && island.nn2[y][x]==1){
                R[_0x1].populate(x,y);
            }
            //R[6]
            if (island.nn1[y][x]==1 && island.nn2[y][x]==1){
                R[_1x1].populate(x,y);
            }
            //R[7]
            if (island.nn1[y][x]==2 && island.nn2[y][x]==1){
                R[_2x1].populate(x,y);
            }
            //R[8]
            if (island.nn1[y][x]==3 && island.nn2[y][x]==1){
                R[_3x1].populate(x,y);
            }
            //R[9]
            if (island.nn1[y][x]==4 && island.nn2[y][x]==1){
                R[_4x1].populate(x,y);
            }
            //R[10]
            if (island.nn1[y][x]==0 && island.nn2[y][x]==2){
                R[_0x2].populate(x,y);
            }
            //R[11]
            if (island.nn1[y][x]==1 && island.nn2[y][x]==2){
                R[_1x2].populate(x,y);
            }
            //R[12]
            if (island.nn1[y][x]==2 && island.nn2[y][x]==2){
                R[_2x2].populate(x,y);
            }
            //R[13]
            if (island.nn1[y][x]==3 && island.nn2[y][x]==2){
                R[_3x2].populate(x,y);
            }
            //R[14]
            if (island.nn1[y][x]==4 && island.nn2[y][x]==2){
                R[_4x2].populate(x,y);
            }
            //R[15]
            if (island.nn1[y][x]==0 && island.nn2[y][x]==3){
                R[_0x3].populate(x,y);
            }
            //R[16]
            if (island.nn1[y][x]==1 && island.nn2[y][x]==3){
                R[_1x3].populate(x,y);
            }
            //R[17]
            if (island.nn1[y][x]==2 && island.nn2[y][x]==3){
                R[_2x3].populate(x,y);
            }
            //R[18]
            if (island.nn1[y][x]==3 && island.nn2[y][x]==3){
                R[_3x3].populate(x,y);
            }
            //R[19]
            if (island.nn1[y][x]==4 && island.nn2[y][x]==3){
                R[_4x3].populate(x,y);
            }
            //R[20]
            if (island.nn1[y][x]==0 && island.nn2[y][x]==4){
                R[_0x4].populate(x,y);
            }
            //R[21]
            if (island.nn1[y][x]==1 && island.nn2[y][x]==4){
                R[_1x4].populate(x,y);
            }
            //R[22]
            if (island.nn1[y][x]==2 && island.nn2[y][x]==4){
                R[_2x4].populate(x,y);
            }
            //R[23]
            if (island.nn1[y][x]==3 && island.nn2[y][x]==4){
                R[_3x4].populate(x,y);
            }
            //4X4 no active site---------------------------------

            for(int k=0; k<adatom.matrix[y][x]; k++){
                // I consider possibility of multiple adatoms on top of each other 
                //R[25]
                R[diffusion].populate(x,y);
                //R24
                if(is_attSite(x,y)) R[attachment].populate(x,y);
            }
        }
    }   
}
 

double KMC ::  det_rate(const int nn1, const int nn2 ) const{

    double rate;
    rate = exp((-J*(nn1+BR*nn2) - A)/current_T);//" attachment strenght"A=E_A-E_D with E_D: diff energy, E_A: attachment energy

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
    (island.matrix[bottom][x] | island.matrix[top][x] | island.matrix[y][left] | island.matrix[y][right]
    |island.matrix[bottom][right] | island.matrix[bottom][left] | island.matrix[top][right] | island.matrix[top][left] );
	
	return contact;
}

double KMC :: get_concentration() const {
    return concentration;
}

int* KMC :: get_nevents() const {
    static int counter[n_classes];
    for (int i = 0; i < n_classes; i++)
    {
        counter[i] =event_counter[i]; 
    }
    return counter;
}


void KMC :: print (int frame, int flag) const{

    auto path = "plots" + (std::to_string(proc_ID));
    auto name_a = path + "/adatom" + std::to_string(frame) + ".txt";
    auto name_b = path + "/island" + std::to_string(frame) + ".txt";

    if (flag==0){
        adatom.print(name_a, R[attachment].mask,current_T,concentration);
        island.print(name_b, R[_1x0].mask,R[_2x0].mask,R[_3x0].mask,current_T,concentration);
    }
    else if (flag ==1){
        island.print(name_b, R[_1x0].mask,R[_2x0].mask,R[_3x0].mask,current_T, concentration);
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



bool KMC :: update_nn1DetachmentClasses(const int x, const int y){

    // IMPORTANT: to make sense neighbours must be updated later!!

    bool error =0;

    int top = y+1;
    if(top==L) top = 0;
    int bottom = y-1;
    if(bottom==-1) bottom = L-1;
    int right = x+1;
    if (right ==L) right = 0;
    int left = x -1;	
    if(left == -1) left = L-1; 

    // ----------- LEFT -----------------
    if(island.nn1[y][left] ==1){
        if(island.nn2[y][left]==0){
            R[_0x0].populate(left,y);
            error= R[_1x0].destroy_coordinates(left,y);
        }
        else if(island.nn2[y][left]==1){
        R[_0x1].populate(left,y); 
        error=R[_1x1].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==2){
            R[_0x2].populate(left,y); 
            error=R[_1x2].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==3){
            R[_0x3].populate(left,y); 
            error=R[_1x3].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==4){
            R[_0x4].populate(left,y); 
            error=R[_1x4].destroy_coordinates(left,y); 
        }
    }
    else if(island.nn1[y][left] ==2){
        if(island.nn2[y][left]==0){
            R[_1x0].populate(left,y); 
            error=R[_2x0].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==1){
            R[_1x1].populate(left,y); 
            error=R[_2x1].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==2){
            R[_1x2].populate(left,y); 
            error=R[_2x2].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==3){
            R[_1x3].populate(left,y); 
            error=R[_2x3].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==4){
            R[_1x4].populate(left,y); 
            error=R[_2x4].destroy_coordinates(left,y); 
        }
    }
    else if(island.nn1[y][left] ==3){
        if(island.nn2[y][left]==0){
        R[_2x0].populate(left,y); 
        error=R[_3x0].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==1){
            R[_2x1].populate(left,y); 
            error=R[_3x1].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==2){
            R[_2x2].populate(left,y); 
            error=R[_3x2].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==3){
            R[_2x3].populate(left,y); 
            error=R[_3x3].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==4){
            R[_2x4].populate(left,y); 
            error=R[_3x4].destroy_coordinates(left,y); 
        }
    }
    else if(island.nn1[y][left] ==4){
        if(island.nn2[y][left]==0){
        R[_3x0].populate(left,y); 
        error=R[_4x0].destroy_coordinates(left,y);
        }
        else if(island.nn2[y][left]==1){
        R[_3x1].populate(left,y); 
        error=R[_4x1].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==2){
            R[_3x2].populate(left,y); 
            error=R[_4x2].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==3){
            R[_3x3].populate(left,y); 
            error=R[_4x3].destroy_coordinates(left,y); 
        }
        else if(island.nn2[y][left]==4){
            R[_3x4].populate(left,y); 
        }
    }
    //--------- RIGHT ----------------

    if(island.nn1[y][right] ==1){
        if(island.nn2[y][right]==0){
            R[_0x0].populate(right,y);
            error= R[_1x0].destroy_coordinates(right,y);
        }
        else if(island.nn2[y][right]==1){
            R[_0x1].populate(right,y); 
            error=R[_1x1].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==2){
            R[_0x2].populate(right,y); 
            error=R[_1x2].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==3){
            R[_0x3].populate(right,y); 
            error=R[_1x3].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==4){
            R[_0x4].populate(right,y); 
            error=R[_1x4].destroy_coordinates(right,y); 
        }
    }
    else if(island.nn1[y][right] ==2){
        if(island.nn2[y][right]==0){
            R[_1x0].populate(right,y); 
            error=R[_2x0].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==1){
            R[_1x1].populate(right,y); 
            error=R[_2x1].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==2){
            R[_1x2].populate(right,y); 
            error=R[_2x2].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==3){
            R[_1x3].populate(right,y); 
            error=R[_2x3].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==4){
            R[_1x4].populate(right,y); 
            error=R[_2x4].destroy_coordinates(right,y); 
        }
    }
    else if(island.nn1[y][right] ==3){
        if(island.nn2[y][right]==0){
            R[_2x0].populate(right,y); 
            error=R[_3x0].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==1){
            R[_2x1].populate(right,y); 
            error=R[_3x1].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==2){
            R[_2x2].populate(right,y); 
            error=R[_3x2].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==3){
            R[_2x3].populate(right,y); 
            error=R[_3x3].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==4){
            R[_2x4].populate(right,y); 
            error=R[_3x4].destroy_coordinates(right,y); 
        }
    }
    else if(island.nn1[y][right] ==4){
        if(island.nn2[y][right]==0){
            R[_3x0].populate(right,y); 
            error=R[_4x0].destroy_coordinates(right,y);
        }
        else if(island.nn2[y][right]==1){
            R[_3x1].populate(right,y); 
            error=R[_4x1].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==2){
            R[_3x2].populate(right,y); 
            error=R[_4x2].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==3){
            R[_3x3].populate(right,y); 
            error=R[_4x3].destroy_coordinates(right,y); 
        }
        else if(island.nn2[y][right]==4){
            R[_3x4].populate(right,y); 
        }
    }

    //--------- TOP ----------

    if(island.nn1[top][x] ==1){
        if(island.nn2[top][x]==0){
            R[_0x0].populate(x,top);
            error= R[_1x0].destroy_coordinates(x,top);
        }
        else if(island.nn2[top][x]==1){
        R[_0x1].populate(x,top); 
        error=R[_1x1].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==2){
            R[_0x2].populate(x,top); 
            error=R[_1x2].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==3){
            R[_0x3].populate(x,top); 
            error=R[_1x3].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==4){
            R[_0x4].populate(x,top); 
            error=R[_1x4].destroy_coordinates(x,top); 
        }
    }
    else if(island.nn1[top][x] ==2){
        if(island.nn2[top][x]==0){
            R[_1x0].populate(x,top); 
            error=R[_2x0].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==1){
            R[_1x1].populate(x,top); 
            error=R[_2x1].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==2){
            R[_1x2].populate(x,top); 
            error=R[_2x2].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==3){
            R[_1x3].populate(x,top); 
            error=R[_2x3].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==4){
            R[_1x4].populate(x,top); 
            error=R[_2x4].destroy_coordinates(x,top); 
        }

    }
    else if(island.nn1[top][x] ==3){
        if(island.nn2[top][x]==0){
            R[_2x0].populate(x,top); 
            error=R[_3x0].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==1){
            R[_2x1].populate(x,top); 
            error=R[_3x1].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==2){
            R[_2x2].populate(x,top); 
            error=R[_3x2].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==3){
            R[_2x3].populate(x,top); 
            error=R[_3x3].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==4){
            R[_2x4].populate(x,top); 
            error=R[_3x4].destroy_coordinates(x,top); 
        }
    }
    else if(island.nn1[top][x] ==4){
        if(island.nn2[top][x]==0){
            R[_3x0].populate(x,top); 
            error=R[_4x0].destroy_coordinates(x,top);
        }
        else if(island.nn2[top][x]==1){
            R[_3x1].populate(x,top); 
            error=R[_4x1].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==2){
            R[_3x2].populate(x,top); 
            error=R[_4x2].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==3){
            R[_3x3].populate(x,top); 
            error=R[_4x3].destroy_coordinates(x,top); 
        }
        else if(island.nn2[top][x]==4){
            R[_3x4].populate(x,top); 
        }
    }

    //--------- BOTTOM --------------

    if(island.nn1[bottom][x] ==1){
        if(island.nn2[bottom][x]==0){
            R[_0x0].populate(x,bottom);
            error= R[_1x0].destroy_coordinates(x,bottom);
        }
        else if(island.nn2[bottom][x]==1){
        R[_0x1].populate(x,bottom); 
        error=R[_1x1].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==2){
            R[_0x2].populate(x,bottom); 
            error=R[_1x2].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==3){
            R[_0x3].populate(x,bottom); 
            error=R[_1x3].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==4){
            R[_0x4].populate(x,bottom); 
            error=R[_1x4].destroy_coordinates(x,bottom); 
        }
    }
    else if(island.nn1[bottom][x] ==2){
        if(island.nn2[bottom][x]==0){
            R[_1x0].populate(x,bottom); 
            error=R[_2x0].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==1){
            R[_1x1].populate(x,bottom); 
            error=R[_2x1].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==2){
            R[_1x2].populate(x,bottom); 
            error=R[_2x2].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==3){
            R[_1x3].populate(x,bottom); 
            error=R[_2x3].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==4){
            R[_1x4].populate(x,bottom); 
            error=R[_2x4].destroy_coordinates(x,bottom); 
        }

    }
    else if(island.nn1[bottom][x] ==3){
        if(island.nn2[bottom][x]==0){
            R[_2x0].populate(x,bottom); 
            error=R[_3x0].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==1){
            R[_2x1].populate(x,bottom); 
            error=R[_3x1].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==2){
            R[_2x2].populate(x,bottom); 
            error=R[_3x2].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==3){
            R[_2x3].populate(x,bottom); 
            error=R[_3x3].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==4){
            R[_2x4].populate(x,bottom); 
            error=R[_3x4].destroy_coordinates(x,bottom); 
        }
    }
    else if(island.nn1[bottom][x] ==4){
        if(island.nn2[bottom][x]==0){
            R[_3x0].populate(x,bottom); 
            error=R[_4x0].destroy_coordinates(x,bottom);
        }
        else if(island.nn2[bottom][x]==1){
            R[_3x1].populate(x,bottom); 
            error=R[_4x1].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==2){
            R[_3x2].populate(x,bottom); 
            error=R[_4x2].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==3){
            R[_3x3].populate(x,bottom); 
            error=R[_4x3].destroy_coordinates(x,bottom); 
        }
        else if(island.nn2[bottom][x]==4){
            R[_3x4].populate(x,bottom); 
        }
    }



    island.nn1[y][x] = 0;
    island.nn1[top][x] = island.get_neighbours1(x,top); 
    island.nn1[bottom][x] = island.get_neighbours1(x,bottom); 
    island.nn1[y][left] =island.get_neighbours1(left,y); 
    island.nn1[y][right] = island.get_neighbours1(right,y);

    return error;
}

bool KMC :: update_nn2DetachmentClasses(const int x, const int y){

    // IMPORTANT: to make sense neighbours must be updated later!!
    
    bool error =0;

    int top = y+1;
    if(top==L) top = 0;
    int bottom = y-1;
    if(bottom==-1) bottom = L-1;
    int right = x+1;
    if (right ==L) right = 0;
    int left = x -1;	
    if(left == -1) left = L-1; 

    // ----------- TOP LEFT -----------------
    if(island.nn2[top][left] ==1){
        if(island.nn1[top][left]==0){
            R[_0x0].populate(left,top);
            error= R[_0x1].destroy_coordinates(left,top);
        }
        else if(island.nn1[top][left]==1){
        R[_1x0].populate(left,top); 
        error=R[_1x1].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==2){
            R[_2x0].populate(left,top); 
            error=R[_2x1].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==3){
            R[_3x0].populate(left,top); 
            error=R[_3x1].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==4){
            R[_4x0].populate(left,top); 
            error=R[_4x1].destroy_coordinates(left,top); 
        }
    }
    else if(island.nn2[top][left] ==2){
        if(island.nn1[top][left]==0){
            R[_0x1].populate(left,top); 
            error=R[_0x2].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==1){
            R[_1x1].populate(left,top); 
            error=R[_1x2].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==2){
            R[_2x1].populate(left,top); 
            error=R[_2x2].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==3){
            R[_3x1].populate(left,top); 
            error=R[_3x2].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==4){
            R[_4x1].populate(left,top); 
            error=R[_4x2].destroy_coordinates(left,top); 
        }
    }
    else if(island.nn2[top][left] ==3){
        if(island.nn1[top][left]==0){
        R[_0x2].populate(left,top); 
        error=R[_0x3].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==1){
            R[_1x2].populate(left,top); 
            error=R[_1x3].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==2){
            R[_2x2].populate(left,top); 
            error=R[_2x3].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==3){
            R[_3x2].populate(left,top); 
            error=R[_3x3].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==4){
            R[_4x2].populate(left,top); 
            error=R[_4x3].destroy_coordinates(left,top); 
        }
    }
    else if(island.nn2[top][left] ==4){
        if(island.nn1[top][left]==0){
        R[_0x3].populate(left,top); 
        error=R[_0x4].destroy_coordinates(left,top);
        }
        else if(island.nn1[top][left]==1){
        R[_1x3].populate(left,top); 
        error=R[_1x4].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==2){
            R[_2x3].populate(left,top); 
            error=R[_2x4].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==3){
            R[_3x3].populate(left,top); 
            error=R[_3x4].destroy_coordinates(left,top); 
        }
        else if(island.nn1[top][left]==4){
            R[_4x3].populate(left,top); 
        }
    }
    //--------- TOP RIGHT ----------------

    if(island.nn2[top][right] ==1){
        if(island.nn1[top][right]==0){
            R[_0x0].populate(right,top);
            error= R[_0x1].destroy_coordinates(right,top);
        }
        else if(island.nn1[top][right]==1){
        R[_1x0].populate(right,top); 
        error=R[_1x1].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==2){
            R[_2x0].populate(right,top); 
            error=R[_2x1].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==3){
            R[_3x0].populate(right,top); 
            error=R[_3x1].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==4){
            R[_4x0].populate(right,top); 
            error=R[_4x1].destroy_coordinates(right,top); 
        }
    }
    else if(island.nn2[top][right] ==2){
        if(island.nn1[top][right]==0){
            R[_0x1].populate(right,top); 
            error=R[_0x2].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==1){
            R[_1x1].populate(right,top); 
            error=R[_1x2].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==2){
            R[_2x1].populate(right,top); 
            error=R[_2x2].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==3){
            R[_3x1].populate(right,top); 
            error=R[_3x2].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==4){
            R[_4x1].populate(right,top); 
            error=R[_4x2].destroy_coordinates(right,top); 
        }
    }
    else if(island.nn2[top][right] ==3){
        if(island.nn1[top][right]==0){
        R[_0x2].populate(right,top); 
        error=R[_0x3].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==1){
            R[_1x2].populate(right,top); 
            error=R[_1x3].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==2){
            R[_2x2].populate(right,top); 
            error=R[_2x3].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==3){
            R[_3x2].populate(right,top); 
            error=R[_3x3].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==4){
            R[_4x2].populate(right,top); 
            error=R[_4x3].destroy_coordinates(right,top); 
        }
    }
    else if(island.nn2[top][right] ==4){
        if(island.nn1[top][right]==0){
        R[_0x3].populate(right,top); 
        error=R[_0x4].destroy_coordinates(right,top);
        }
        else if(island.nn1[top][right]==1){
        R[_1x3].populate(right,top); 
        error=R[_1x4].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==2){
            R[_2x3].populate(right,top); 
            error=R[_2x4].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==3){
            R[_3x3].populate(right,top); 
            error=R[_3x4].destroy_coordinates(right,top); 
        }
        else if(island.nn1[top][right]==4){
            R[_4x3].populate(right,top); 
        }
    }

    //--------- BOTTOM LEFT ----------

    if(island.nn2[bottom][left] ==1){
        if(island.nn1[bottom][left]==0){
            R[_0x0].populate(left,bottom);
            error= R[_0x1].destroy_coordinates(left,bottom);
        }
        else if(island.nn1[bottom][left]==1){
        R[_1x0].populate(left,bottom); 
        error=R[_1x1].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==2){
            R[_2x0].populate(left,bottom); 
            error=R[_2x1].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==3){
            R[_3x0].populate(left,bottom); 
            error=R[_3x1].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==4){
            R[_4x0].populate(left,bottom); 
            error=R[_4x1].destroy_coordinates(left,bottom); 
        }
    }
    else if(island.nn2[bottom][left] ==2){
        if(island.nn1[bottom][left]==0){
            R[_0x1].populate(left,bottom); 
            error=R[_0x2].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==1){
            R[_1x1].populate(left,bottom); 
            error=R[_1x2].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==2){
            R[_2x1].populate(left,bottom); 
            error=R[_2x2].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==3){
            R[_3x1].populate(left,bottom); 
            error=R[_3x2].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==4){
            R[_4x1].populate(left,bottom); 
            error=R[_4x2].destroy_coordinates(left,bottom); 
        }
    }
    else if(island.nn2[bottom][left] ==3){
        if(island.nn1[bottom][left]==0){
        R[_0x2].populate(left,bottom); 
        error=R[_0x3].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==1){
            R[_1x2].populate(left,bottom); 
            error=R[_1x3].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==2){
            R[_2x2].populate(left,bottom); 
            error=R[_2x3].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==3){
            R[_3x2].populate(left,bottom); 
            error=R[_3x3].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==4){
            R[_4x2].populate(left,bottom); 
            error=R[_4x3].destroy_coordinates(left,bottom); 
        }
    }
    else if(island.nn2[bottom][left] ==4){
        if(island.nn1[bottom][left]==0){
        R[_0x3].populate(left,bottom); 
        error=R[_0x4].destroy_coordinates(left,bottom);
        }
        else if(island.nn1[bottom][left]==1){
        R[_1x3].populate(left,bottom); 
        error=R[_1x4].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==2){
            R[_2x3].populate(left,bottom); 
            error=R[_2x4].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==3){
            R[_3x3].populate(left,bottom); 
            error=R[_3x4].destroy_coordinates(left,bottom); 
        }
        else if(island.nn1[bottom][left]==4){
            R[_4x3].populate(left,bottom); 
        }
    }

    //--------- BOTTOM RIGHT--------------

     if(island.nn2[bottom][right] ==1){
        if(island.nn1[bottom][right]==0){
            R[_0x0].populate(right,bottom);
            error= R[_0x1].destroy_coordinates(right,bottom);
        }
        else if(island.nn1[bottom][right]==1){
        R[_1x0].populate(right,bottom); 
        error=R[_1x1].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==2){
            R[_2x0].populate(right,bottom); 
            error=R[_2x1].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==3){
            R[_3x0].populate(right,bottom); 
            error=R[_3x1].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==4){
            R[_4x0].populate(right,bottom); 
            error=R[_4x1].destroy_coordinates(right,bottom); 
        }
    }
    else if(island.nn2[bottom][right] ==2){
        if(island.nn1[bottom][right]==0){
            R[_0x1].populate(right,bottom); 
            error=R[_0x2].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==1){
            R[_1x1].populate(right,bottom); 
            error=R[_1x2].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==2){
            R[_2x1].populate(right,bottom); 
            error=R[_2x2].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==3){
            R[_3x1].populate(right,bottom); 
            error=R[_3x2].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==4){
            R[_4x1].populate(right,bottom); 
            error=R[_4x2].destroy_coordinates(right,bottom); 
        }
    }
    else if(island.nn2[bottom][right] ==3){
        if(island.nn1[bottom][right]==0){
        R[_0x2].populate(right,bottom); 
        error=R[_0x3].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==1){
            R[_1x2].populate(right,bottom); 
            error=R[_1x3].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==2){
            R[_2x2].populate(right,bottom); 
            error=R[_2x3].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==3){
            R[_3x2].populate(right,bottom); 
            error=R[_3x3].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==4){
            R[_4x2].populate(right,bottom); 
            error=R[_4x3].destroy_coordinates(right,bottom); 
        }
    }
    else if(island.nn2[bottom][right] ==4){
        if(island.nn1[bottom][right]==0){
        R[_0x3].populate(right,bottom); 
        error=R[_0x4].destroy_coordinates(right,bottom);
        }
        else if(island.nn1[bottom][right]==1){
        R[_1x3].populate(right,bottom); 
        error=R[_1x4].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==2){
            R[_2x3].populate(right,bottom); 
            error=R[_2x4].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==3){
            R[_3x3].populate(right,bottom); 
            error=R[_3x4].destroy_coordinates(right,bottom); 
        }
        else if(island.nn1[bottom][right]==4){
            R[_4x3].populate(right,bottom); 
        }
    }


    island.nn2[y][x] = 0;
    island.nn2[bottom][right] = island.get_neighbours2(right,bottom);
    island.nn2[bottom][left] = island.get_neighbours2(left,bottom);
    island.nn2[top][right] = island.get_neighbours2(right,top);
    island.nn2[top][left] = island.get_neighbours2(left,top);

return error;


}


bool KMC :: update_AttachmentClasses(const int x, const int y){

    bool error =0;

    int top = y+1;
    if(top==L) top = 0;
    int bottom = y-1;
    if(bottom==-1) bottom = L-1;
    int right = x+1;
    if (right ==L) right = 0;
    int left = x -1;	
    if(left == -1) left = L-1; 



    if(R[attachment].exist(right,y)&& !is_attSite(right,y)){
            error=R[attachment].destroy_coordinates(right,y);
    }
    if(R[attachment].exist(left,y)&& !is_attSite(left,y)){
            error=R[attachment].destroy_coordinates(left,y);
    }

    if(R[attachment].exist(x,top)&& !is_attSite(x,top)){
            error= R[attachment].destroy_coordinates(x,top);
    }

    if(R[attachment].exist(x,bottom)&& !is_attSite(x,bottom)){
            error= R[attachment].destroy_coordinates(x,bottom);
    }

    // ------------ DIAGONAL UPDATE --------------

    if(R[attachment].exist(right,top)&& !is_attSite(right,top)){
            error=R[attachment].destroy_coordinates(right,top);
    }
    if(R[attachment].exist(left,top)&& !is_attSite(left,top)){
            error=R[attachment].destroy_coordinates(left,top);
    }

    if(R[attachment].exist(right,bottom)&& !is_attSite(right,bottom)){
            error= R[attachment].destroy_coordinates(right,bottom);
    }

    if(R[attachment].exist(left,bottom)&& !is_attSite(left,bottom)){
            error= R[attachment].destroy_coordinates(left,bottom);
    }

    return error;    
}

		
	




   
 void KMC:: step(const double T, const bool debug_mode){

    
    static int step=0;
    int who=-1;
   
	int index,i_rand;
	double R_sum,d_rand;
	int x,y;
    double r [n_classes+1];
    bool error =0;
    

    // update_rate(T)
    current_T = T;

 	R_sum = cumulative(r); 

    //std :: cout << "\n \n \n \n "<< r[0] <<  "\t" << r[1] << "\t"<< r[2] << "\t"<< r[3] << "\t"<< r[4] << "\t"<< r[5] << "\n";


 	d_rand =  ((double) rand() / (RAND_MAX)) * R_sum;

   // std :: cout << "\n \n  " <<d_rand;


/*===================================
DETACHEMENT EVENT AT A NN1=0, NN2 =0 SITE
===================================
 */


if (r[_0x0]<d_rand && d_rand<r[_1x0]){

    who = _0x0;
    event_counter[_0x0] +=1;
    index = extract(R[_0x0].N);//simple uniform random generator
    x = R[_0x0].where(index)[0];
    y = R[_0x0].where(index)[1];

    R[_0x0].destroy(index);//delete from list
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

//NO NEIGHBOR DETACHMENT CLASS UPDATE SINCE ISOLATED ELEMENT

// update attachment list


	error = update_AttachmentClasses(x,y);


    //Diffusion list update

    R[diffusion].populate(x,y);

    //*********************************
}
 

/*===================================
DETACHEMENT EVENT AT NN1 =1 and  NN2 =0 SITE (SPECIAL CASE)
===================================
 */
else if (r[_1x0]<d_rand && d_rand<r[_2x0]){

         who = _1x0;

        event_counter[_1x0] +=1;
		index = extract(R[_1x0].N);//simple uniform random generator
		// locate event
	     
        // std:: cout << "\n counter 1    " << count_dnn1 << "\n";

 		x = R[_1x0].where(index)[0];
		y = R[_1x0].where(index)[1];

     // std :: cout << "\n DETACHMENT event 1 , coordinate="<< x << " , " << y<< "\n \n";

		R[_1x0].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;


    // ****************************************

		island.nn1[y][x] = 0;

		int top = y+1;
		if(top==L) top = 0;
		int bottom = y-1;
		if(bottom==-1) bottom = L-1;
		int right = x+1;
		if (right ==L) right = 0;
		int left = x -1;	
		if(left == -1) left = L-1; 
		
/*In this scenario I can be more specific and use else if for exclusive situations 
(rest of island touch detaching element only top, lef, bottom, or right)*/
		if(island.matrix[y][right]){
			island.nn1[y][right] -= 1;
            if(island.nn1[y][right] ==0){
                if(island.nn2[y][right]==0){
                    R[_0x0].populate(right,y);
                    error= R[_1x0].destroy_coordinates(right,y);
                }
                else if(island.nn2[y][right]==1){
                R[_0x1].populate(right,y); 
                error=R[_1x1].destroy_coordinates(right,y); 
                }
                else if(island.nn2[y][right]==2){
                    R[_0x2].populate(right,y); 
                    error=R[_1x2].destroy_coordinates(right,y); 
                }
                // else if(island.nn2[y][right]==3){
                //     R[_0x3].populate(right,y); 
                //     error=R[_1x3].destroy_coordinates(right,y); 
                // }
                // else if(island.nn2[y][right]==4){
                //     R[_0x4].populate(right,y); 
                //     error=R[_1x4].destroy_coordinates(right,y); 
                // }
            }
            else if(island.nn1[y][right] ==1){
                if(island.nn2[y][right]==0){
                    R[_1x0].populate(right,y); 
                    error=R[_2x0].destroy_coordinates(right,y); 
                }
                else if(island.nn2[y][right]==1){
                    R[_1x1].populate(right,y); 
                    error=R[_2x1].destroy_coordinates(right,y); 
                }
                else if(island.nn2[y][right]==2){
                    R[_1x2].populate(right,y); 
                    error=R[_2x2].destroy_coordinates(right,y); 
                }
                // else if(island.nn2[y][right]==3){
                //     R[_1x3].populate(right,y); 
                //     error=R[_2x3].destroy_coordinates(right,y); 
                // }
                // else if(island.nn2[y][right]==4){
                //     R[_1x4].populate(right,y); 
                //     error=R[_2x4].destroy_coordinates(right,y); 
                // }

			}
			else if(island.nn1[y][right] ==2){
                if(island.nn2[y][right]==0){
                    R[_2x0].populate(right,y); 
                    error=R[_3x0].destroy_coordinates(right,y); 
                }
                else if(island.nn2[y][right]==1){
                    R[_2x1].populate(right,y); 
                    error=R[_3x1].destroy_coordinates(right,y); 
                }
                else if(island.nn2[y][right]==2){
                    R[_2x2].populate(right,y); 
                    error=R[_3x2].destroy_coordinates(right,y); 
                }
                // else if(island.nn2[y][right]==3){
                //     R[_2x3].populate(right,y); 
                //     error=R[_3x3].destroy_coordinates(right,y); 
                // }
                // else if(island.nn2[y][right]==4){
                //     R[_2x4].populate(right,y); 
                //     error=R[_3x4].destroy_coordinates(right,y); 
                // }
            }
			else if(island.nn1[y][right] ==3){
                if(island.nn2[y][right]==0){
                    R[_3x0].populate(right,y); 
                    error=R[_4x0].destroy_coordinates(right,y);
                }
                else if(island.nn2[y][right]==1){
                    R[_3x1].populate(right,y); 
                    error=R[_4x1].destroy_coordinates(right,y); 
                }
                else if(island.nn2[y][right]==2){
                    R[_3x2].populate(right,y); 
                    error=R[_4x2].destroy_coordinates(right,y); 
                }
                // else if(island.nn2[y][right]==3){
                //     R[_3x3].populate(right,y); 
                //     error=R[_4x3].destroy_coordinates(right,y); 
                // }
                // else if(island.nn2[y][right]==4){
                //     R[_3x4].populate(right,y); 
                // }
			}
		}

        // ---------- LEFT ------------------
        
		else if(island.matrix[y][left]){
			island.nn1[y][left] -= 1; 
            if(island.nn1[y][left] ==0){
                if(island.nn2[y][left]==0){
                    R[_0x0].populate(left,y);
                    error= R[_1x0].destroy_coordinates(left,y);
                }
                else if(island.nn2[y][left]==1){
                R[_0x1].populate(left,y); 
                error=R[_1x1].destroy_coordinates(left,y); 
                }
                else if(island.nn2[y][left]==2){
                    R[_0x2].populate(left,y); 
                    error=R[_1x2].destroy_coordinates(left,y); 
                }
                // else if(island.nn2[y][left]==3){
                //     R[_0x3].populate(left,y); 
                //     error=R[_1x3].destroy_coordinates(left,y); 
                // }
                // else if(island.nn2[y][left]==4){
                //     R[_0x4].populate(left,y); 
                //     error=R[_1x4].destroy_coordinates(left,y); 
                // }
            }
            else if(island.nn1[y][left] ==1){
                if(island.nn2[y][left]==0){
                    R[_1x0].populate(left,y); 
                    error=R[_2x0].destroy_coordinates(left,y); 
                }
                else if(island.nn2[y][left]==1){
                    R[_1x1].populate(left,y); 
                    error=R[_2x1].destroy_coordinates(left,y); 
                }
                else if(island.nn2[y][left]==2){
                    R[_1x2].populate(left,y); 
                    error=R[_2x2].destroy_coordinates(left,y); 
                }
                // else if(island.nn2[y][left]==3){
                //     R[_1x3].populate(left,y); 
                //     error=R[_2x3].destroy_coordinates(left,y); 
                // }
                // else if(island.nn2[y][left]==4){
                //     R[_1x4].populate(left,y); 
                //     error=R[_2x4].destroy_coordinates(left,y); 
                // }
			 }
			 else if(island.nn1[y][left] ==2){
                 if(island.nn2[y][left]==0){
                    R[_2x0].populate(left,y); 
                    error=R[_3x0].destroy_coordinates(left,y); 
                }
                else if(island.nn2[y][left]==1){
                    R[_2x1].populate(left,y); 
                    error=R[_3x1].destroy_coordinates(left,y); 
                }
                else if(island.nn2[y][left]==2){
                    R[_2x2].populate(left,y); 
                    error=R[_3x2].destroy_coordinates(left,y); 
                }
                // else if(island.nn2[y][left]==3){
                //     R[_2x3].populate(left,y); 
                //     error=R[_3x3].destroy_coordinates(left,y); 
                // }
                // else if(island.nn2[y][left]==4){
                //     R[_2x4].populate(left,y); 
                //     error=R[_3x4].destroy_coordinates(left,y); 
                // }
			 }
			 else if(island.nn1[y][left] ==3){
                 if(island.nn2[y][left]==0){
				    R[_3x0].populate(left,y); 
                    error=R[_4x0].destroy_coordinates(left,y);
                 }
                 else if(island.nn2[y][left]==1){
                    R[_3x1].populate(left,y); 
                    error=R[_4x1].destroy_coordinates(left,y); 
                }
                else if(island.nn2[y][left]==2){
                    R[_3x2].populate(left,y); 
                    error=R[_4x2].destroy_coordinates(left,y); 
                }
                // else if(island.nn2[y][left]==3){
                //     R[_3x3].populate(left,y); 
                //     error=R[_4x3].destroy_coordinates(left,y); 
                // }
                // else if(island.nn2[y][left]==4){
                //     R[_3x4].populate(left,y); 
                // }
			 }
			}

            // --------- TOP ---------------------
		else if(island.matrix[top][x]) {
			island.nn1[top][x] -= 1; 
            if(island.nn1[top][x] ==0){
                if(island.nn2[top][x]==0){
                    R[_0x0].populate(x,top);
                    error= R[_1x0].destroy_coordinates(x,top);
                }
                else if(island.nn2[top][x]==1){
                R[_0x1].populate(x,top); 
                error=R[_1x1].destroy_coordinates(x,top); 
                }
                else if(island.nn2[top][x]==2){
                    R[_0x2].populate(x,top); 
                    error=R[_1x2].destroy_coordinates(x,top); 
                }
                // else if(island.nn2[top][x]==3){
                //     R[_0x3].populate(x,top); 
                //     error=R[_1x3].destroy_coordinates(x,top); 
                // }
                // else if(island.nn2[top][x]==4){
                //     R[_0x4].populate(x,top); 
                //     error=R[_1x4].destroy_coordinates(x,top); 
                // }
            }
            else if(island.nn1[top][x] ==1){
                if(island.nn2[top][x]==0){
                    R[_1x0].populate(x,top); 
                    error=R[_2x0].destroy_coordinates(x,top); 
                }
                else if(island.nn2[top][x]==1){
                    R[_1x1].populate(x,top); 
                    error=R[_2x1].destroy_coordinates(x,top); 
                }
                else if(island.nn2[top][x]==2){
                    R[_1x2].populate(x,top); 
                    error=R[_2x2].destroy_coordinates(x,top); 
                }
                // else if(island.nn2[top][x]==3){
                //     R[_1x3].populate(x,top); 
                //     error=R[_2x3].destroy_coordinates(x,top); 
                // }
                // else if(island.nn2[top][x]==4){
                //     R[_1x4].populate(x,top); 
                //     error=R[_2x4].destroy_coordinates(x,top); 
                // }
			}
			else if(island.nn1[top][x] ==2){
                if(island.nn2[top][x]==0){
                    R[_2x0].populate(x,top); 
                    error=R[_3x0].destroy_coordinates(x,top); 
                }
                else if(island.nn2[top][x]==1){
                    R[_2x1].populate(x,top); 
                    error=R[_3x1].destroy_coordinates(x,top); 
                }
                else if(island.nn2[top][x]==2){
                    R[_2x2].populate(x,top); 
                    error=R[_3x2].destroy_coordinates(x,top); 
                }
                // else if(island.nn2[top][x]==3){
                //     R[_2x3].populate(x,top); 
                //     error=R[_3x3].destroy_coordinates(x,top); 
                // }
                // else if(island.nn2[top][x]==4){
                //     R[_2x4].populate(x,top); 
                //     error=R[_3x4].destroy_coordinates(x,top); 
                // }
			}
			else if(island.nn1[top][x] ==3){
                if(island.nn2[top][x]==0){
				    R[_3x0].populate(x,top); 
                    error=R[_4x0].destroy_coordinates(x,top);
                }
                else if(island.nn2[top][x]==1){
                    R[_3x1].populate(x,top); 
                    error=R[_4x1].destroy_coordinates(x,top); 
                }
                else if(island.nn2[top][x]==2){
                    R[_3x2].populate(x,top); 
                    error=R[_4x2].destroy_coordinates(x,top); 
                }
                // else if(island.nn2[top][x]==3){
                //     R[_3x3].populate(x,top); 
                //     error=R[_4x3].destroy_coordinates(x,top); 
                // }
                // else if(island.nn2[top][x]==4){
                //     R[_3x4].populate(x,top); 
                // }
			}
		}

            //----------- BOTTOM -------------
		else if(island.matrix[bottom][x]){
			island.nn1[bottom][x] -= 1;
            if(island.nn1[bottom][x] ==0){
                if(island.nn2[bottom][x]==0){
                    R[_0x0].populate(x,bottom);
                    error= R[_1x0].destroy_coordinates(x,bottom);
                }
                else if(island.nn2[bottom][x]==1){
                R[_0x1].populate(x,bottom); 
                error=R[_1x1].destroy_coordinates(x,bottom); 
                }
                else if(island.nn2[bottom][x]==2){
                    R[_0x2].populate(x,bottom); 
                    error=R[_1x2].destroy_coordinates(x,bottom); 
                }
                // else if(island.nn2[bottom][x]==3){
                //     R[_0x3].populate(x,bottom); 
                //     error=R[_1x3].destroy_coordinates(x,bottom); 
                // }
                // else if(island.nn2[bottom][x]==4){
                //     R[_0x4].populate(x,bottom); 
                //     error=R[_1x4].destroy_coordinates(x,bottom); 
                // }
            }
            else if(island.nn1[bottom][x] ==1){
                if(island.nn2[bottom][x]==0){
                    R[_1x0].populate(x,bottom); 
                    error=R[_2x0].destroy_coordinates(x,bottom); 
                }
                else if(island.nn2[bottom][x]==1){
                    R[_1x1].populate(x,bottom); 
                    error=R[_2x1].destroy_coordinates(x,bottom); 
                }
                else if(island.nn2[bottom][x]==2){
                    R[_1x2].populate(x,bottom); 
                    error=R[_2x2].destroy_coordinates(x,bottom); 
                }
                // else if(island.nn2[bottom][x]==3){
                //     R[_1x3].populate(x,bottom); 
                //     error=R[_2x3].destroy_coordinates(x,bottom); 
                // }
                // else if(island.nn2[bottom][x]==4){
                //     R[_1x4].populate(x,bottom); 
                //     error=R[_2x4].destroy_coordinates(x,bottom); 
                // }

			}
			else if(island.nn1[bottom][x] ==2){
                if(island.nn2[bottom][x]==0){
                    R[_2x0].populate(x,bottom); 
                    error=R[_3x0].destroy_coordinates(x,bottom); 
                }
                else if(island.nn2[bottom][x]==1){
                    R[_2x1].populate(x,bottom); 
                    error=R[_3x1].destroy_coordinates(x,bottom); 
                }
                else if(island.nn2[bottom][x]==2){
                    R[_2x2].populate(x,bottom); 
                    error=R[_3x2].destroy_coordinates(x,bottom); 
                }
                // else if(island.nn2[bottom][x]==3){
                //     R[_2x3].populate(x,bottom); 
                //     error=R[_3x3].destroy_coordinates(x,bottom); 
                // }
                // else if(island.nn2[bottom][x]==4){
                //     R[_2x4].populate(x,bottom); 
                //     error=R[_3x4].destroy_coordinates(x,bottom); 
                // }
			}
			else if(island.nn1[bottom][x] ==3){
                if(island.nn2[bottom][x]==0){
				    R[_3x0].populate(x,bottom); 
                    error=R[_4x0].destroy_coordinates(x,bottom);
                }
                else if(island.nn2[bottom][x]==1){
                    R[_3x1].populate(x,bottom); 
                    error=R[_4x1].destroy_coordinates(x,bottom); 
                }
                else if(island.nn2[bottom][x]==2){
                    R[_3x2].populate(x,bottom); 
                    error=R[_4x2].destroy_coordinates(x,bottom); 
                }
                // else if(island.nn2[bottom][x]==3){
                //     R[_3x3].populate(x,bottom); 
                //     error=R[_4x3].destroy_coordinates(x,bottom); 
                // }
                // else if(island.nn2[bottom][x]==4){
                //     R[_3x4].populate(x,bottom); 
                // }
			}
		}

        //provvisional, just to update correctly 2nn
        // USELESS : NO EFFECT ON DIAGONAL ELEMENTS BECAUSE I KNOW ALREADY I DO NOT HAVE ANY DIAGONAL NEIGHBOUR
        //island.nn2[x][y] =0;
        // island.nn2[top][right] = island.get_neighbours2(right,top);
        // island.nn2[top][left] = island.get_neighbours2(left,top);
        // island.nn2[bottom][right] = island.get_neighbours2(right,bottom);
        // island.nn2[bottom][left] = island.get_neighbours2(left,bottom);


    // ********* Attachment neighbour list update


    error =update_AttachmentClasses(x,y);


    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }
   
    //Diffusion list update

    R[diffusion].populate(x,y);

    //*********************************
        


	}
/*===================================
DETACHEMENT EVENT AT A NN1 =2 SITE AND NN2 =0 SITE
===================================
 */

 	else if (r[_2x0]<d_rand && d_rand<r[_3x0]){

        who = _2x0;

	    event_counter[_2x0] +=1;

		index = extract(R[_2x0].N);
	
		// std:: cout << "\n counter 2    " << count_dnn2 << "\n";


		x = R[_2x0].where(index)[0];
		y = R[_2x0].where(index)[1];

       //std :: cout << "\n DETACHMENT event 2 , coordinate="<< x << " , " << y<< "\n \n";

		R[_2x0].destroy(index);
        adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;



        //*******************************
		//update detachment classes and nn1

		
	    error =update_nn1DetachmentClasses(x,y);

        // island.nn1[y][x] = 0;
        // island.nn1[top][x] = island.get_neighbours1(x,top); 
        // island.nn1[bottom][x] = island.get_neighbours1(x,bottom); 
        // island.nn1[y][left] =island.get_neighbours1(left,y); 
        // island.nn1[y][right] = island.get_neighbours1(right,y);
       // ********* Attachment neighbour list update


        error =update_AttachmentClasses(x,y);
        
        

        // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
        for(int i =0; i< adatom.matrix[y][x];i++){
            R[attachment].populate(x,y);
        }
   
        //Diffusion list update

        R[diffusion].populate(x,y);

    //*********************************

	}

/*===================================
DETACHMENT EVENT AT A NN1 =3 SITE AND NN2 =0 site
===================================
 */

	else if (r[_3x0]<d_rand && d_rand<r[_4x0]){

        who = _3x0;

	    index = extract(R[_3x0].N);
		
        event_counter[_3x0] +=1;

		x = R[_3x0].where(index)[0];
		y = R[_3x0].where(index)[1];

		R[_3x0].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;

		//*******************************
		//update detachment classes and nn

	    error =update_nn1DetachmentClasses(x,y);
 

        error =update_AttachmentClasses(x,y);
        
        // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
        for(int i =0; i< adatom.matrix[y][x];i++){
            R[attachment].populate(x,y);
        }
   
        //Diffusion list update

        R[diffusion].populate(x,y);
    }

 /*===================================
DETACHMENT EVENT AT A NN1 =4 SITE AND NN2 =0 site
===================================
 */
else if (r[_4x0]<d_rand && d_rand<r[_0x1]){

        who = _4x0;

	    index = extract(R[_4x0].N);
		
        event_counter[_4x0] +=1;

		x = R[_4x0].where(index)[0];
		y = R[_4x0].where(index)[1];

		R[_4x0].destroy(index);
		
		adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;

		//*******************************
		//update detachment classes and nn

	    error =update_nn1DetachmentClasses(x,y);

        error =update_AttachmentClasses(x,y);
        
        // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
        for(int i =0; i< adatom.matrix[y][x];i++){
            R[attachment].populate(x,y);
        }
   
        //Diffusion list update

        R[diffusion].populate(x,y);
}


/*===================================
DETACHMENT EVENT AT NN1= 0 SITE AND NN2=1 (SPECIAL CASE)
===================================
 */
else if (r[_0x1]<d_rand && d_rand<r[_1x1]){

    who = _0x1;

    event_counter[_0x1] +=1;
    index = extract(R[_0x1].N);//simple uniform random generator
    // locate event
        
    // std:: cout << "\n counter 1    " << count_dnn1 << "\n";

    x = R[_0x1].where(index)[0];
    y = R[_0x1].where(index)[1];

    // std :: cout << "\n DETACHMENT event 1 , coordinate="<< x << " , " << y<< "\n \n";

    R[_0x1].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;


// ****************************************

    island.nn2[y][x] = 0;

    int top = y+1;
    if(top==L) top = 0;
    int bottom = y-1;
    if(bottom==-1) bottom = L-1;
    int right = x+1;
    if (right ==L) right = 0;
    int left = x -1;	
    if(left == -1) left = L-1; 
    
/*In this scenario I can be more specific and use else if for exclusive situations 
*/

// ============ TOP RIGHT =========================
    if(island.matrix[top][right]){
        island.nn2[top][right] -= 1;
        // In the following neighbour already updated
        if(island.nn2[top][right] ==0){
            if(island.nn1[top][right]==0){
                R[_0x0].populate(right,top);
                error= R[_0x1].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right]==1){
            R[_1x0].populate(right,top); 
            error=R[_1x1].destroy_coordinates(right,top); 
            }
            else if(island.nn1[top][right]==2){
                R[_2x0].populate(right,top); 
                error=R[_2x1].destroy_coordinates(right,top); 
            }
            // THIS IS IMPOSSIBLE IF THE DETACHING ELEMENT IS 0x1
            // else if(island.nn1[top][right]==3){
            //     R[_3x0].populate(right,top); 
            //     error=R[_3x1].destroy_coordinates(right,top); 
            // }
            // else if(island.nn1[top][right]==4){
            //     R[_4x0].populate(right,top); 
            //     error=R[_4x1].destroy_coordinates(right,top); 
            // }
        }
        else if(island.nn2[top][right] ==1){
            if(island.nn1[top][right]==0){
                R[_0x1].populate(right,top); 
                error=R[_0x2].destroy_coordinates(right,top); 
            }
            else if(island.nn1[top][right]==1){
                R[_1x1].populate(right,top); 
                error=R[_1x2].destroy_coordinates(right,top); 
            }
            else if(island.nn1[top][right]==2){
                R[_2x1].populate(right,top); 
                error=R[_2x2].destroy_coordinates(right,top); 
            }
        }
        else if(island.nn2[top][right] ==2){
            if(island.nn1[top][right]==0){
                R[_0x2].populate(right,top); 
                error=R[_0x3].destroy_coordinates(right,top); 
            }
            else if(island.nn1[top][right]==1){
                R[_1x2].populate(right,top); 
                error=R[_1x3].destroy_coordinates(right,top); 
            }
            else if(island.nn1[top][right]==2){
                R[_2x2].populate(right,top); 
                error=R[_2x3].destroy_coordinates(right,top); 
            }
        }
        else if(island.nn2[top][right] ==3){
            if(island.nn1[top][right]==0){
                R[_0x3].populate(right,top); 
                error=R[_0x4].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right]==1){
                R[_1x3].populate(right,top); 
                error=R[_1x4].destroy_coordinates(right,top); 
            }
            else if(island.nn1[top][right]==2){
                R[_2x3].populate(right,top); 
                error=R[_2x4].destroy_coordinates(right,top); 
            }
        }
    }

// ============ TOP LEFT =========================

    if(island.matrix[top][left]){
        island.nn2[top][left] -= 1;
        // In the following neighbour already updated
        if(island.nn2[top][left] ==0){
            if(island.nn1[top][left]==0){
                R[_0x0].populate(left,top);
                error= R[_0x1].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left]==1){
            R[_1x0].populate(left,top); 
            error=R[_1x1].destroy_coordinates(left,top); 
            }
            else if(island.nn1[top][left]==2){
                R[_2x0].populate(left,top); 
                error=R[_2x1].destroy_coordinates(left,top); 
            }
        }
        else if(island.nn2[top][left] ==1){
            if(island.nn1[top][left]==0){
                R[_0x1].populate(left,top); 
                error=R[_0x2].destroy_coordinates(left,top); 
            }
            else if(island.nn1[top][left]==1){
                R[_1x1].populate(left,top); 
                error=R[_1x2].destroy_coordinates(left,top); 
            }
            else if(island.nn1[top][left]==2){
                R[_2x1].populate(left,top); 
                error=R[_2x2].destroy_coordinates(left,top); 
            }
        }
        else if(island.nn2[top][left] ==2){
            if(island.nn1[top][left]==0){
                R[_0x2].populate(left,top); 
                error=R[_0x3].destroy_coordinates(left,top); 
            }
            else if(island.nn1[top][left]==1){
                R[_1x2].populate(left,top); 
                error=R[_1x3].destroy_coordinates(left,top); 
            }
            else if(island.nn1[top][left]==2){
                R[_2x2].populate(left,top); 
                error=R[_2x3].destroy_coordinates(left,top); 
            }
        }
        else if(island.nn2[top][left] ==3){
            if(island.nn1[top][left]==0){
                R[_0x3].populate(left,top); 
                error=R[_0x4].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left]==1){
                R[_1x3].populate(left,top); 
                error=R[_1x4].destroy_coordinates(left,top); 
            }
            else if(island.nn1[top][left]==2){
                R[_2x3].populate(left,top); 
                error=R[_2x4].destroy_coordinates(left,top); 
            }
        }
    }

// ============ BOTTOM LEFT =========================

    if(island.matrix[bottom][left]){
        island.nn2[bottom][left] -= 1;
        // In the following neighbour already updated
        if(island.nn2[bottom][left] ==0){
            if(island.nn1[bottom][left]==0){
                R[_0x0].populate(left,bottom);
                error= R[_0x1].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left]==1){
            R[_1x0].populate(left,bottom); 
            error=R[_1x1].destroy_coordinates(left,bottom); 
            }
            else if(island.nn1[bottom][left]==2){
                R[_2x0].populate(left,bottom); 
                error=R[_2x1].destroy_coordinates(left,bottom); 
            }
        }
        else if(island.nn2[bottom][left] ==1){
            if(island.nn1[bottom][left]==0){
                R[_0x1].populate(left,bottom); 
                error=R[_0x2].destroy_coordinates(left,bottom); 
            }
            else if(island.nn1[bottom][left]==1){
                R[_1x1].populate(left,bottom); 
                error=R[_1x2].destroy_coordinates(left,bottom); 
            }
            else if(island.nn1[bottom][left]==2){
                R[_2x1].populate(left,bottom); 
                error=R[_2x2].destroy_coordinates(left,bottom); 
            }
        }
        else if(island.nn2[bottom][left] ==2){
            if(island.nn1[bottom][left]==0){
                R[_0x2].populate(left,bottom); 
                error=R[_0x3].destroy_coordinates(left,bottom); 
            }
            else if(island.nn1[bottom][left]==1){
                R[_1x2].populate(left,bottom); 
                error=R[_1x3].destroy_coordinates(left,bottom); 
            }
            else if(island.nn1[bottom][left]==2){
                R[_2x2].populate(left,bottom); 
                error=R[_2x3].destroy_coordinates(left,bottom); 
            }
        }
        else if(island.nn2[bottom][left] ==3){
            if(island.nn1[bottom][left]==0){
                R[_0x3].populate(left,bottom); 
                error=R[_0x4].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left]==1){
                R[_1x3].populate(left,bottom); 
                error=R[_1x4].destroy_coordinates(left,bottom); 
            }
            else if(island.nn1[bottom][left]==2){
                R[_2x3].populate(left,bottom); 
                error=R[_2x4].destroy_coordinates(left,bottom); 
            }
        }
    }
    

// ============ BOTTOM RIGHT =========================

    if(island.matrix[bottom][right]){
        island.nn2[bottom][right] -= 1;
        // In the following neighbour already updated
        if(island.nn2[bottom][right] ==0){
            if(island.nn1[bottom][right]==0){
                R[_0x0].populate(right,bottom);
                error= R[_0x1].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right]==1){
            R[_1x0].populate(right,bottom); 
            error=R[_1x1].destroy_coordinates(right,bottom); 
            }
            else if(island.nn1[bottom][right]==2){
                R[_2x0].populate(right,bottom); 
                error=R[_2x1].destroy_coordinates(right,bottom); 
            }
        }
        else if(island.nn2[bottom][right] ==1){
            if(island.nn1[bottom][right]==0){
                R[_0x1].populate(right,bottom); 
                error=R[_0x2].destroy_coordinates(right,bottom); 
            }
            else if(island.nn1[bottom][right]==1){
                R[_1x1].populate(right,bottom); 
                error=R[_1x2].destroy_coordinates(right,bottom); 
            }
            else if(island.nn1[bottom][right]==2){
                R[_2x1].populate(right,bottom); 
                error=R[_2x2].destroy_coordinates(right,bottom); 
            }
        }
        else if(island.nn2[bottom][right] ==2){
            if(island.nn1[bottom][right]==0){
                R[_0x2].populate(right,bottom); 
                error=R[_0x3].destroy_coordinates(right,bottom); 
            }
            else if(island.nn1[bottom][right]==1){
                R[_1x2].populate(right,bottom); 
                error=R[_1x3].destroy_coordinates(right,bottom); 
            }
            else if(island.nn1[bottom][right]==2){
                R[_2x2].populate(right,bottom); 
                error=R[_2x3].destroy_coordinates(right,bottom); 
            }
        }
        else if(island.nn2[bottom][right] ==3){
            if(island.nn1[bottom][right]==0){
                R[_0x3].populate(right,bottom); 
                error=R[_0x4].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right]==1){
                R[_1x3].populate(right,bottom); 
                error=R[_1x4].destroy_coordinates(right,bottom); 
            }
            else if(island.nn1[bottom][right]==2){
                R[_2x3].populate(right,bottom); 
                error=R[_2x4].destroy_coordinates(right,bottom); 
            }
        }
    }

// ********* Attachment neighbour list update


    error = update_AttachmentClasses(x,y);


    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }
   
    //Diffusion list update

    R[diffusion].populate(x,y);




}
/*===================================
DETACHMENT EVENT AT NN1 =1 AND NN2 = 1
===================================
 */
else if (r[_1x1]<d_rand && d_rand<r[_2x1]){

    who = _1x1;

    index = extract(R[_1x1].N);
    
    event_counter[_1x1] +=1;

    x = R[_1x1].where(index)[0];
    y = R[_1x1].where(index)[1];

    R[_1x1].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); //need to update also elements on diagonal because is a mixed case

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }
   
        //Diffusion list update

    R[diffusion].populate(x,y);
}


/*===================================
DETACHMENT EVENT AT NN1 =2 AND NN2 = 1
===================================*/
else if (r[_2x1]<d_rand && d_rand<r[_3x1]){

    who = _2x1;

    index = extract(R[_2x1].N);
    
    event_counter[_2x1] +=1;

    x = R[_2x1].where(index)[0];
    y = R[_2x1].where(index)[1];

    R[_2x1].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =3 AND NN2 = 1
===================================*/
else if (r[_3x1]<d_rand && d_rand<r[_4x1]){

    who = _3x1;

    index = extract(R[_3x1].N);
    
    event_counter[_3x1] +=1;

    x = R[_3x1].where(index)[0];
    y = R[_3x1].where(index)[1];

    R[_3x1].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =4 AND NN2 = 1
===================================*/
else if (r[_4x1]<d_rand && d_rand<r[_0x2]){

    who = _4x1;

    index = extract(R[_4x1].N);
    
    event_counter[_4x1] +=1;

    x = R[_4x1].where(index)[0];
    y = R[_4x1].where(index)[1];

    R[_4x1].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =0 AND NN2 = 2
===================================*/
else if (r[_0x2]<d_rand && d_rand<r[_1x2]){

    who = _0x2;

    index = extract(R[_0x2].N);
    
    event_counter[_0x2] +=1;

    x = R[_0x2].where(index)[0];
    y = R[_0x2].where(index)[1];

    R[_0x2].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    //update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); //need only to update this since I have 0 nn1

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =1 AND NN2 = 2
===================================*/
else if (r[_1x2]<d_rand && d_rand<r[_2x2]){

    who = _1x2;

    index = extract(R[_1x2].N);
    
    event_counter[_1x2] +=1;

    x = R[_1x2].where(index)[0];
    y = R[_1x2].where(index)[1];

    R[_1x2].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error = update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =2 AND NN2 = 2
===================================*/
else if (r[_2x2]<d_rand && d_rand<r[_3x2]){

    who = _2x2;

    index = extract(R[_2x2].N);
    
    event_counter[_2x2] +=1;

    x = R[_2x2].where(index)[0];
    y = R[_2x2].where(index)[1];

    R[_2x2].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =3 AND NN2 = 2
===================================*/
else if (r[_3x2]<d_rand && d_rand<r[_4x2]){

    who = _3x2;

    index = extract(R[_3x2].N);
    
    event_counter[_3x2] +=1;

    x = R[_3x2].where(index)[0];
    y = R[_3x2].where(index)[1];

    R[_3x2].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =4 AND NN2 = 2
===================================*/
else if (r[_4x2]<d_rand && d_rand<r[_0x3]){

    who = _4x2;

    index = extract(R[_4x2].N);
    
    event_counter[_4x2] +=1;

    x = R[_4x2].where(index)[0];
    y = R[_4x2].where(index)[1];

    R[_4x2].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =0 AND NN2 = 3
===================================*/
else if (r[_0x3]<d_rand && d_rand<r[_1x3]){

    who = _0x3;

    index = extract(R[_0x3].N);
    
    event_counter[_0x3] +=1;

    x = R[_0x3].where(index)[0];
    y = R[_0x3].where(index)[1];

    R[_0x3].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    //update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =1 AND NN2 = 3
===================================*/
else if (r[_1x3]<d_rand && d_rand<r[_2x3]){

    who = _1x3;

    index = extract(R[_1x3].N);
    
    event_counter[_1x3] +=1;

    x = R[_1x3].where(index)[0];
    y = R[_1x3].where(index)[1];

    R[_1x3].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =2 AND NN2 = 3
===================================*/
else if (r[_2x3]<d_rand && d_rand<r[_3x3]){

    who = _2x3;

    index = extract(R[_2x3].N);
    
    event_counter[_2x3] +=1;

    x = R[_2x3].where(index)[0];
    y = R[_2x3].where(index)[1];

    R[_2x3].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =3 AND NN2 = 3
===================================*/
else if (r[_3x3]<d_rand && d_rand<r[_4x3]){

    who = _3x3;

    index = extract(R[_3x3].N);
    
    event_counter[_3x3] +=1;

    x = R[_3x3].where(index)[0];
    y = R[_3x3].where(index)[1];

    R[_3x3].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =4 AND NN2 = 3
===================================*/
else if (r[_4x3]<d_rand && d_rand<r[_0x4]){

    who = _4x3;

    index = extract(R[_4x3].N);
    
    event_counter[_4x3] +=1;

    x = R[_4x3].where(index)[0];
    y = R[_4x3].where(index)[1];

    R[_4x3].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =0 AND NN2 = 4
===================================*/
else if (r[_0x4]<d_rand && d_rand<r[_1x4]){

    who = _0x4;

    index = extract(R[_0x4].N);
    
    event_counter[_0x4] +=1;

    x = R[_0x4].where(index)[0];
    y = R[_0x4].where(index)[1];

    R[_0x4].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    //update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =1 AND NN2 = 4
===================================*/
else if (r[_1x4]<d_rand && d_rand<r[_2x4]){

    who = _1x4;

    index = extract(R[_1x4].N);
    
    event_counter[_1x4] +=1;

    x = R[_1x4].where(index)[0];
    y = R[_1x4].where(index)[1];

    R[_1x4].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =2 AND NN2 = 4
===================================*/
else if (r[_2x4]<d_rand && d_rand<r[_3x4]){

    who = _2x4;

    index = extract(R[_2x4].N);
    
    event_counter[_2x4] +=1;

    x = R[_2x4].where(index)[0];
    y = R[_2x4].where(index)[1];

    R[_2x4].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
DETACHMENT EVENT AT NN1 =3 AND NN2 = 4
===================================*/
else if (r[_3x4]<d_rand && d_rand<r[attachment]){

    who = _3x4;

    index = extract(R[_3x4].N);
    
    event_counter[_3x4] +=1;

    x = R[_3x4].where(index)[0];
    y = R[_3x4].where(index)[1];

    R[_3x4].destroy(index);
    
    adatom.matrix[y][x] +=1;
    adatom.N +=1;
    island.matrix[y][x] =0;

    //*******************************
    //update detachment classes and nn

    error =update_nn1DetachmentClasses(x,y);
    error =update_nn2DetachmentClasses(x,y); 

    error =update_AttachmentClasses(x,y);
    
    // ************* ON SITE ATTACHMENT AND DIFFUSION LIST UPDATE 
    for(int i =0; i< adatom.matrix[y][x];i++){
        R[attachment].populate(x,y);
    }

    //Diffusion list update

    R[diffusion].populate(x,y);
}

/*===================================
ATTACHMENT EVENT
===================================
 */

	else if (r[attachment]<d_rand && d_rand<r[diffusion]){

        who = attachment;

        event_counter[attachment]+=1;
        index = extract(R[attachment].N);
        
        // std:: cout << "\n counter att    " << count_a << "\n";
        

        x = R[attachment].where(index)[0];
		y = R[attachment].where(index)[1];

     //   std :: cout << "\n ATTACHMENT event, coordinate="<< x << " , " << y << "\n \n";


		
        error=R[attachment].destroy_coordinates(x,y);//I have to remove in this case all multiples adatoms on the site
        error=R[diffusion].destroy_singleCoordinate(x,y);// remove a single diffusion event on same coordinate

        adatom.matrix[y][x]-=1;
        adatom.N -=1;
        island.matrix[y][x] = 1;

// ****************************

        // UPDATES OF CLASSES, HERE I NEED DIFFERENT RULES THAN PREVIOUS ONES
        //Displace in dedicated function for readibility ? 

    
        int top = y+1;
		if(top==L) top = 0;
		int bottom = y-1;
		if(bottom==-1) bottom = L-1;
		int right = x+1;
		if (right ==L) right = 0;
		int left = x -1;	
		if(left == -1) left = L-1; 

       
        // If it already existed (it contains one or mutiple time) the coordinate: doesn't need to be added 
        if(! R[attachment].exist(right,y)&& !island.matrix[y][right]){ //or use is_attSite
            for(int i =0; i < adatom.matrix[y][right];i++){
             R[attachment].populate(right,y);
             }
        }
        if(!R[attachment].exist(left,y)&& !island.matrix[y][left]){
            for(int i =0; i < adatom.matrix[y][left];i++){
                R[attachment].populate(left,y);
            }
        }


        if(!R[attachment].exist(x,top)&& !island.matrix[top][x]){
            for(int i =0; i < adatom.matrix[top][x];i++){
                R[attachment].populate(x,top);
            }
        }

        if(!R[attachment].exist(x,bottom)&& !island.matrix[bottom][x]){
            for(int i =0; i < adatom.matrix[bottom][x];i++){
                R[attachment].populate(x,bottom);
            }
        }

        // =========== UPDATE ON DIAGONAL NEIGHBOURS (NN2) ====================
        if(! R[attachment].exist(right,top)&& !island.matrix[top][right]){ //or use is_attSite
            for(int i =0; i < adatom.matrix[top][right];i++){
             R[attachment].populate(right,top);
             }
        }

        if(! R[attachment].exist(right,bottom)&& !island.matrix[bottom][right]){ 
            for(int i =0; i < adatom.matrix[bottom][right];i++){
             R[attachment].populate(right,bottom);
             }
        }

        if(! R[attachment].exist(left,top)&& !island.matrix[top][left]){ 
            for(int i =0; i < adatom.matrix[top][left];i++){
             R[attachment].populate(left,top);
             }
        }

        if(! R[attachment].exist(left,bottom)&& !island.matrix[bottom][left]){ 
            for(int i =0; i < adatom.matrix[bottom][left];i++){
             R[attachment].populate(left,bottom);
             }
        }

        

    //*******************************

        // Island classes :

        //number of neighbours not updated yet
        //changes in nn1 only on first neighbours 
        // changes in nn2 only in second neighbours (diagonal)

        //-----------    nn1  class changes ------------------

        if(island.nn1[y][right]==3){
            if(island.nn2[y][right] ==0){
                R[_4x0].populate(right,y);
                error=R[_3x0].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==1){
                R[_4x1].populate(right,y);
                error=R[_3x1].destroy_coordinates(right,y);
            }  
            else if(island.nn2[y][right] ==2){
                R[_4x2].populate(right,y);
                error=R[_3x2].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==3){
                R[_4x3].populate(right,y);
                error=R[_3x3].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==4){
                error=R[_3x4].destroy_coordinates(right,y);
            }      
        }
        else if(island.nn1[y][right]==2){
            if(island.nn2[y][right] ==0){
                R[_3x0].populate(right,y);
                error=R[_2x0].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==1){
                R[_3x1].populate(right,y);
                error=R[_2x1].destroy_coordinates(right,y);
            }  
            else if(island.nn2[y][right] ==2){
                R[_3x2].populate(right,y);
                error=R[_2x2].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==3){
                R[_3x3].populate(right,y);
                error=R[_2x3].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==4){
                R[_3x4].populate(right,y);
                error=R[_2x4].destroy_coordinates(right,y);
            }      
        }

        else if(island.nn1[y][right]==1){
            if(island.nn2[y][right] ==0){
                R[_2x0].populate(right,y);
                error=R[_1x0].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==1){
                R[_2x1].populate(right,y);
                error=R[_1x1].destroy_coordinates(right,y);
            }  
            else if(island.nn2[y][right] ==2){
                R[_2x2].populate(right,y);
                error=R[_1x2].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==3){
                R[_2x3].populate(right,y);
                error=R[_1x3].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==4){
                R[_2x4].populate(right,y);
                error=R[_1x4].destroy_coordinates(right,y);
            }  
        }

        else if(island.nn1[y][right]==0 ){
            if(island.nn2[y][right] ==0 && island.matrix[y][right]){
                R[_1x0].populate(right,y);
                error=R[_0x0].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==1){
                R[_1x1].populate(right,y);
                error=R[_0x1].destroy_coordinates(right,y);
            }  
            else if(island.nn2[y][right] ==2){
                R[_1x2].populate(right,y);
                error=R[_0x2].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==3){
                R[_1x3].populate(right,y);
                error=R[_0x3].destroy_coordinates(right,y);
            }
            else if(island.nn2[y][right] ==4){
                R[_1x4].populate(right,y);
                error=R[_0x4].destroy_coordinates(right,y);
            } 
        }
        //+++++++++++++++++ LEFT ++++++++++++++
        if(island.nn1[y][left]==3){
            if(island.nn2[y][left] ==0){
                R[_4x0].populate(left,y);
                error=R[_3x0].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==1){
                R[_4x1].populate(left,y);
                error=R[_3x1].destroy_coordinates(left,y);
            }  
            else if(island.nn2[y][left] ==2){
                R[_4x2].populate(left,y);
                error=R[_3x2].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==3){
                R[_4x3].populate(left,y);
                error=R[_3x3].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==4){
                error=R[_3x4].destroy_coordinates(left,y);
            }      
        }
        else if(island.nn1[y][left]==2){
            if(island.nn2[y][left] ==0){
                R[_3x0].populate(left,y);
                error=R[_2x0].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==1){
                R[_3x1].populate(left,y);
                error=R[_2x1].destroy_coordinates(left,y);
            }  
            else if(island.nn2[y][left] ==2){
                R[_3x2].populate(left,y);
                error=R[_2x2].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==3){
                R[_3x3].populate(left,y);
                error=R[_2x3].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==4){
                R[_3x4].populate(left,y);
                error=R[_2x4].destroy_coordinates(left,y);
            }      
        }

        else if(island.nn1[y][left]==1){
            if(island.nn2[y][left] ==0){
                R[_2x0].populate(left,y);
                error=R[_1x0].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==1){
                R[_2x1].populate(left,y);
                error=R[_1x1].destroy_coordinates(left,y);
            }  
            else if(island.nn2[y][left] ==2){
                R[_2x2].populate(left,y);
                error=R[_1x2].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==3){
                R[_2x3].populate(left,y);
                error=R[_1x3].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==4){
                R[_2x4].populate(left,y);
                error=R[_1x4].destroy_coordinates(left,y);
            }  
        }

        else if(island.nn1[y][left]==0 ){
            if(island.nn2[y][left] ==0 && island.matrix[y][left]){
                R[_1x0].populate(left,y);
                error=R[_0x0].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==1){
                R[_1x1].populate(left,y);
                error=R[_0x1].destroy_coordinates(left,y);
            }  
            else if(island.nn2[y][left] ==2){
                R[_1x2].populate(left,y);
                error=R[_0x2].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==3){
                R[_1x3].populate(left,y);
                error=R[_0x3].destroy_coordinates(left,y);
            }
            else if(island.nn2[y][left] ==4){
                R[_1x4].populate(left,y);
                error=R[_0x4].destroy_coordinates(left,y);
            } 
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++    TOP     +++++++++++++++++++++++++

        if(island.nn1[top][x]==3){
            if(island.nn2[top][x] ==0){
                R[_4x0].populate(x,top);
                error=R[_3x0].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==1){
                R[_4x1].populate(x,top);
                error=R[_3x1].destroy_coordinates(x,top);
            }  
            else if(island.nn2[top][x] ==2){
                R[_4x2].populate(x,top);
                error=R[_3x2].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==3){
                R[_4x3].populate(x,top);
                error=R[_3x3].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==4){
                error=R[_3x4].destroy_coordinates(x,top);
            }      
        }
        else if(island.nn1[top][x]==2){
            if(island.nn2[top][x] ==0){
                R[_3x0].populate(x,top);
                error=R[_2x0].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==1){
                R[_3x1].populate(x,top);
                error=R[_2x1].destroy_coordinates(x,top);
            }  
            else if(island.nn2[top][x] ==2){
                R[_3x2].populate(x,top);
                error=R[_2x2].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==3){
                R[_3x3].populate(x,top);
                error=R[_2x3].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==4){
                R[_3x4].populate(x,top);
                error=R[_2x4].destroy_coordinates(x,top);
            }      
        }

        else if(island.nn1[top][x]==1){
            if(island.nn2[top][x] ==0){
                R[_2x0].populate(x,top);
                error=R[_1x0].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==1){
                R[_2x1].populate(x,top);
                error=R[_1x1].destroy_coordinates(x,top);
            }  
            else if(island.nn2[top][x] ==2){
                R[_2x2].populate(x,top);
                error=R[_1x2].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==3){
                R[_2x3].populate(x,top);
                error=R[_1x3].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==4){
                R[_2x4].populate(x,top);
                error=R[_1x4].destroy_coordinates(x,top);
            }  
        }

        else if(island.nn1[top][x]==0 ){
            if(island.nn2[top][x] ==0 && island.matrix[top][x]){
                R[_1x0].populate(x,top);
                error=R[_0x0].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==1){
                R[_1x1].populate(x,top);
                error=R[_0x1].destroy_coordinates(x,top);
            }  
            else if(island.nn2[top][x] ==2){
                R[_1x2].populate(x,top);
                error=R[_0x2].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==3){
                R[_1x3].populate(x,top);
                error=R[_0x3].destroy_coordinates(x,top);
            }
            else if(island.nn2[top][x] ==4){
                R[_1x4].populate(x,top);
                error=R[_0x4].destroy_coordinates(x,top);
            } 
        }
//+++++++++++++++++++++++++++++++++++++++
        //++++++++++++ BOTTOM ++++++++++++++++++++++++++
       if(island.nn1[bottom][x]==3){
            if(island.nn2[bottom][x] ==0){
                R[_4x0].populate(x,bottom);
                error=R[_3x0].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==1){
                R[_4x1].populate(x,bottom);
                error=R[_3x1].destroy_coordinates(x,bottom);
            }  
            else if(island.nn2[bottom][x] ==2){
                R[_4x2].populate(x,bottom);
                error=R[_3x2].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==3){
                R[_4x3].populate(x,bottom);
                error=R[_3x3].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==4){
                error=R[_3x4].destroy_coordinates(x,bottom);
            }      
        }
        else if(island.nn1[bottom][x]==2){
            if(island.nn2[bottom][x] ==0){
                R[_3x0].populate(x,bottom);
                error=R[_2x0].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==1){
                R[_3x1].populate(x,bottom);
                error=R[_2x1].destroy_coordinates(x,bottom);
            }  
            else if(island.nn2[bottom][x] ==2){
                R[_3x2].populate(x,bottom);
                error=R[_2x2].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==3){
                R[_3x3].populate(x,bottom);
                error=R[_2x3].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==4){
                R[_3x4].populate(x,bottom);
                error=R[_2x4].destroy_coordinates(x,bottom);
            }      
        }

        else if(island.nn1[bottom][x]==1){
            if(island.nn2[bottom][x] ==0){
                R[_2x0].populate(x,bottom);
                error=R[_1x0].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==1){
                R[_2x1].populate(x,bottom);
                error=R[_1x1].destroy_coordinates(x,bottom);
            }  
            else if(island.nn2[bottom][x] ==2){
                R[_2x2].populate(x,bottom);
                error=R[_1x2].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==3){
                R[_2x3].populate(x,bottom);
                error=R[_1x3].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==4){
                R[_2x4].populate(x,bottom);
                error=R[_1x4].destroy_coordinates(x,bottom);
            }  
        }

        else if(island.nn1[bottom][x]==0 ){
            if(island.nn2[bottom][x] ==0 && island.matrix[bottom][x]){
                R[_1x0].populate(x,bottom);
                error=R[_0x0].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==1){
                R[_1x1].populate(x,bottom);
                error=R[_0x1].destroy_coordinates(x,bottom);
            }  
            else if(island.nn2[bottom][x] ==2){
                R[_1x2].populate(x,bottom);
                error=R[_0x2].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==3){
                R[_1x3].populate(x,bottom);
                error=R[_0x3].destroy_coordinates(x,bottom);
            }
            else if(island.nn2[bottom][x] ==4){
                R[_1x4].populate(x,bottom);
                error=R[_0x4].destroy_coordinates(x,bottom);
            } 
        }


        //--------------------
        //--------------- nn2 class changes --------------------


        // --------- TOP RIGHT -----------------------------

        if(island.nn2[top][right]==3){
            if(island.nn1[top][right] ==0){
                R[_0x4].populate(right,top);
                error=R[_0x3].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==1){
                R[_1x4].populate(right,top);
                error=R[_1x3].destroy_coordinates(right,top);
            }  
            else if(island.nn1[top][right] ==2){
                R[_2x4].populate(right,top);
                error=R[_2x3].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==3){
                R[_3x4].populate(right,top);
                error=R[_3x3].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==4){
                error=R[_4x3].destroy_coordinates(right,top);
            }      
        }
        else if(island.nn2[top][right]==2){
            if(island.nn1[top][right] ==0){
                R[_0x3].populate(right,top);
                error=R[_0x2].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==1){
                R[_1x3].populate(right,top);
                error=R[_1x2].destroy_coordinates(right,top);
            }  
            else if(island.nn1[top][right] ==2){
                R[_2x3].populate(right,top);
                error=R[_2x2].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==3){
                R[_3x3].populate(right,top);
                error=R[_3x2].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==4){
                R[_4x3].populate(right,top);
                error=R[_4x2].destroy_coordinates(right,top);
            }      
        }

        else if(island.nn2[top][right]==1){
            if(island.nn1[top][right] ==0){
                R[_0x2].populate(right,top);
                error=R[_0x1].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==1){
                R[_1x2].populate(right,top);
                error=R[_1x1].destroy_coordinates(right,top);
            }  
            else if(island.nn1[top][right] ==2){
                R[_2x2].populate(right,top);
                error=R[_2x1].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==3){
                R[_3x2].populate(right,top);
                error=R[_3x1].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==4){
                R[_4x2].populate(right,top);
                error=R[_4x1].destroy_coordinates(right,top);
            }      
        }

        else if(island.nn2[top][right]==0){
            if(island.nn1[top][right] ==0 && island.matrix[top][right]){
                R[_0x1].populate(right,top);
                error=R[_0x0].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==1){
                R[_1x1].populate(right,top);
                error=R[_1x0].destroy_coordinates(right,top);
            }  
            else if(island.nn1[top][right] ==2){
                R[_2x1].populate(right,top);
                error=R[_2x0].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==3){
                R[_3x1].populate(right,top);
                error=R[_3x0].destroy_coordinates(right,top);
            }
            else if(island.nn1[top][right] ==4){
                R[_4x1].populate(right,top);
                error=R[_4x0].destroy_coordinates(right,top);
            }      
        }
        //+++++++++++++++++ TOP LEFT ++++++++++++++
        if(island.nn2[top][left]==3){
            if(island.nn1[top][left] ==0){
                R[_0x4].populate(left,top);
                error=R[_0x3].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==1){
                R[_1x4].populate(left,top);
                error=R[_1x3].destroy_coordinates(left,top);
            }  
            else if(island.nn1[top][left] ==2){
                R[_2x4].populate(left,top);
                error=R[_2x3].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==3){
                R[_3x4].populate(left,top);
                error=R[_3x3].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==4){
                error=R[_4x3].destroy_coordinates(left,top);
            }      
        }
        else if(island.nn2[top][left]==2){
            if(island.nn1[top][left] ==0){
                R[_0x3].populate(left,top);
                error=R[_0x2].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==1){
                R[_1x3].populate(left,top);
                error=R[_1x2].destroy_coordinates(left,top);
            }  
            else if(island.nn1[top][left] ==2){
                R[_2x3].populate(left,top);
                error=R[_2x2].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==3){
                R[_3x3].populate(left,top);
                error=R[_3x2].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==4){
                R[_4x3].populate(left,top);
                error=R[_4x2].destroy_coordinates(left,top);
            }      
        }

        else if(island.nn2[top][left]==1){
            if(island.nn1[top][left] ==0){
                R[_0x2].populate(left,top);
                error=R[_0x1].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==1){
                R[_1x2].populate(left,top);
                error=R[_1x1].destroy_coordinates(left,top);
            }  
            else if(island.nn1[top][left] ==2){
                R[_2x2].populate(left,top);
                error=R[_2x1].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==3){
                R[_3x2].populate(left,top);
                error=R[_3x1].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==4){
                R[_4x2].populate(left,top);
                error=R[_4x1].destroy_coordinates(left,top);
            }      
        }

        else if(island.nn2[top][left]==0){
            if(island.nn1[top][left] ==0 && island.matrix[top][left]){
                R[_0x1].populate(left,top);
                error=R[_0x0].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==1){
                R[_1x1].populate(left,top);
                error=R[_1x0].destroy_coordinates(left,top);
            }  
            else if(island.nn1[top][left] ==2){
                R[_2x1].populate(left,top);
                error=R[_2x0].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==3){
                R[_3x1].populate(left,top);
                error=R[_3x0].destroy_coordinates(left,top);
            }
            else if(island.nn1[top][left] ==4){
                R[_4x1].populate(left,top);
                error=R[_4x0].destroy_coordinates(left,top);
            }      
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++    BOTTOM RIGHT     +++++++++++++++++++++++++

        if(island.nn2[bottom][right]==3){
            if(island.nn1[bottom][right] ==0){
                R[_0x4].populate(right,bottom);
                error=R[_0x3].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==1){
                R[_1x4].populate(right,bottom);
                error=R[_1x3].destroy_coordinates(right,bottom);
            }  
            else if(island.nn1[bottom][right] ==2){
                R[_2x4].populate(right,bottom);
                error=R[_2x3].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==3){
                R[_3x4].populate(right,bottom);
                error=R[_3x3].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==4){
                error=R[_4x3].destroy_coordinates(right,bottom);
            }      
        }
        else if(island.nn2[bottom][right]==2){
            if(island.nn1[bottom][right] ==0){
                R[_0x3].populate(right,bottom);
                error=R[_0x2].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==1){
                R[_1x3].populate(right,bottom);
                error=R[_1x2].destroy_coordinates(right,bottom);
            }  
            else if(island.nn1[bottom][right] ==2){
                R[_2x3].populate(right,bottom);
                error=R[_2x2].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==3){
                R[_3x3].populate(right,bottom);
                error=R[_3x2].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==4){
                R[_4x3].populate(right,bottom);
                error=R[_4x2].destroy_coordinates(right,bottom);
            }      
        }

        else if(island.nn2[bottom][right]==1){
            if(island.nn1[bottom][right] ==0){
                R[_0x2].populate(right,bottom);
                error=R[_0x1].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==1){
                R[_1x2].populate(right,bottom);
                error=R[_1x1].destroy_coordinates(right,bottom);
            }  
            else if(island.nn1[bottom][right] ==2){
                R[_2x2].populate(right,bottom);
                error=R[_2x1].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==3){
                R[_3x2].populate(right,bottom);
                error=R[_3x1].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==4){
                R[_4x2].populate(right,bottom);
                error=R[_4x1].destroy_coordinates(right,bottom);
            }      
        }

        else if(island.nn2[bottom][right]==0){
            if(island.nn1[bottom][right] ==0 && island.matrix[bottom][right]){
                R[_0x1].populate(right,bottom);
                error=R[_0x0].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==1){
                R[_1x1].populate(right,bottom);
                error=R[_1x0].destroy_coordinates(right,bottom);
            }  
            else if(island.nn1[bottom][right] ==2){
                R[_2x1].populate(right,bottom);
                error=R[_2x0].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==3){
                R[_3x1].populate(right,bottom);
                error=R[_3x0].destroy_coordinates(right,bottom);
            }
            else if(island.nn1[bottom][right] ==4){
                R[_4x1].populate(right,bottom);
                error=R[_4x0].destroy_coordinates(right,bottom);
            }      
        }
//+++++++++++++++++++++++++++++++++++++++
        //++++++++++++ BOTTOM LEFT++++++++++++++++++++++++++

       if(island.nn2[bottom][left]==3){
            if(island.nn1[bottom][left] ==0){
                R[_0x4].populate(left,bottom);
                error=R[_0x3].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==1){
                R[_1x4].populate(left,bottom);
                error=R[_1x3].destroy_coordinates(left,bottom);
            }  
            else if(island.nn1[bottom][left] ==2){
                R[_2x4].populate(left,bottom);
                error=R[_2x3].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==3){
                R[_3x4].populate(left,bottom);
                error=R[_3x3].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==4){
                error=R[_4x3].destroy_coordinates(left,bottom);
            }      
        }
        else if(island.nn2[bottom][left]==2){
            if(island.nn1[bottom][left] ==0){
                R[_0x3].populate(left,bottom);
                error=R[_0x2].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==1){
                R[_1x3].populate(left,bottom);
                error=R[_1x2].destroy_coordinates(left,bottom);
            }  
            else if(island.nn1[bottom][left] ==2){
                R[_2x3].populate(left,bottom);
                error=R[_2x2].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==3){
                R[_3x3].populate(left,bottom);
                error=R[_3x2].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==4){
                R[_4x3].populate(left,bottom);
                error=R[_4x2].destroy_coordinates(left,bottom);
            }      
        }

        else if(island.nn2[bottom][left]==1){
            if(island.nn1[bottom][left] ==0){
                R[_0x2].populate(left,bottom);
                error=R[_0x1].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==1){
                R[_1x2].populate(left,bottom);
                error=R[_1x1].destroy_coordinates(left,bottom);
            }  
            else if(island.nn1[bottom][left] ==2){
                R[_2x2].populate(left,bottom);
                error=R[_2x1].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==3){
                R[_3x2].populate(left,bottom);
                error=R[_3x1].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==4){
                R[_4x2].populate(left,bottom);
                error=R[_4x1].destroy_coordinates(left,bottom);
            }      
        }

        else if(island.nn2[bottom][left]==0){
            if(island.nn1[bottom][left] ==0 && island.matrix[bottom][left]){
                R[_0x1].populate(left,bottom);
                error=R[_0x0].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==1){
                R[_1x1].populate(left,bottom);
                error=R[_1x0].destroy_coordinates(left,bottom);
            }  
            else if(island.nn1[bottom][left] ==2){
                R[_2x1].populate(left,bottom);
                error=R[_2x0].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==3){
                R[_3x1].populate(left,bottom);
                error=R[_3x0].destroy_coordinates(left,bottom);
            }
            else if(island.nn1[bottom][left] ==4){
                R[_4x1].populate(left,bottom);
                error=R[_4x0].destroy_coordinates(left,bottom);
            }      
        }

       

//***************************************************************
//      ASSIGN CORRECT DETACHMENT CLASS TO NEWLY ATTACHED ADATOM AND UPDATE NEIGHBOURS
//----------------------------------------------------------------

        int s = island.get_neighbours1(x,y); // also update neighbours matrix at given point
        int d = island.get_neighbours2(x,y);
        if((s+d)<8){
            R[s+ 5*d].populate(x,y);
        }
        if(s==0&&d==0){
            std :: cout << s<<"\t" << d<<std::flush; 
            std :: cout << "\n Problem, tried to attach adatom in("<< x << ", " << y << ") not on att site"<<std::flush;
            exit(EXIT_FAILURE);}

       
           // INSERIRE SOPRA NEGLI IF SOPRA QUESTI UPDATES CASO PER CASO ?
            island.nn1[top][x] = island.get_neighbours1(x,top); 
            island.nn1[bottom][x] = island.get_neighbours1(x,bottom); 
            island.nn1[y][left] =island.get_neighbours1(left,y); 
            island.nn1[y][right] = island.get_neighbours1(right,y); 
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

	else if (r[diffusion]<d_rand && d_rand<r[n_classes]){

        who = diffusion;

        event_counter[diffusion]+=1;
        index = extract(R[diffusion].N);

        x = R[diffusion].where(index)[0];
		y = R[diffusion].where(index)[1];

     //  std :: cout << "\n\n DIFFUSION event , coordinate="<< x << " , " << y << "\n \n";

        R[diffusion].destroy(index);
        adatom.matrix[y][x] -=1;
        i_rand = rand() % 4 +1;

        if(is_attSite(x,y)){
        //remove if it was (in the previous position) on an attachment site
            error=R[attachment].destroy_singleCoordinate(x,y);
        }
        
        if(i_rand ==1){

            int top = y+1;
            if(top==L) top = 0;
            R[diffusion].populate(x,top);
            adatom.matrix[top][x] += 1;

            if(is_attSite(x,top)){
                R[attachment].populate(x,top);
            }
        }
        else if(i_rand ==2){

		    int bottom = y-1;
		    if(bottom==-1) bottom = L-1;
            R[diffusion].populate(x,bottom);
            adatom.matrix[bottom][x] += 1;
            
            if(is_attSite(x,bottom)){
                R[attachment].populate(x,bottom);
            }
       }
        else if(i_rand ==3){

		    int right = x+1;
		    if (right ==L) right = 0;

            R[diffusion].populate(right,y);
            adatom.matrix[y][right] += 1;

            if(is_attSite(right,y)){
                R[attachment].populate(right,y);
            }
       }
        else if(i_rand ==4){

		    int left = x -1;	
		    if(left == -1) left = L-1; 

            R[diffusion].populate(left,y);
            adatom.matrix[y][left] += 1;

            if(is_attSite(left,y)){

                R[attachment].populate(left,y);
            }
       }

    }

    concentration = static_cast<double>(adatom.N)/(L*L);//update average concentration of adatoms

  

step++;
if((proc_ID==root_process && debug_mode)||(proc_ID==root_process && error==true)){

//CHECKS*********************
    
	std :: cout<< "\n  " << step <<"  KMC step \n";
    if(who==_0x0){
        std :: cout << "\n DETACHMENT event 0x0, coordinate="<< x << " , " << y << "\n \n";
    }

    else if(who==_1x0){
        std :: cout << "\n DETACHMENT event 1x0, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_2x0){
        std :: cout << "\n DETACHMENT event 2x0, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_3x0){
        std :: cout << "\n DETACHMENT event 3x0, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_4x0){
        std :: cout << "\n DETACHMENT event 4x0, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_0x1){
        std :: cout << "\n DETACHMENT event 0x1, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_1x1){
        std :: cout << "\n DETACHMENT event 1x1, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_2x1){
        std :: cout << "\n DETACHMENT event 2x1, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_3x1){
        std :: cout << "\n DETACHMENT event 3x1, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_4x1){
        std :: cout << "\n DETACHMENT event 4x1, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_0x2){
        std :: cout << "\n DETACHMENT event 0x2, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_1x2){
        std :: cout << "\n DETACHMENT event 1x2, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_2x2){
        std :: cout << "\n DETACHMENT event 2x2, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_3x2){
        std :: cout << "\n DETACHMENT event 3x2, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_4x2){
        std :: cout << "\n DETACHMENT event 4x2, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_0x3){
        std :: cout << "\n DETACHMENT event 0x3, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_1x3){
        std :: cout << "\n DETACHMENT event 1x3, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_2x3){
        std :: cout << "\n DETACHMENT event 2x3, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_3x3){
        std :: cout << "\n DETACHMENT event 3x3, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_4x3){
        std :: cout << "\n DETACHMENT event 4x3, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_0x4){
        std :: cout << "\n DETACHMENT event 0x4, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_1x4){
        std :: cout << "\n DETACHMENT event 1x4, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_2x4){
        std :: cout << "\n DETACHMENT event 2x4, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==_3x4){
        std :: cout << "\n DETACHMENT event 3x4, coordinate="<< x << " , " << y << "\n \n";
    }
    else if(who==attachment){
        std :: cout << "\n ATTACHMENT event, coordinate="<< x << " , " << y << "\n \n";
    }    
    else if(who==diffusion){
        std :: cout << "\n DIFFUSION event, coordinate="<< x << " , " << y << "\n \n";
    }

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
	// std :: cout << "\n Elements in det class nn1 = 0, nn2=0:  "<< R[_0x0].N <<"\t elements in det class nn1 = 1, nn2=0:  " << R[_1x0].N << 
    // "\t elements in det class nn1 = 2, nn2=0:  " << R[_2x0].N << "\t elements in det class nn1 = 3, nn2=0:  " << R[_3x0].N  << "\t elements in det class nn1 = 4, nn2=0:  " << R[_4x0].N<<
    // "\n Elements in det class nn1 = 0, nn2=1:  "<< R[_0x1].N <<"\t elements in det class nn1 = 1, nn2=1:  " << R[_1x1].N << 
    // "\t elements in det class nn1 = 2, nn2=1:  " << R[_2x1].N << "\t elements in det class nn1 = 3, nn2=1:  " << R[_3x1].N  << "\t elements in det class nn1 = 4, nn2=1:  " << R[_4x1].N<<
    // "\n Elements in det class nn1 = 0, nn2=2:  "<< R[_0x2].N <<"\t elements in det class nn1 = 1, nn2=2:  " << R[_1x2].N << 
    // "\t elements in det class nn1 = 2, nn2=2:  " << R[_2x2].N << "\t elements in det class nn1 = 3, nn2=2:  " << R[_3x2].N  << "\t elements in det class nn1 = 4, nn2=2:  " << R[_4x2].N<<
    // "\n Elements in det class nn1 = 0, nn2=3:  "<< R[_0x3].N <<"\t elements in det class nn1 = 1, nn2=3:  " << R[_1x3].N << 
    // "\t elements in det class nn1 = 2, nn2=3:  "<< R[_2x3].N << "\t elements in det class nn1 = 3, nn2=3:  " << R[_3x3].N  << "\t elements in det class nn1 = 4, nn2=3:  " << R[_4x3].N<<
    // "\n Elements in det class nn1 = 0, nn2=4:  "<< R[_0x4].N <<"\t elements in det class nn1 = 1, nn2=4:  " << R[_1x4].N << 
    // "\t elements in det class nn1 = 2, nn2=4:  "<< R[_2x4].N << "\t elements in det class nn1 = 3, nn2=3:  " << R[_3x4].N  << "\n";


	std :: cout << "\n Elements in att class:" << R[attachment].N << "\t  elements in diffusion class" << R[diffusion].N <<"\n";

    std :: cout << "\n " << R[_0x0].N << " detachment 0nn1, 0nn2 list \n";
	 for (int i = 0; i < R[_0x0].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[_0x0].where(i)[0]<< ","<< R[_0x0].where(i)[1] << ")\n";
 	 }
	std :: cout << "\n " << R[_1x0].N << " detachment 1nn1, 0nn2 list \n";
	 for (int i = 0; i < R[_1x0].N; i++)
 	  {
	 	std :: cout << i <<"\t(" << R[_1x0].where(i)[0]<< ","<< R[_1x0].where(i)[1] << ")\n";
 	 }

	std :: cout << "\n " << R[_2x0].N << " detachment 2nn1, 0nn2 list \n";
	 for (int i = 0; i < R[_2x0].N; i++){
	 	std :: cout << i <<"\t(" << R[_2x0].where(i)[0]<< ","<< R[_2x0].where(i)[1] << ")\n";
 	}

	std :: cout << "\n " << R[_3x0].N << " detachment 3nn1, 0nn2 list \n";
	for (int i = 0; i < R[_3x0].N; i++){
	    std :: cout << i <<"\t(" << R[_3x0].where(i)[0]<< ","<< R[_3x0].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_4x0].N << " detachment 4nn1, 0nn2 list \n";
	for (int i = 0; i < R[_4x0].N; i++){
	 	std :: cout << i <<"\t(" << R[_4x0].where(i)[0]<< ","<< R[_4x0].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_0x1].N << "  detachment 0nn1, 1nn2 list \n";
	for (int i = 0; i < R[_0x1].N; i++){
	 	std :: cout << i <<"\t(" << R[_0x1].where(i)[0]<< ","<< R[_0x1].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_1x1].N << " detachment 1nn1, 1nn2 list \n";
	for (int i = 0; i < R[_1x1].N; i++){
	 	std :: cout << i <<"\t(" << R[_1x1].where(i)[0]<< ","<< R[_1x1].where(i)[1] << ")\n";
 	}
    std :: cout << "\n" << R[_2x1].N << "  detachment 2nn1, 1nn2 list \n";
	for (int i = 0; i < R[_2x1].N; i++){
	 	std :: cout << i <<"\t(" << R[_2x1].where(i)[0]<< ","<< R[_2x1].where(i)[1] << ")\n";
 	}
    std :: cout << "\n" << R[_3x1].N << " detachment 3nn1, 1nn2 list \n";
	for (int i = 0; i < R[_3x1].N; i++){
	 	std :: cout << i <<"\t(" << R[_3x1].where(i)[0]<< ","<< R[_3x1].where(i)[1] << ")\n";
 	}
    std :: cout << "\n" << R[_4x1].N << "  detachment 4nn1, 1nn2 list \n";
	for (int i = 0; i < R[_4x1].N; i++){
	 	std :: cout << i <<"\t(" << R[_4x1].where(i)[0]<< ","<< R[_4x1].where(i)[1] << ")\n";
 	}
    std :: cout << "\n" << R[_0x2].N << "  detachment 0nn1, 2nn2 list \n";
	for (int i = 0; i < R[_0x2].N; i++){
	 	std :: cout << i <<"\t(" << R[_0x2].where(i)[0]<< ","<< R[_0x2].where(i)[1] << ")\n";
 	}
    std :: cout << "\n" << R[_1x2].N << " detachment 1nn1, 2nn2 list \n";
	for (int i = 0; i < R[_1x2].N; i++){
	 	std :: cout << i <<"\t(" << R[_1x2].where(i)[0]<< ","<< R[_1x2].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_2x2].N << " detachment 2nn1, 2nn2 list \n";
	for (int i = 0; i < R[_2x2].N; i++){
	 	std :: cout << i <<"\t(" << R[_2x2].where(i)[0]<< ","<< R[_2x2].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_3x2].N << " detachment 3nn1, 2nn2 list \n";
	for (int i = 0; i < R[_3x2].N; i++){
	 	std :: cout << i <<"\t(" << R[_3x2].where(i)[0]<< ","<< R[_3x2].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_4x2].N << " detachment 4nn1, 2nn2 list \n";
	for (int i = 0; i < R[_4x2].N; i++){
	 	std :: cout << i <<"\t(" << R[_4x2].where(i)[0]<< ","<< R[_4x2].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_0x3].N << " detachment 0nn1, 3nn2 list \n";
	for (int i = 0; i < R[_0x3].N; i++){
	    std :: cout << i <<"\t(" << R[_0x3].where(i)[0]<< ","<< R[_0x3].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_1x3].N << " detachment 1nn1, 3nn2 list \n";
	for (int i = 0; i < R[_1x3].N; i++){
	 	std :: cout << i <<"\t(" << R[_1x3].where(i)[0]<< ","<< R[_1x3].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_2x3].N << " detachment 2nn1, 3nn2 list \n";
	for (int i = 0; i < R[_2x3].N; i++){
	 	std :: cout << i <<"\t(" << R[_2x3].where(i)[0]<< ","<< R[_2x3].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_3x3].N << " detachment 3nn1, 3nn2 list \n";
	for (int i = 0; i < R[_3x3].N; i++){
	 	std :: cout << i <<"\t(" << R[_3x3].where(i)[0]<< ","<< R[_3x3].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_4x3].N << " detachment 4nn1, 3nn2 list \n";
	for (int i = 0; i < R[_4x3].N; i++){
	 	std :: cout << i <<"\t(" << R[_4x3].where(i)[0]<< ","<< R[_4x3].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_0x4].N << " detachment 0nn1, 4nn2 list \n";
	for (int i = 0; i < R[_0x4].N; i++){
	 	std :: cout << i <<"\t(" << R[_0x4].where(i)[0]<< ","<< R[_0x4].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_1x4].N << " detachment 1nn1, 4nn2 list \n";
	for (int i = 0; i < R[_1x4].N; i++){
	 	std :: cout << i <<"\t(" << R[_1x4].where(i)[0]<< ","<< R[_1x4].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_2x4].N << " detachment 2nn1, 4nn2 list \n";
	for (int i = 0; i < R[_2x4].N; i++){
	 	std :: cout << i <<"\t(" << R[_2x4].where(i)[0]<< ","<< R[_2x4].where(i)[1] << ")\n";
 	}
    std :: cout << "\n " << R[_3x4].N << " detachment 3nn1, 4nn2 list \n";
	for (int i = 0; i < R[_3x4].N; i++){
	 	std :: cout << i <<"\t(" << R[_3x4].where(i)[0]<< ","<< R[_3x4].where(i)[1] << ")\n";
 	}


	std :: cout << "\n attachment list \n";
	 for (int i = 0; i < R[attachment].N; i++){
	 	std :: cout << i <<"\t(" << R[attachment].where(i)[0]<< ","<< R[attachment].where(i)[1] << ")\n";
 	 }
}
 //std :: cout << "\n error = " << error;
if (error){
    std:: cout << "Error = "<< error << "\n" << std :: flush ; 
    exit(EXIT_FAILURE);
}

}
