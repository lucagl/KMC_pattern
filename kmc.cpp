
#include "kmc.h"



omp_lock_t writelock1, writelock2;
static unsigned step_counter=0;

enum det_classes {
    // #nn1 X #nn2
    // N = nn1 + 5 * nn2
    _0x0 = 0, _1x0=1, _2x0=2, _3x0=3, _4x0=4,
    _0x1 = 5,_1x1 = 6, _2x1 = 7, _3x1 = 8, _4x1 = 9,
    _0x2 = 10, _1x2 = 11, _2x2 = 12, _3x2 = 13, _4x2 = 14,
    _0x3=15, _1x3 = 16, _2x3 = 17, _3x3 = 18, _4x3 = 19,
    _0x4 = 20, _1x4 = 21, _2x4 = 22, _3x4 = 23,  attachment=24, diffusion =25
};


KMC:: KMC(const double J_read, const double BR_read, const double A_read, const double E_read,const int L_read, const bool is_circle, const int radius, const double conc_read, const double T0){

    if(proc_ID == root_process){
            std :: cout << "\n Starting KMC with " << n_classes << " classes of events \n" << std :: flush;
    }
    J = J_read; // link strenght
    BR = BR_read; // Ratio between first neighbours energy and second neighbours one
    A =A_read; // attachment over diffusion parameter >0 =few attachement, <0 many attachement (diffusion dominated)
    E_shift = E_read; //Shift in detachment energy to increase equilibrium concentration
    L = L_read;
    current_T = T0;// initial temperature
    c0 = conc_read; //initial concentration
    concentration = c0;//set current concentration
    init_isCircle=is_circle;
    init_radius = radius;
    
    island = Island(L,is_circle,radius);
    adatom = Adatom(L,c0);
    // std :: cout << "\nHERE "<< std :: flush;
    // std :: cout << island.nn1[0][4]<< std :: flush;
    // std :: cout << "\nHERE "<< std :: flush;
    for (int i = 0; i < n_classes-1; i++)
    {
        R[i].init();
    }
    R[diffusion].init(true);

    
}

// void KMC :: initConv_island(double sigma){

//    // std :: cout << "\n Initialising convolution routines \n" << std :: flush;

//     island.initConv(sigma);
    
// }

// void KMC :: initConv_adatom(double sigma){

//    adatom.initConv(sigma);
    
// }


void KMC :: reset(){
    if(proc_ID == root_process){
        std :: cout << "\n RESETTING \n" << std :: flush;
    }
    island = Island(L,init_isCircle,init_radius);
    adatom = Adatom(L,c0);
    
    for (int i = 0; i < n_classes; i++)
    {
        R[i].clear();
    }

    KMC :: init();
    
    
}

//BETA NEVER TESTED------------------------------------
void KMC :: read (const std :: string filename){      
        if(proc_ID == root_process){
            std :: cout << "\n Overwriting initial configuration \n";
        }
        std :: string line;

        // std :: cin >> filename;
       // std :: cout << filename;

	    std :: ifstream finput(filename);
        
        // island.~Island();
        // adatom.~Adatom();
	    if (finput.is_open()){

            std :: getline(finput,line, '\t'); //extracts until character delimiter '\t'
            L= std::stoi(line); //parse to integer
            
            island = Island(L);
            adatom = Adatom(L);

            std :: getline(finput,line, '\t'); //skip
            std :: getline(finput,line, '\t'); 
            current_T = std::stod(line);
            //std :: cout << current_T<< "\n";

            std :: getline(finput,line,'\t'); 
            concentration = std::stod(line);
            //std :: cout << concentration<< "\n";

            std :: cout << "L = " << L << " T =" << current_T << "initial concentration = " << concentration << "\n"<< std::flush ;
            std :: getline(finput,line);//skip one line
            int dummy;
            
            for (int i = 0; i < L; i++){
                std :: getline(finput,line);
                std::istringstream ss(line);
                for (int j = 0; j < L; j++){              
                    ss >> dummy;
                    island.matrix[i][j] =(dummy ? true : false);//a bit involved way to convert int to bool
                }
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
        std :: cout << "Input file not found. Abort \n";
		exit(EXIT_FAILURE);
        }
}


void KMC :: init(){
    
// assign correct rate per class element 
    if(proc_ID == root_process){
        std :: cout << "\n Initialising classes \n"<< std :: flush;
    }
    R[_0x0].setRate(det_rate(0,0)); 
    R[_1x0].setRate(det_rate(1,0));
    R[_2x0].setRate(det_rate(2,0));
    R[_3x0].setRate(det_rate(3,0));
    R[_4x0].setRate(det_rate(4,0));
    R[_0x1].setRate(det_rate(0,1));
    R[_1x1].setRate(det_rate(1,1));
    R[_2x1].setRate(det_rate(2,1));
    R[_3x1].setRate(det_rate(3,1));
    R[_4x1].setRate(det_rate(4,1));
    R[_0x2].setRate(det_rate(0,2));
    R[_1x2].setRate(det_rate(1,2));
    R[_2x2].setRate(det_rate(2,2));
    R[_3x2].setRate(det_rate(3,2));
    R[_4x2].setRate(det_rate(4,2));
    R[_0x3].setRate(det_rate(0,3));
    R[_1x3].setRate(det_rate(1,3));
    R[_2x3].setRate(det_rate(2,3));
    R[_3x3].setRate(det_rate(3,3));
    R[_4x3].setRate(det_rate(4,3));
    R[_0x4].setRate(det_rate(0,4));
    R[_1x4].setRate(det_rate(1,4));
    R[_2x4].setRate(det_rate(2,4));
    R[_3x4].setRate(det_rate(3,4));
    
    R[attachment].setRate(att_rate());
    R[diffusion].setRate(4.0 *1); //4 possible movements, diffusion constant =1 (normalization)
    //std :: cout << "\n HERE "<< std :: flush;
    //compute all neighbour list
    island.init_neighbours(); 
   // std :: cout << "\n HERE2 "<< std :: flush;
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

    // std :: cout << "\n Initial rate per diffusing adatom  " << R[diffusion].getPeratomRate() <<" \n";
    // std :: cout << "\n total rate for diff adatom  " << R[diffusion].getRate() <<" \n";
    // std :: cout << "\n Initial rate per attachment adatom  " << R[attachment].getPeratomRate() << "  # in the class= "<< R[attachment].N <<" \n";
    // std :: cout << "\n total rate for diff adatom  " << R[attachment].getRate() <<" \n";


}
 

inline double KMC :: det_rate(const int nn1, const int nn2 ) const{

    double rate;
    rate = exp((-J*(nn1+BR*nn2) - A + E_shift)/current_T);//" attachment strenght"A=E_A-E_D with E_D: diff energy, E_A: attachment energy

    return rate;
}

inline double KMC ::  att_rate() const{

    double rate;
    rate = exp(-A/current_T);

    return rate;
}


double KMC :: cumulative (double* r){
    
    double R_sum =0;
    r[0] = 0;
    for (int i = 1; i <=n_classes; i++)
    {
        r[i] = r[i-1] + R[(i - 1)].getRate();
    }
    R_sum = r[n_classes];
    return R_sum;
}



inline bool KMC :: is_attSite(const int x,const int y) const{
	
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

int* KMC :: get_classN() const {
    static int members[n_classes];
    for (int i = 0; i < n_classes; i++)
    {
        members[i] =R[i].N; 
    }
    return members;
}


void KMC :: saveTxt (const std:: string path, int frame,  bool isConv, bool flag) const{


    if (flag==0){
        if(!isConv){
                auto name_a = path + "/adatom" + std::to_string(frame) + ".txt";
                auto name_b = path + "/island" + std::to_string(frame) + ".txt";
                adatom.saveTxt(name_a,current_T,concentration);
                island.saveTxt(name_b,current_T,concentration);
        }
        else
        {
            auto name_a = path + "/adatom_conv" + std::to_string(frame) + ".txt";
            auto name_b = path + "/island_conv" + std::to_string(frame) + ".txt";
            adatom.saveTxt_conv(name_a,current_T,concentration);
            island.saveTxt_conv(name_b,current_T,concentration);
        }
        
    }
    else if (flag ==1){
        if(!isConv){
           
            auto name_b = path + "/island" + std::to_string(frame) + ".txt";
            island.saveTxt(name_b,current_T,concentration);
        }
        else
        {
            auto name_b = path + "/island_conv" + std::to_string(frame) + ".txt";
            island.saveTxt_conv(name_b,current_T,concentration);
        }
    }

}   


void KMC :: print_final (const int n_frames, bool isConv=0){


    auto path = "plots" + (std::to_string(proc_ID));
	std :: ofstream outfile1 (path+"/configuration.txt");

    //KMC :: initConv(sigma,1);//1 initialises both island and adatom convolution if not already initialised


	if (outfile1.is_open()){
        outfile1 << L << "\t" << n_frames << "\t" << current_T << "\t" << concentration << "\t" << "\n";
        for ( int i = 0;i < L;i++){
            for(int j =0;j<L;j++){
                outfile1 << "\t"<< island.matrix[i][j]; 
            }
            outfile1 <<"\n";
	    }
        outfile1 <<"\n";
        for ( int i = 0;i < L;i++){
            for(int j =0;j<L;j++){
                outfile1 << "\t"<< adatom.matrix[i][j]; 
            }
            outfile1 <<"\n";
	    }
    }
    outfile1.close();

    if (isConv){

        double** outIsland;
        double** outAdatom;
        std :: ofstream outfile2 (path+"/final_islandConv.txt");
        std :: ofstream outfile3 (path+"/final_adatomConv.txt");

        if(!adatom.isConv()){
            std :: cout <<"Adatom convolution not initialized, using defaut sigma"<< std :: endl;
            double sigma = L/50;
            adatom.initConv(sigma);
            }
        if(!island.isConv()){
            std :: cout <<"Island convolution not initialized, using defaut sigma"<< std :: endl;
            double sigma = L/100;
            island.initConv(sigma);
            }
    
        outIsland = island.gaussianConv();
        outAdatom = adatom.gaussianConv();

        if (outfile2.is_open()&&outfile3.is_open()){
                for ( int i = 0;i < L;i++){
                    for(int j =0;j<L;j++){
                        outfile2 <<  outIsland[i][j] << "\n"; 
                        outfile3 <<  outAdatom[i][j] << "\n"; 
                    }
                }
            outfile2.close();
            outfile3.close();
        }

    }
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

		
	




   
 double KMC:: step(const double T, const bool debug_mode){

    unsigned who = 100;//dummy value
   
	int index;
    int i_rand;
	double R_sum,d_rand;
	int x,y;
    double r [n_classes+1];
    bool error = 0;

    current_T = T;

 	R_sum = cumulative(r); 

 	d_rand= ((double) rand() / (RAND_MAX)) * R_sum;

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
			}
		}

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


		x = R[_2x0].where(index)[0];
		y = R[_2x0].where(index)[1];


		R[_2x0].destroy(index);
        adatom.matrix[y][x] +=1;
        adatom.N +=1;
		island.matrix[y][x] =0;

        //*******************************
		//update detachment classes and nn1

		
	    error =update_nn1DetachmentClasses(x,y);


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

    x = R[_0x1].where(index)[0];
    y = R[_0x1].where(index)[1];


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

        
        
        omp_init_lock(&writelock1);
        omp_init_lock(&writelock2);


        who = diffusion;

        event_counter[diffusion]+=1;
   

/* OBSERVATION
Explicit locking and atomic clause have same efficiency..
Each random number sequence independent for each trheads. Indeed rand() not thread safe since it relies on a shared static seed (advancing at each call).
Simultaneaous call might lead to repetition od same random number..  
*/

/*PROBLEMS
The problem is that the change of content of an index with list can problems since list are not stored in continuous memory locations.
Therefore I cannot have random access like vector or arrays.. I have to iterate on a pointer which is shared!
*/
std :: vector <std :: tuple<int, int>> Diff_adatoms{};
for (int i = 0; i < R[diffusion].N; i++){
            Diff_adatoms.push_back(std :: make_tuple (R[diffusion].where(i)[0],R[diffusion].where(i)[1]));
        }
    
    #pragma omp parallel private (i_rand,x,y) num_threads(n_threads)
    {   
         //OBS: 
        /* Localseed must be set externally or at every call new memory location assigned and rand_r will re-initialise the sequence.
        Observe that I cannot have same seed or every thread has exactly same sequence.
        The use of time() for seeding might be dangerous since the variation might be very small if loop is fast. 
        However, numbers within the same sequence are random but not same index among different sequences.
        See: https://blogs.unity3d.com/2015/01/07/a-primer-on-repeatable-random-numbers 
        This problem persist, but I think is less severe in the new implementation. 
        Indeed, before I had this correlation buildig up at every call of a diffusion event on random sequneces long as the subloop per thread. 
        Now "only" among threads whose indexing is shuffling at every for loop*/

        int id = omp_get_thread_num();


        #pragma omp for schedule(dynamic) nowait//further could prevent correlations in random sequnces
            for (unsigned long int i = 0; i < Diff_adatoms.size(); i++){
                //std :: cout << "\n thread "<< omp_get_thread_num() << "loop index" << i << std:: flush;

                x= std :: get<0>(Diff_adatoms[i]);
                y= std :: get<1>(Diff_adatoms[i]);
                
                i_rand = rand_r(&localseed[id]) % 4 +1;//independent seed for each thread, if seed does not change sequence keep going from that seed
               
               #pragma omp atomic
                    adatom.matrix[y][x] -= 1;
               
                if(i_rand ==1){      
                    int top = y+1;
                    if(top==L) top = 0;

                    #pragma omp atomic
                        adatom.matrix[top][x] += 1;

                }
                else if(i_rand ==2){

                    int bottom = y-1;
                    if(bottom==-1) bottom = L-1;

                   #pragma omp atomic
                        adatom.matrix[bottom][x] += 1;
                }
                else if(i_rand ==3){

                    int right = x+1;
                    if (right ==L) right = 0;

                    #pragma omp atomic
                        adatom.matrix[y][right] += 1;
                }
                else if(i_rand ==4){

                    int left = x-1;
                    if (left == -1) left = L-1;

                    #pragma omp atomic
                        adatom.matrix[y][left] += 1;
                }
            }

        // Refill diffusion and attachment classes 
        #pragma omp single  
        {
        R[diffusion].clear();
        R[attachment].clear();
        }
        
        
        #pragma omp for 
        for (int i = 0; i < L*L; i++){
            for(int k=0; k<adatom.matrix[i/L][i%L]; k++){//(y,x)[i][j] j advances at contiguous mem locations
                omp_set_lock(&writelock1);
                    R[diffusion].populate(i%L,i/L);
                omp_unset_lock(&writelock1);
                
                if(is_attSite(i%L,i/L)) {
                    omp_set_lock(&writelock2);
                        R[attachment].populate(i%L,i/L);
                    omp_unset_lock(&writelock2);
                }
            }
        }
    }
omp_destroy_lock(&writelock1);  
omp_destroy_lock(&writelock2); 

}

concentration = static_cast<double>(adatom.N)/(L*L);//update average concentration of adatoms

step_counter++;


//############# Time computation ########################

//d_ransetRate(((double) rand() / (RAND_MAX));
//double time = -log(d_rand)/R_sum;

//########################################################


if((proc_ID==root_process && debug_mode)||(proc_ID==root_process && error==true)){
    debug(who,x,y,error);
}


return 0;


}
 
 
 
void KMC :: debug (const unsigned who, const int x, const int y,  const bool error) const {

	std :: cout<< "\n  " << step_counter <<"  KMC step \n";
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

    else{
        std:: cout << "Nothing happenesetRate("<< who << "\n" << std :: flush ; 
        exit(EXIT_FAILURE);
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

    std :: cout << "\n diffusion list \n";
	 for (int i = 0; i < R[diffusion].N; i++){
	 	std :: cout << i <<"\t(" << R[diffusion].where(i)[0]<< ","<< R[diffusion].where(i)[1] << ")\n";
 	 }
 
if (error){
    std:: cout << "Error = "<< error << "\n" << std :: flush ; 
    exit(EXIT_FAILURE);
}

}