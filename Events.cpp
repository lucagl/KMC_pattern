//#include "global.h"
#include "Events.h"
#include "functions.h"


void Events :: populate(const int  x,const int y){
    element.push_back(std::make_tuple(x,y));
   // index +=1;//useless
	//mask[y][x] +=1; //keep trck of having more than 1 adatom on a site: useful only for adatom class

    N = element.size(); //= N of the class
    // std :: cout << "\n"<< index;
    // if(N!=index) N=0;
}
void Events :: change(const int i,const int x, const int y){

    auto it=element.begin();
    advance(it,i);
   // index -=1;
	//mask[std :: get<1>(*it)][std :: get<0>(*it)] -=1;
    *it = (std::make_tuple(x,y));

  //  mask[y][x] +=1;
}
void Events :: destroy(const int i){

    auto it=element.begin();
    advance(it,i);
   // index -=1;
	//mask[std :: get<1>(*it)][std :: get<0>(*it)] -=1;
    element.erase(it);
    N = element.size();// updates N
}
bool Events :: destroy_coordinates(const int x, const int y){
    bool error =0;
    int oldN;
    
    oldN = element.size();
    auto coordinate = std :: make_tuple(x,y);
    element.remove(coordinate);//O(n) Tkink better way?
	//mask[std :: get<1>(coordinate)][std :: get<0>(coordinate)] =0;

    N = element.size();// updates N


    if (oldN == N) {
        std :: cout << "\n Old N =" <<oldN << "\t" << "new N =" << N;  
        std :: cout << "\n Unsucesful attempt to remove"<<" x = " << x <<" y = " << y << "unexisting MULTIPLE OR SINGLE element in class list. \n"<< std::flush ;
        //exit (EXIT_FAILURE);
        error=true;
    }
    return error;
}
bool Events :: destroy_singleCoordinate(const int x, const int y){
    bool error =0;
    
    
    auto coordinate = std :: make_tuple(x,y);
    auto it = std::find(element.begin(), element.end(), coordinate);

    if (it == element.end()) {
        std :: cout << "\n Unsucesful attempt to remove"<<" x = " << x <<" y = " << y << "unexisting SINGLE element in class list. \n"<< std::flush ;
        // exit (EXIT_FAILURE);
        error = true;
    }
    element.erase(it);
	//mask[std :: get<1>(coordinate)][std :: get<0>(coordinate)] -=1;

    N = element.size();// updates N
    return error;
}
bool Events :: exist (const int x,const int y) const{
    // return true if coordinate already present. More than one could be present
    bool answ;
    
    auto coordinate = std :: make_tuple(x,y);
    answ = std::find(element.begin(), element.end(), coordinate) != element.end();
    return answ;
}



//method to return human readable coordinate

int * Events :: where (const int i) const{

    auto it=element.begin();

    static int coordinate[2];

    advance(it,i);

    coordinate[0] = std :: get<0>(*it);
    coordinate[1] = std :: get<1>(*it);
    return coordinate;
}
// void Events :: clear () {
//     element.clear();    
//     N=0;
// }




void Events :: setRateTimeDep(const double T_new , const double T_old) {
    r = pow(r, T_old/T_new);
    return;
}