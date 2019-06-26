#include "global.h"
#include "Events.h"
#include<list>
#include<tuple>
#include<algorithm>

void Events :: populate(int  x, int y){
    element.push_back(std::make_tuple(x,y));
   // index +=1;//useless
	mask[y][x] +=1; //keep trck of having more than 1 adatom on a site: useful only for adatom class

    N = element.size(); //= N of the class
    // std :: cout << "\n"<< index;
    // if(N!=index) N=0;
}
void Events :: change(int i,int x, int y){

    std :: list<std :: tuple<int, int>>::iterator it=element.begin();
    advance(it,i);
   // index -=1;
	mask[std :: get<1>(*it)][std :: get<0>(*it)] -=1;
    *it = (std::make_tuple(x,y));
    //N = element.size();//also updates N
    mask[y][x] +=1;
}
void Events :: destroy(int i){

    std :: list<std :: tuple<int, int>>::iterator it=element.begin();
    advance(it,i);
   // index -=1;
	mask[std :: get<1>(*it)][std :: get<0>(*it)] -=1;
    element.erase(it);
    N = element.size();//also updates N
}
void Events :: destroy_byPosition(int x, int y){
    int oldN;
    
    oldN = element.size();
    std :: tuple <int,int> coordinate;
    coordinate = std :: make_tuple(x,y);
    element.remove(coordinate);//O(n) Tkink better way?
	mask[std :: get<1>(coordinate)][std :: get<0>(coordinate)] -=1;

    N = element.size();//also updates N
    if (oldN == N) {
        std :: cout << "\n Unsucesful attempt to remove"<<" x = " << x <<" y = " << y << "unexisting element in class list. \n";
        exit (EXIT_FAILURE);
    }
}
bool Events :: exist (int x,int y){
    bool answ;
    std :: tuple <int,int> coordinate;
    coordinate = std :: make_tuple(x,y);
    answ = std::find(element.begin(), element.end(), coordinate) != element.end();
    return answ;
}



//method to return human readable coordinate

int * Events :: where (int i){
    std :: list<std :: tuple<int, int>>::iterator it=element.begin();
    static int coordinate[2];


    advance(it,i);

    coordinate[0] = std :: get<0>(*it);
    coordinate[1] = std :: get<1>(*it);
    return coordinate;
}


// method coordinate in


//---------------

double Events :: rate(){
    double r;
    r = D*double(N);
    return r;
}