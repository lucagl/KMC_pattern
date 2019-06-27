// Observable Class contains relevant  observables
#ifndef EVENTS_H
#define EVENTS_H

#include "global.h"
#include<list>
#include<tuple>

class Events {
	

	private:
		int L;
		//int index =0;// just to check, can be erased in funture versions
		std::list<std::tuple<int, int>> element{}; // list of tuple (i-th element touple of coordinates)
	//list more efficient than vector for deletion and insertion

	public:
		double D; // constant rate
		std::tuple<int, int> coordinate(int i); //return a tuple
		int N=0; // elements in the class
		int ** mask; //better to make it boolean
		void coordinate_in (int i, int x, int y); //change coordinate element, to keep track changes. Maybe input a touple
		void destroy(int i); //must take care of reduction of elements. When called mask[][] must also set to 0 concerned element
		void destroy_coordinates (int x,int y);
		void destroy_singleCoordinate (int x,int y);
		void populate (int x, int y);
		void change (int i, int x, int y);
		bool exist (int x, int y);
		int* where (int i);
		double rate();

		Events (int box_size){

		//L is the box size, to allocate the mask matrix
			L = box_size;
			mask = new int*[L];
		
			for (int i =0;i<L;i++){
				mask[i] = new int[L] (); 
			}



		}





		~Events(){
			element.clear();
			for(int i = 0; i < L; ++i) {
				delete [] mask[i];
			}
		delete [] mask;
		}
};

#endif
