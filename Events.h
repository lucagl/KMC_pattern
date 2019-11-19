// Observable Class contains relevant  observables
#ifndef EVENTS_H
#define EVENTS_H

#include "global.h"
#include<list>
#include<vector>
#include<tuple>

class Events {
	

	private:
		int L;
		std::list<std::tuple<int, int>> element{}; // list of tuple (i-th element touple of coordinates)
	//list more efficient than vector for deletion and insertion
		bool is_diff;
		double rate_diffusion();
		double rate_general();
		double (Events :: *tRate) ();

	public:
		double D; // constant rate
		std::tuple<int, int> coordinate(int i); //return a tuple
		unsigned long int N=0; // elements in the class
		//unsigned short ** mask; 
		
		void destroy(int i); //destroy ith element
		bool destroy_coordinates (const int ,const int );//destroy ALL elements with that coordinate
		bool  destroy_singleCoordinate (const int , const int ); // destroy only one occurrence of that coordinate
		void populate (const int , const int );
		void change (const int , const int , const int );// could be useful: changes coordinate of a specific element in the list
		bool exist (const int ,const int ) const;// checks existence of the element in the list based on its coordinate
		int* where (const int ) const; //returns coordinate of the element
		void clear();
		//double rate();
		
		double getRate() {return (this->*tRate) ();}

		void init (const int box_size, const bool isDiff =0){
			L = box_size;
			is_diff = isDiff;
			if(is_diff) tRate = & Events ::  rate_diffusion;
			else tRate = & Events :: rate_general;
			//mask = new unsigned short*[L];
			// for (int i =0;i<L;i++){
			// 	mask[i] = new unsigned short[L] (); 
			// }
		}
		~Events(){
			element.clear();
		// 	for(int i = 0; i < L; ++i) {
		// 		delete [] mask[i];
		// 	}
		// delete [] mask;
		}
		};


#endif
