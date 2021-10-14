// Observable Class contains relevant  observables
#ifndef EVENTS_H
#define EVENTS_H

#include "functions.h"

class Events {
	

	private:
		std::list<std::tuple<int, int>> element{}; // list of tuple (i-th element touple of coordinates)
	//list more efficient than vector for deletion and insertion
		bool is_diff;
		double r;
		inline double rate_general(){return r*double(N);};
		double rate_diffusion(){return r;}; 
		double (Events :: *tRate) ();
		 // elements in the class
		

	public:
		void setRate(double in_rate){r=in_rate;}; // constant rate
		double getPeratomRate(){return r;};
		//std::tuple<int, int> coordinate(int i); //return a tuple
		//unsigned getNumber() const {return N;};
		int  N=0;
		//unsigned short ** mask; 
		
		void destroy(int i); //destroy ith element
		bool destroy_coordinates (const int ,const int );//destroy ALL elements with that coordinate
		bool  destroy_singleCoordinate (const int , const int ); // destroy only one occurrence of that coordinate
		void populate (const int , const int );
		void change (const int , const int , const int );// could be useful: changes coordinate of a specific element in the list
		bool exist (const int ,const int ) const;// checks existence of the element in the list based on its coordinate
		int* where (const int ) const; //returns coordinate of the element
		void clear() {element.clear(); N=0;};
		std::list<std::tuple<int, int>> :: iterator begin() {return element.begin();};
		std::list<std::tuple<int, int>> :: iterator end() {return element.end();};
		std::list<std::tuple<int, int>> :: iterator erase(std::list<std::tuple<int, int>> :: iterator i){return element.erase(i);};
		int size() {N=element.size(); return N;};
		void merge(Events& other) {N+=other.element.size(); element.merge(other.element);};
		//double rate();
		
		double getRate() {return (this->*tRate) ();}

		void init (const bool isDiff =0){
			//L = box_size;
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
		// Events(){};
		Events& operator=( const Events& other ) {
		N = other.N;
		element = other.element;
		// NOTE : keep original def of rate not the one from the assigned event class..
		if(is_diff) tRate = & Events ::  rate_diffusion;
		else tRate = & Events :: rate_general;
		
		return *this;
  }
		};

// Events& Events::operator=( const Events& other ) {
//       N = other.N;
//       element = other.element;
//       return *this;
//   }
#endif
