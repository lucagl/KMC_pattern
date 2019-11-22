#ifndef	KMC_H
#define KMC_H
#include "Events.h"
#include "Island.h"
#include "Adatom.h"


const int n_classes = 26 ;

class KMC {

public:

    KMC(const double, const double, const double, const double,const int, const bool, const int, const double, const double T0);

    void reset ();

    void init();

    void read(const std :: string);

    double step(const double , const bool debug_mode  = false);

    void saveTxt(const std :: string,const int,const  int flag =0 ) const;
    void print_final(int) const;//CHECK and CHANGE

    double get_concentration() const;
    int* get_nevents() const;
    int* get_classN() const;
    double c0; //initial concentration

    unsigned** getIsland();
    unsigned** getAdatom();


private: 
    bool init_isCircle;
    int init_radius;
    int event_counter[n_classes] {};
    int L;
    double  concentration;
    double J;
    double BR;
    double A;//change to treat it as the temperature
    double E_shift;
    double  current_T;

    Adatom adatom;
    Island island;
    Events R[n_classes];
    
    double det_rate(const int, const int ) const;
    double att_rate() const;
    double cumulative ( double * ) ;// calls rate evaluation, so is not a constant member function
    int extract (const int ) const;// extracts event
    bool is_attSite(const int , const int ) const;
    void debug(const unsigned, const int, const int, const bool) const;

    bool update_nn1DetachmentClasses(const int, const int);
    bool update_nn2DetachmentClasses(const int, const int);
    bool update_AttachmentClasses(const int, const int);

};


#endif
