#ifndef	KMC_H
#define KMC_H

#include "Events.h"
#include "Island.h"
#include "Adatom.h"


const int n_classes = 26 ;

class KMC {

public:

    KMC(const double, const double, const double, const double,const int, const bool, const int, const double, const double);

    void initConv_island (double sigma){island.initConv(sigma);};
    void initConv_adatom (double sigma){adatom.initConv(sigma);};

    void reset ();

    void init();

    void read(const std :: string);

    double step(const double , const bool debug_mode  = false);

    void saveTxt(const std :: string,const int, const double , const bool isConv=false,const  bool flag =false ) const;
    void print_final(const int, const bool) ;//CHECK and CHANGE not very elegant

    double get_concentration() const;
    int* get_nevents() const;
    int* get_classN() const;

    double c0; //initial concentration

    unsigned short** getIsland() const { 
        unsigned short** matrix;//like this I return by copy and not the pointer: this is safe
        matrix = new unsigned short*[L];
        for (int i = 0; i < L; i++)
        {
            matrix[i] = new unsigned short[L];
            for (int j = 0; j < L; j++)
            {
                matrix[i][j] = island.matrix[i][j];
            }
            
        }
        return matrix;
    };
    //unsigned** getIsland_conv(double sigma) const { return &island.gaussianConv(signma);};
    unsigned short ** getAdatom()const {
        unsigned short** matrix;//like this I return by copy and not the pointer: this is safe
        matrix = new unsigned short *[L];
        for (int i = 0; i < L; i++)
        {
            matrix[i] = new unsigned short [L];
            for (int j = 0; j < L; j++)
            {
                matrix[i][j] = adatom.matrix[i][j];
            }
            
        }
        return matrix;
    };

    double ** getIslandConv()const{
        return island.gaussianConv();
    };
    double ** getAdatomConv()const{
        return adatom.gaussianConv();
    };

private: 
    bool init_isCircle;
    int init_radius;
    int event_counter[n_classes] {};
    int L;
    double  concentration;
    double J;
    double BR;
    double A;//change to treat it as the temperature ( a parameter passed from main)
    double E_shift;
    double  current_T;

    Adatom adatom;
    Island island;
    Events R[n_classes];
    
    double det_rate(const int, const int ) const;
    double att_rate() const;
    double cumulative ( double * ) ;// calls rate evaluation, so is not a constant member function
    bool is_attSite(const int , const int ) const;
    void debug(const unsigned, const int, const int, const bool) const;

    bool update_nn1DetachmentClasses(const int, const int);
    bool update_nn2DetachmentClasses(const int, const int);
    bool update_AttachmentClasses(const int, const int);

};


#endif
