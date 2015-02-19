#include <vector>

#include "DATA.h"
#include "EVAP.h"
#include "etr.h"
#include "smr.h"

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * STASH.h
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-17 -- created
 * 2015-02-19 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This is the header file for the C++ STASH class.
 * 
 * ----------
 * changelog:
 * ----------
 * 01. moved header includes from cpp file here [15.02.19]
 * 02. created header guard [15.02.19]
 * 03. added etr structure [15.02.19]
 * 04. added smr structure [15.02.19]
 * 05. added DATA header to include list [15.02.19]
 * 06. added quick_run & spin_up functions [15.02.19]
 * 
 * //////////////////////////////////////////////////////////////////////// */
#ifndef STASH_H
#define STASH_H
class STASH {
    private:
        // Constants:
        double kCw;      // supply constant, mm/hr 
        double kWm;      // soil moisture capacity, mm 
        
        // Static variables:
        double lat;                       // latitude, degrees
        double elv;                       // elevation, meters
        
        // Daily status variables:
        etr dvap;                         // daily etr struct
        smr dsoil;                        // daily smr struct
        double precip;                    // daily precipitation, mm
        
    public:
        // Constructors:
        STASH(double latitude, double elevation);
        
        // Functions:
        void quick_run(int n, int y, double wn, double sf, double tc, 
                       double pn, smr &dsm); 
        void run_one_day(int n, int y, double wn, double sf, double tc, 
                         double pn);
        void spin_up(DATA &d);
        double get_lat();
        double get_elv();
        void print_vals();
};
#endif
