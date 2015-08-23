#include <vector>

#include "DATA.h"
#include "EVAP.h"
#include "etr.h"
#include "smr.h"

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SPLASH.h
 * 
 * 2015-02-17 -- created
 * 2015-08-22 -- last updated
 * 
 * ~~~~~~~~~
 * citation:
 * ~~~~~~~~~
 * T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
 * Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
 * algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
 * transpiration and plant-available moisture, Geoscientific Model Development, 
 * 2015 (in progress)
 * 
 * ~~~~~~~~~~~~
 * description:
 * ~~~~~~~~~~~~
 * This is the header file for the C++ SPLASH class.
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
#ifndef SPLASH_H
#define SPLASH_H
class SPLASH {
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
        SPLASH(double latitude, double elevation);
        
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
