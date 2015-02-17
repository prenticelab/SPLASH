#include <cstdlib>
#include <vector>

#include "EVAP.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * STASH.cpp
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-17 -- created
 * 2015-02-17 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This class updates daily quantities of radiation, evapotranspiration, soil 
 * moisture and runoff based on the STASH methodology.
 * 
 * ----------
 * changelog:
 * ----------
 * 
 * 
 * //////////////////////////////////////////////////////////////////////// */

class STASH {
    private:
        // Constants:
        double kCw;      // supply constant, mm/hr 
        double kWm;      // soil moisture capacity, mm 
        
        // Static variables:
        double lat;                       // latitude, degrees
        double elv;                       // elevation, meters
        
        // Daily status variables:
        double ho;                        // daily solar irradiation, J/m2
        double hn;                        // daily net radiation, J/m2
        double ppfd;                      // daily PPFD, mol/m2
        double cond;                      // daily condensation water, mm
        double wn;                        // daily soil moisture, mm
        double precip;                    // daily precipitation, mm
        double ro;                        // daily runoff, mm
        double eet;                       // daily equilibrium ET, mm
        double pet;                       // daily potential ET, mm
        double aet;                       // daily actual ET, mm
        
    public:
        // Constructors:
        STASH(double latitude, double elevation);
        
        // Functions:
        void run_one_day(int n, int y, double wn, double sf, double tc, 
                         double pn);
        double get_lat();
        double get_elv();
};

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Constructors: 
// ////////////////////////////////////////////////////////////////////////
STASH::STASH(double latitude, double elevation)
    : lat(latitude), elv(elevation), ho(0.0), hn(0.0), ppfd(0.0), cond(0.0),
      wn(0.0), precip(0.0), ro(0.0), eet(0.0), pet(0.0), aet(0.0)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialize constants:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    kCw = 1.05;       // (Federer, 1982)
    kWm = 150.0;      // (Cramer & Prentice, 1988)
}

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Function Definitions
// ////////////////////////////////////////////////////////////////////////
void STASH::run_one_day(int n, int y, double wn, double sf, double tc, 
                        double pn) {
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 0. Set meteorological variable:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    precip = pn;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate evaporative supply rate (sw), mm/h
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sw = kCw*(wn/kWm);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate radiation and evaporation quantities
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    ho = my_evap.get_ho();
    hn = my_evap.get_hn();
    ppfd = my_evap.get_ppfd();
    cond = my_evap.get_cond();
    eet = my_evap.get
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sm = wn + pn + my_evap.get_cond() - my_evap.get_aet();
    
    // ~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    if (sm > kWm){
        // Bucket is too full:
        //   allocate excess water to runoff
        //   set soil moisture to capacity
        ro = (sm - kWm);
        sm = kWm;
    } else if (sm < 0){
        // Bucket is too empty:
        //   set soil moisture & runoff to zero
        ro = 0.0;
        sm = 0.0;
    }
}

double STASH::get_elv(){
    /* ***********************************************************************
    Name:     STASH.get_elv
    Input:    None
    Output:   double
    Features: Returns the elevation, meters.
    *********************************************************************** */
    return elv;
}

double STASH::get_lat(){
    /* ***********************************************************************
    Name:     STASH.get_lat
    Input:    None
    Output:   double
    Features: Returns the latitude, degrees.
    *********************************************************************** */
    return lat;
}
