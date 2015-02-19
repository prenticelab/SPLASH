#include "STASH.h"
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * STASH.cpp
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
 * This class updates daily quantities of radiation, evapotranspiration, soil 
 * moisture and runoff based on the STASH methodology.
 * 
 * ----------
 * changelog:
 * ----------
 * 01. added STASH header to include list [15.02.19]
 * 02. added iostream to include list [15.02.19]
 * 03. created print_vals function [15.02.19]
 * 04. added smr struct (dsoil) [15.02.19]
 * 05. added vector to include list [15.02.19]
 * 06. created quick_run & spin_up functions [15.02.19]
 * 
 * //////////////////////////////////////////////////////////////////////// */

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Constructors: 
// ////////////////////////////////////////////////////////////////////////
STASH::STASH(double latitude, double elevation)
    : lat(latitude), elv(elevation), precip(0.0)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialize constants:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    kCw = 1.05;       // (Federer, 1982)
    kWm = 150.0;      // (Cramer & Prentice, 1988)
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialize structure values:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dvap.ho = 0.0;
    dvap.hn = 0.0;
    dvap.ppfd = 0.0;
    dvap.cond = 0.0;
    dvap.eet = 0.0;
    dvap.pet = 0.0;
    dvap.aet = 0.0;
    
    dsoil.sm = 0.0;
    dsoil.ro = 0.0;
}

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Function Definitions
// ////////////////////////////////////////////////////////////////////////
void STASH::quick_run(int n, int y, double wn, double sf, double tc, 
                       double pn, smr &dsm) {
    /* ***********************************************************************
    Name:     STASH.quick_run
    Input:    - int, day of year (n)
              - int, year (y)
              - double, previous day's soil moisture, mm (wn)
              - double, daily sunshine fraction (sf)
              - double, daily air temperature, deg C (tc)
              - double, daily precipitation, mm (pn)
              - smr, daily soil moisture & runoff
    Output:   None.
    Features: Calculates daily soil moisture and runoff based on STASH 
              methods.
    *********************************************************************** */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate evaporative supply rate (sw), mm/h
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sw = kCw*(wn/kWm);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate radiation and evaporation quantities
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    etr dn = my_evap.get_vals();
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sm = wn + pn + dn.cond - dn.aet;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    double ro = 0.0;
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
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Update soil moisture & runoff
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dsm.sm = sm;
    dsm.ro = ro;
}

void STASH::run_one_day(int n, int y, double wn, double sf, double tc, 
                        double pn) {
    /* ***********************************************************************
    Name:     STASH.run_one_day
    Input:    - int, day of year (n)
              - int, year (y)
              - double, previous day's soil moisture, mm (wn)
              - double, daily sunshine fraction (sf)
              - double, daily air temperature, deg C (tc)
              - double, daily precipitation, mm (pn)
    Output:   None.
    Features: Calculates daily soil moisture and runoff based on STASH 
              methods and updates class variables accordingly.
    *********************************************************************** */
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
    dvap = my_evap.get_vals();
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sm = wn + pn + dvap.cond - dvap.aet;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    double ro = 0.0;
    if (sm > kWm){
        // Bucket is too full:
        //   allocate excess water to runoff
        //   set soil moisture to capacity
        ro = (sm - kWm);
        sm = kWm;
    } else if (sm < 0){
        // Bucket is too empty:
        //   reduce actual ET by discrepany amout
        //   set soil moisture & runoff to zero
        dvap.aet += sm;
        ro = 0.0;
        sm = 0.0;
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Update soil moisture & runoff
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dsoil.sm = sm;
    dsoil.ro = ro;
}

void STASH::spin_up(DATA &d){
    /* ***********************************************************************
    Name:     STASH.spin_up
    Input:    DATA class (d)
    Output:   None
    Features: Spins up the daily soil moisture
    Depends:  quick_run
    *********************************************************************** */
    double wn;  // previous day's soil moisture, mm
    smr dsm;    // daily soil moisture and runoff
    
    // Create daily soil moisture vector
    int n = d.nlines();
    int y = d.get_year();
    vector<double> wn_vec(n, 0.0);
    
    // Run one year:
    for (int i=0; i<n; i++){
        // Get preceeding soil moisture status:
        if (i == 0){
            wn = wn_vec[n-1];
        } else {
            wn = wn_vec[(i-1)];
        }
        
        // Calculate soil moisture and runoff
        quick_run((i+1), y, wn, d.get_one_sf(i), d.get_one_tair(i), 
                  d.get_one_pn(i), dsm);
        
        wn_vec[i] = dsm.sm;
    }
    
    // Calculate change in starting soil moisture:
    double start_sm = wn_vec[0];
    quick_run(1, y, wn_vec[n-1], d.get_one_sf(0), d.get_one_tair(0), 
              d.get_one_pn(0), dsm);
    double end_sm = dsm.sm;
    double diff_sm = (end_sm - start_sm);
    if (diff_sm < 0){
        diff_sm = (start_sm - end_sm);
    }
    
    // Equilibrate
    int spin_count = 1;
    while (diff_sm > 1.0){
        for (int i=0; i<n; i++){
            // Get preceeding soil moisture status:
            if (i == 0){
                wn = wn_vec[n-1];
            } else {
                wn = wn_vec[(i-1)];
            }
            
            // Calculate soil moisture and runoff
            quick_run((i+1), y, wn, d.get_one_sf(i), d.get_one_tair(i), 
                      d.get_one_pn(i), dsm);
            
            wn_vec[i] = dsm.sm;
        }
        
        // Calculate difference
        start_sm = wn_vec[0];
        quick_run(1, y, wn_vec[n-1], d.get_one_sf(0), d.get_one_tair(0), 
              d.get_one_pn(0), dsm);
        end_sm = dsm.sm;
        diff_sm = (end_sm - start_sm);
        if (diff_sm < 0){
            diff_sm = (start_sm - end_sm);
        }
        
        spin_count++;
    }
    
    // Save initial soil moisture condition:
    dsoil.sm = wn_vec[n-1];
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

void STASH::print_vals(){
    /* ***********************************************************************
    Name:     STASH.print_vals
    Input:    None
    Output:   None
    Features: Prints the current dvap values.
    *********************************************************************** */
    cout<<"Daily values:" << endl;
    cout<<"  Ho: " << (1.0e-6)*dvap.ho << " MJ/m^2" << endl;
    cout<<"  HN: " << (1.0e-6)*dvap.hn << " MJ/m^2" << endl;
    cout<<"  PAR: " << dvap.ppfd << " mol/m^2" << endl;
    cout<<"  Cn: " << dvap.cond << " mm" << endl;
    cout<<"  EET: " << dvap.eet << " mm" << endl;
    cout<<"  PET: " << dvap.pet << " mm" << endl;
    cout<<"  AET: " << dvap.aet << " mm" << endl;
    cout<<"  Wn: " << dsoil.sm << " mm" << endl;
    cout<<"  RO: " << dsoil.ro << " mm" << endl;
}
