#include <cstdlib>
#include <iostream>
#include <vector>

#include "global.h"
#include "SPLASH.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SPLASH.cpp
 *
 * 2015-02-17 -- created
 * 2016-01-22 -- last updated
 *
 * ~~~~~~~~~
 * citation:
 * ~~~~~~~~~
 * T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
 * Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
 * algorithms for simulating habitats (SPLASH): Robust indices of radiation
 * evapotranspiration and plant-available moisture, Geoscientific Model
 * Development, 2016 (in progress)
 *
 * ~~~~~~~~~~~~
 * description:
 * ~~~~~~~~~~~~
 * This class updates daily quantities of radiation, evapotranspiration, soil
 * moisture and runoff based on the SPLASH methodology.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. added header to include list [15.02.19]
 * 02. added iostream to include list [15.02.19]
 * 03. created print_vals function [15.02.19]
 * 04. added smr struct (dsoil) [15.02.19]
 * 05. added vector to include list [15.02.19]
 * 06. created quick_run & spin_up functions [15.02.19]
 * 07. included global.h [16.01.22]
 *
 * //////////////////////////////////////////////////////////////////////// */

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Constructors:
   ///////////////////////////////////////////////////////////////////// */
SPLASH::SPLASH(double a, double b)
    : lat(a), elv(b), precip(0.0)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialize class & structure values:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    evap = EVAP(a, b);

    dvap.cond = 0.0;
    dvap.eet = 0.0;
    dvap.pet = 0.0;
    dvap.aet = 0.0;

    dsoil.sm = 0.0;
    dsoil.ro = 0.0;
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Function Definitions
   ///////////////////////////////////////////////////////////////////// */
void SPLASH::quick_run(int n, int y, double wn, double sf, double tc,
                       double pn, smr &dsm) {
    /* ***********************************************************************
    Name:     SPLASH.quick_run
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
    double sw = Global::Cw*(wn/Global::Wm);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate radiation and evaporation quantities
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    evap.calculate_daily_fluxes(sw, n, y, sf, tc);
    etr dn = evap.get_vals();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sm = wn + pn + dn.cond - dn.aet;

    // ~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    double ro = 0.0;
    if (sm > Global::Wm){
        // Bucket is too full:
        //   allocate excess water to runoff
        //   set soil moisture to capacity
        ro = (sm - Global::Wm);
        sm = Global::Wm;
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

void SPLASH::run_one_day(int n, int y, double wn, double sf, double tc,
                         double pn) {
    /* ***********************************************************************
    Name:     SPLASH.run_one_day
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
    double sw = Global::Cw*(wn/Global::Wm);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate radiation and evaporation quantities
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    evap.calculate_daily_fluxes(sw, n, y, sf, tc);
    dvap = evap.get_vals();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double sm = wn + pn + dvap.cond - dvap.aet;

    // ~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    double ro = 0.0;
    if (sm > Global::Wm){
        // Bucket is too full:
        //   allocate excess water to runoff
        //   set soil moisture to capacity
        ro = (sm - Global::Wm);
        sm = Global::Wm;
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

void SPLASH::spin_up(DATA &d){
    /* ***********************************************************************
    Name:     SPLASH.spin_up
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

double SPLASH::get_elv(){
    /* ***********************************************************************
    Name:     SPLASH.get_elv
    Input:    None
    Output:   double
    Features: Returns the elevation, meters.
    *********************************************************************** */
    return elv;
}

double SPLASH::get_lat(){
    /* ***********************************************************************
    Name:     SPLASH.get_lat
    Input:    None
    Output:   double
    Features: Returns the latitude, degrees.
    *********************************************************************** */
    return lat;
}

void SPLASH::print_vals(){
    /* ***********************************************************************
    Name:     SPLASH.print_vals
    Input:    None
    Output:   None
    Features: Prints the current dvap values.
    *********************************************************************** */
    evap.display();
    cout<<"Daily SPLASH values:" << endl;
    cout<<"  Cn: " << dvap.cond << " mm" << endl;
    cout<<"  EET: " << dvap.eet << " mm" << endl;
    cout<<"  PET: " << dvap.pet << " mm" << endl;
    cout<<"  AET: " << dvap.aet << " mm" << endl;
    cout<<"  Wn: " << dsoil.sm << " mm" << endl;
    cout<<"  RO: " << dsoil.ro << " mm" << endl;
}
