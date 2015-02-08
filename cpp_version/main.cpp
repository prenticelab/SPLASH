#include <iostream>
#include <string>
#include <vector>

#include "DATA.h"
#include "EVAP.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * main.cpp
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-07 -- created
 * 2015-02-08 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This script serves as the main function to run the STASH model.
 * 
 * ----------
 * changelog:
 * ----------
 * 01. added DATA class [15.02.08]
 * 
 * -----
 * todo:
 * -----
 * 01. create STASH class
 * 
 * //////////////////////////////////////////////////////////////////////// */

int main() {
    // Open input data:
    string fname = "example_data.csv";
    DATA my_data(fname);
    cout << "Read " << my_data.nlines() << " lines of data." << endl;
    
    // Calculate radiation & evaporation terms:
    double lat = 37.7;
    double elv = 142.0;
    
    int n = 172;
    int y = 2001; 
    double sf = my_data.get_one_sf(n-1); 
    double tc = my_data.get_one_tair(n-1); 
    double sw = 0.5;
    EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    my_evap.display();
}
