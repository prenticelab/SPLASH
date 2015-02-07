#include "EVAP.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * main.cpp
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-07 -- created
 * 2015-02-07 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This script serves as the main function to run the STASH model.
 * 
 * //////////////////////////////////////////////////////////////////////// */

int main() {
    double lat = 51.4;
    int n = 172;
    double elv = 74.0; 
    int y = 2001; 
    float sf = 0.43; 
    float tc = 17.3; 
    float sw = 0.5;
    EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    my_evap.display();
}
