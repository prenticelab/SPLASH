#include <iostream>
#include <string>
#include <vector>

#include "DATA.h"
#include "EVAP.h"
#include "STASH.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * main.cpp
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-07 -- created
 * 2015-02-17 -- last updated
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
 * 02. updated DATA class for csv and txt data input [15.02.17]
 * 03. added STASH class [15.02.17]
 * 
 * -----
 * todo:
 * -----
 * 01. create STASH class
 * 
 * //////////////////////////////////////////////////////////////////////// */

int main() {
    // Example 1: read all data from CSV file:
    DATA my_data;
    string fname = "example_data.csv";
    my_data.read_csv(fname);
    
    // Example 2: read individual variables from plain text files:
    //DATA my_data;
    //string sf_file = "daily_sf_2000_cruts.txt";
    //string pn_file = "daily_pn_2000_wfdei.txt";
    //string tc_file = "daily_tair_2000_wfdei.txt";
    //my_data.read_txt(pn_file, "pn", 2000);
    //my_data.read_txt(sf_file, "sf", 2000);
    //my_data.read_txt(tc_file, "tair", 2000);
    
    vector<double> my_sf = my_data.get_all_sf();
    vector<double> my_pn = my_data.get_all_pn();
    vector<double> my_tc = my_data.get_all_tair();
    
    cout << "Year: " << my_data.get_year() << endl;
    cout << "Saved " << my_sf.size() << " data points to sf_vec" << endl;
    cout << "Saved " << my_pn.size() << " data points to pn_vec" << endl;
    cout << "Saved " << my_tc.size() << " data points to tair_vec" << endl;
    
    // Calculate radiation & evaporation terms:
    double lat = 37.7;
    double elv = 142.0;
    STASH my_stash(lat, elv);
    
    cout << "Created STASH class at " << my_stash.get_elv() << " m" << endl;
    
    //int n = 172;
    //int y = 2001; 
    //double sf = my_data.get_one_sf(n-1); 
    //double tc = my_data.get_one_tair(n-1);
    //double sw = 0.5;
    //EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    //my_evap.display();
}
