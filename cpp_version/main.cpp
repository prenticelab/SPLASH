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
 * 
 * -----
 * todo:
 * -----
 * 01. create STASH class
 * 
 * //////////////////////////////////////////////////////////////////////// */

int main() {
    // Open input data:
    DATA my_data;
    
    // Example 1: read all data from CSV file:
    string fname = "example_data.csv";
    my_data.read_csv(fname);
    //cout << "Read " << my_data.nlines() << " lines of data." << endl;
    //cout << "Year: " << my_data.get_year() << endl;
    
    // Example 2: read individual variable data from plain text files:
    //string sf_file = "daily_sf_2000_cruts.txt";
    //string pn_file = "daily_pn_2000_wfdei.txt";
    //string tc_file = "daily_tair_2000_wfdei.txt";
    //my_data.read_txt(pn_file, "pn", 2000);
    //my_data.read_txt(sf_file, "sf", 2000);
    //my_data.read_txt(tc_file, "tair", 2000);
    
    vector<double> my_sf = my_data.get_all_sf();
    vector<double> my_pn = my_data.get_all_pn();
    vector<double> my_tc = my_data.get_all_tair();
    
    cout << "Saved " << my_sf.size() << " data points to sf_vec" << endl;
    cout << "Saved " << my_pn.size() << " data points to pn_vec" << endl;
    cout << "Saved " << my_tc.size() << " data points to tair_vec" << endl;
    
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
