#include <iostream>
#include <string>
#include <vector>

#include "DATA.h"
#include "EVAP.h"
#include "SPLASH.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * main.cpp
 *
 * VERSION 1.1-dev
 * LAST UPDATED: 2016-09-01
 *
 * ~~~~~~~~
 * license:
 * ~~~~~~~~
 * Copyright (C) 2016 Prentice Lab
 *
 * This file is part of the SPLASH model.
 *
 * SPLASH is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * SPLASH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
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
 * This script serves as the main function to run the SPLASH model.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. added DATA class [15.02.08]
 * 02. updated DATA class for csv and txt data input [15.02.17]
 * 03. added SPLASH class [15.02.17]
 * 04. added test 3 [16.02.05]
 *
 * //////////////////////////////////////////////////////////////////////// */

int main() {
    // Test 1 & 2: SOLAR and EVAP classes:
    EVAP my_evap = EVAP(37.7, 142.0);
    my_evap.calculate_daily_fluxes(0.9, 172, 2000, 1.0, 23.0);
    cout << "Test 1 & 2: " << endl;
    my_evap.display();

    // Test 3: SPLASH run_one_day
    SPLASH my_splash = SPLASH(37.7, 142.0);
    my_splash.run_one_day(172, 2000, 75, 1.0, 23.0, 5.0);
    cout << "Test 3: " << endl;
    my_splash.print_vals();

    // Test 4: SPIN-UP
    DATA my_data;
    string fname = "../../data/example_data.csv";
    my_data.read_csv(fname);
    double lat = 37.7;
    double elv = 142.0;
    SPLASH my_class(lat, elv);
    my_class.spin_up(my_data);
    my_class.print_daily_wn();

    /*

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
    SPLASH my_class(lat, elv);
    my_class.spin_up(my_data);
    // run_one_day(int n, int y, double wn, double sf, double tc, double pn)
    my_class.run_one_day(172, 2000, 145.0, 1.0, 23.0, 10.0);
    my_class.print_vals();
    //cout << "Created SPLASH class at " << my_class.get_elv() << " m" << endl;

    //int n = 172;
    //int y = 2001;
    //double sf = my_data.get_one_sf(n-1);
    //double tc = my_data.get_one_tair(n-1);
    //double sw = 0.5;
    //EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    //my_evap.display();

    */
}
