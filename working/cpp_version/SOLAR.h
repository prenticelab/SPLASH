#include <vector>

#include "srad.h"

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SOLAR.h
 *
 * VERSION 1.1-dev
 * LAST UPDATED: 2016-02-19
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
 * This is the C++ header for the SOLAR class.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef SOLAR_H
#define SOLAR_H
class SOLAR {
    private:
        double lat;                // latitude, degrees
        double elv;                // elevation, m
        int day;                   // day of year
        int kN;                    // days in year
        double nu;                 // true anomaly, degrees
        double lam;                // true longitude, degrees
        double rho;                // rel. earth-sun distance
        double dr;                 // distance factor
        double delta;              // declination, degrees
        double ru, rv, rw, rx;     // variable substitutions
        double hs, hn, hi;         // hour angles, degrees
        double ra_d;               // daily solar irradiation, J/m^2
        double tau_o;              // surface transmittivity
        double tau;                // elv. corrected transmittivity
        double ppfd_d;             // daily PPFD, mol/m^2
        double rnl;                // longwave radiation flux, W/m^2
        double rn_d;               // daily net radiation, J/m^2
        double rnn_d;              // daily nighttime net radiation, J/m^2
        srad dsr;                  // daily srad struct

        // Functions:
        double dcos(double x);
        double dsin(double x);
        std::vector<double> berger_tls(int n);
        int julian_day(int y, int m, int i);

    public:
        // Constructors
        SOLAR();
        SOLAR(double a, double b);

        // General purpose functions:
        void calculate_daily_fluxes(int n, int y, double sf, double tc);
        srad get_vals();
        void display();
};
#endif
