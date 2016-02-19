#include <vector>

#include "etr.h"
#include "smr.h"
#include "DATA.h"
#include "EVAP.h"

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SPLASH.h
 *
 * VERSION 1.1-dev
 * LAST UPDATED: 2016-02-17
 *
 * ~~~~~~~~
 * license:
 * ~~~~~~~~
 * Copyright (C) 2015, see LICENSE
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
 * This is the header file for the C++ SPLASH class.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. moved header includes from cpp file here [15.02.19]
 * 02. created header guard [15.02.19]
 * 03. added etr structure [15.02.19]
 * 04. added smr structure [15.02.19]
 * 05. added DATA header to include list [15.02.19]
 * 06. added quick_run & spin_up functions [15.02.19]
 * 07. removed constants; now in global.h [16.01.22]
 * 08. made wn_vec a private variable [16.02.06]
 * 09. added print daily wn function [16.02.06]
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef SPLASH_H
#define SPLASH_H
class SPLASH {
    private:
        // Static variables:
        double lat;                       // latitude, degrees
        double elv;                       // elevation, meters

        // Daily status variables:
        EVAP evap;                        // daily EVAP class
        etr dvap;                         // daily etr struct
        smr dsoil;                        // daily smr struct
        double precip;                    // daily precipitation, mm

        // Daily soil moisture
        std::vector<double> wn_vec;

    public:
        // Constructors:
        SPLASH(double a, double b);

        // Functions:
        void quick_run(int n, int y, double wn, double sf, double tc,
                       double pn, smr &dsm);
        void run_one_day(int n, int y, double wn, double sf, double tc,
                         double pn);
        void spin_up(DATA &d);
        double get_lat();
        double get_elv();
        void print_daily_wn();
        void print_vals();
};
#endif
