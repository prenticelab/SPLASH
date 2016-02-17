/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * srad.h
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
 * Header file for the radiation struct, srad.
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef SRAD_H
#define SRAD_H
struct srad {
    double ru;      // variable substitute
    double rv;      // variable substitute
    double rw;      // variable substitute
    double rnl;     // net longwave radiation, W/m^2
    double hn;      // cross-over hour angle, degrees
    double rn_d;    // daytime net radiation, J/m^2
    double rnn_d;   // nighttime net radiation, J/m^2
};
#endif
