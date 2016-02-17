/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * etr.h
 *
 * VERSION 1.0
 * LAST UPDATED: 2016-02-16
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
 * Header file for the etr ('e'vapo't'ranspiration and 'r'adiation) struct.
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef ETR_H
#define ETR_H
struct etr {
    double cond;    // condensation, mm
    double eet;     // equilibrium evapotranspiration, mm
    double pet;     // potential evapotranspiration, mm
    double aet;     // actual evapotranspiration, mm
};
#endif
