/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * smr.h
 * 
 * 2015-02-19 -- created
 * 2015-08-22 -- last updated
 * 
 * ~~~~~~~~~
 * citation:
 * ~~~~~~~~~
 * T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
 * Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
 * algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
 * transpiration and plant-available moisture, Geoscientific Model Development, 
 * 2015 (in progress)
 * 
 * ~~~~~~~~~~~~
 * description:
 * ~~~~~~~~~~~~
 * Header file for the smr ('s'oil 'm'oisture and 'r'unoff) struct.
 * 
 * ~~~~~
 * todo:
 * ~~~~~
 * 1. Consider adding an index (for day or month of year)
 * 2. Consider adding precip and AET for annual water balance checks, i.e.:
 *    SUM_yr(Pn + Cn) = SUM_yr(AET + RO) --- or make another struct for this?
 * 
 * //////////////////////////////////////////////////////////////////////// */
#ifndef SMR_H
#define SMR_H
struct smr {
    double sm;      // soil moisture, mm
    double ro;      // runoff, mm
};
#endif
