/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * smr.h
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-19 -- created
 * 2015-02-19 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * Header file for the smr ('s'oil 'm'oisture and 'r'unoff) struct.
 * 
 * -----
 * todo:
 * -----
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
