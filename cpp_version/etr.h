/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * etr.h
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
 * Header file for the etr ('e'vapo't'ranspiration and 'r'adiation) struct.
 * 
 * //////////////////////////////////////////////////////////////////////// */
#ifndef ETR_H
#define ETR_H
struct etr {
    double ho;      // solar irradiation, J/m^2
    double hn;      // net radiation, J/m^2
    double ppfd;    // PPFD, mol/m^2
    double cond;    // condensation, mm
    double eet;     // equilibrium evapotranspiration, mm
    double pet;     // potential evapotranspiration, mm
    double aet;     // actual evapotranspiration, mm
};
#endif
