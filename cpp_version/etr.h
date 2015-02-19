/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * etr.h
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
