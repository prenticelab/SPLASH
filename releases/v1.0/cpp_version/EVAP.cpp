#include <cmath>
#include <stdio.h>

#include "global.h"
#include "EVAP.h"
#include "SOLAR.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * EVAP.cpp
 *
 * VERSION 1.0-r2
 * LAST UPDATED: 2016-09-11
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
 * This script includes the definitions for the EVAP class, which calculates
 * the daily quantities of radiation, evaporation, and condensation for the
 * SPLASH model.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. finished class function definitions and debugged [15.02.07]
 * 02. removed kCw and kWm from constants (moved to SPLASH class) [15.02.17]
 * 03. changed get_cn to get_cond [15.02.17]
 * 04. added EVAP header file to include list [15.02.19]
 * 05. updated R and To [15.08.22]
 * 06. included global.h [16.01.22]
 * 07. updated documentation [16.05.27]
 * 08. addressed specific heat limitation [16.09.11]
 *
 * //////////////////////////////////////////////////////////////////////// */

/*
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Class Constructors:
/////////////////////////////////////////////////////////////////////
*/
EVAP::EVAP()
    : lat(0.0), elv(0.0), solar(SOLAR())
{
}

EVAP::EVAP(double a, double b){
    lat = a;
    elv = b;
    solar = SOLAR(lat, elv);
}

void EVAP::calculate_daily_fluxes(double sw, int n, int y, double sf,
                                  double tc){
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate radiation values
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    solar.calculate_daily_fluxes(n, y, sf, tc);
    d_sr = solar.get_vals();
    ru = d_sr.ru;
    rv = d_sr.rv;
    rw = d_sr.rw;
    rnl = d_sr.rnl;
    hn = d_sr.hn;
    rn_d = d_sr.rn_d;
    rnn_d = d_sr.rnn_d;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate water-to-energy conversion (econ), m^3/J
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    patm = elv2pres(elv);
    s = sat_slope(tc);
    lv = enthalpy_vap(tc);
    pw = density_h2o(tc, patm);
    g = psychro(tc, patm);
    econ = s/(lv*pw*(s + g));

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate daily condensation (wc), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cn = (1.0e3)*econ*abs(rnn_d);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Estimate daily equilibrium evapotranspiration (eet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eet_d = (1.0e3)*econ*rn_d;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Estimate daily potential evapotranspiration (pet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pet_d = (1.0 + Global::w)*eet_d;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx = (3.6e6)*(1.0 + Global::w)*econ;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7. Calculate the intersection hour angle (hi), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cos_hi = sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv;
    if (cos_hi >= 1.0) {
        // Supply exceeds demand:
        hi = 0.0;
    } else if (cos_hi <= -1.0) {
        // Supply limits demand everywhere:
        hi = 180.0;
    } else {
        hi = acos(cos_hi);
        hi /= Global::pir;
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 8. Estimate daily actual evapotranspiration (aet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aet_d = sw*hi*Global::pir;
    aet_d += rx*rw*rv*(dsin(hn) - dsin(hi));
    aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*Global::pir;
    aet_d *= (24.0/Global::PI);
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Function Definitions
   ///////////////////////////////////////////////////////////////////// */
double EVAP::dcos(double x){
    /* ***********************************************************************
    Name:     EVAP.dcos
    Input:    double, angle, degrees (x)
    Output:   double, cos(x*pi/180)
    Features: Calculates the cosine of an angle given in degrees
    *********************************************************************** */
    double v = cos(x*Global::pir);
    return v;
}

double EVAP::dsin(double x){
    /* ***********************************************************************
    Name:     EVAP.dsin
    Input:    double, angle, degrees (x)
    Output:   double, sin(x*pi/180)
    Features: Calculates the sine of an angle given in degrees
    *********************************************************************** */
    double v = sin(x*Global::pir);
    return v;
}

double EVAP::sat_slope(double tc){
    /* ***********************************************************************
    Name:     EVAP.sat_slope
    Input:    double, air temperature (tc), degrees C
    Output:   double, slope of sat vap press temp curve (s)
    Features: Calculates the slope of the sat pressure temp curve, Pa/K
    Ref:      Eq. 13, Allen et al. (1998)
    *********************************************************************** */
    double s = exp((tc*17.269)/(tc + 237.3));
    s /= pow((tc + 237.3), 2.0);
    s *= (17.269)*(237.3)*(610.78);
    return s;
}

double EVAP::enthalpy_vap(double tc){
    /* ***********************************************************************
    Name:     EVAP.enthalpy_vap
    Input:    double, air temperature (tc), degrees C
    Output:   double, latent heat of vaporization
    Features: Calculates the enthalpy of vaporization, J/kg
    Ref:      Eq. 8, Henderson-Sellers (1984)
    *********************************************************************** */
    double lv = (tc + 273.15)/(tc + 273.15 - 33.91);
    lv = pow(lv, 2.0);
    lv *= 1.91846e6;
    return lv;
}

double EVAP::elv2pres(double z){
    /* ***********************************************************************
    Name:     EVAP.elv2pres
    Input:    double, elevation above sea level (z), m
    Output:   double, atmospheric pressure, Pa
    Features: Calculates atm. pressure for a given elevation
    Depends:  Global constants
              - Po .... base pressure     - To .... base temperature
              - L ..... temp. lapse rate  - Ma .... molecular wt. dry air
              - G ..... standard gravity  - R ..... universal gas constant
    Ref:      Allen et al. (1998)
    *********************************************************************** */
    double ep = (Global::G * Global::Ma)/(Global::R * Global::L);
    double pa = (1.0 - z*Global::L/Global::To);
    pa = pow(pa, ep);
    pa *= Global::Po;
    return pa;
}

double EVAP::density_h2o(double tc, double p){
    /* ***********************************************************************
    Name:     EVAP.density_h2o
    Input:    - double, air temperature (tc), degrees C
              - double, atmospheric pressure (p), Pa
    Output:   double, density of water, kg/m^3
    Features: Calculates density of water at a given temperature and
              pressure
    Ref:      Chen et al. (1977)
    *********************************************************************** */
    // Calculate density at 1 atm:
    double po = 0.99983952;
    po += (6.788260e-5)*tc;
    po += -(9.08659e-6)*tc*tc;
    po += (1.022130e-7)*tc*tc*tc;
    po += -(1.35439e-9)*tc*tc*tc*tc;
    po += (1.471150e-11)*tc*tc*tc*tc*tc;
    po += -(1.11663e-13)*tc*tc*tc*tc*tc*tc;
    po += (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc;
    po += -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc;

    // Calculate bulk modulus at 1 atm:
    double ko = 19652.17;
    ko += 148.1830*tc;
    ko += -2.29995*tc*tc;
    ko += 0.01281*tc*tc*tc;
    ko += -(4.91564e-5)*tc*tc*tc*tc;
    ko += (1.035530e-7)*tc*tc*tc*tc*tc;

    // Calculate temperature dependent coefficients:
    double ca = 3.26138;
    ca += (5.223e-4)*tc;
    ca += (1.324e-4)*tc*tc;
    ca += -(7.655e-7)*tc*tc*tc;
    ca += (8.584e-10)*tc*tc*tc*tc;

    double cb = (7.2061e-5);
    cb += -(5.8948e-6)*tc;
    cb += (8.69900e-8)*tc*tc;
    cb += -(1.0100e-9)*tc*tc*tc;
    cb += (4.3220e-12)*tc*tc*tc*tc;

    // Convert atmospheric pressure to bar (1 bar = 100000 Pa)
    double pbar = (1.0e-5)*p;

    // Calculate water density, kg/m^3
    double pw = (ko + ca*pbar + cb*pow(pbar, 2.0));
    pw /= (ko + ca*pbar + cb*pow(pbar, 2.0) - pbar);
    pw *= (1.0e3)*po;

    return pw;
}

double EVAP::psychro(double tc, double p){
    /* ***********************************************************************
    Name:     EVAP.psychro
    Input:    - double, air temperature (tc), degrees C
              - double, atm. pressure (p), Pa
    Output:   double, psychrometric constant, Pa/K
    Features: Calculates the psychrometric constant for a given temperature
              and pressure
    Depends:  Global constants:
              - kMa
              - kMv
    Refs:     Allen et al. (1998); Tsilingiris (2008)
    *********************************************************************** */
    // Calculate the specific heat capacity of water, J/kg/K
    double cp = specific_heat(tc);

    // Calculate latent heat of vaporization, J/kg
    double lv = enthalpy_vap(tc);

    // Calculate psychrometric constant, Pa/K
    //   Eq. 8, Allen et al. (1998)
    double ps = (Global::Ma*cp*p)/(Global::Mv*lv);

    return ps;
}

double EVAP::specific_heat(double tc){
    /* ***********************************************************************
    Name:     EVAP.specific_heat
    Input:    double, air temperature (tc), degrees C
    Output:   double, specific heat of moist air, J/kg/K
    Features: Calculates the specific heat of moist air for a given temperature
    Refs:     Eq. 47, Tsilingiris (2008); valid only for air temp 0--100 deg C
    *********************************************************************** */
    double cp;
    if (tc < 0) {
        cp = 1004.5714270;
    } else if (tc > 100) {
        cp = 2031.2260590;
    } else {
        cp = 1.0045714270;
        cp += (2.050632750e-3)*tc;
        cp += -(1.631537093e-4)*tc*tc;
        cp += (6.212300300e-6)*tc*tc*tc;
        cp += -(8.830478888e-8)*tc*tc*tc*tc;
        cp += (5.071307038e-10)*tc*tc*tc*tc*tc;
        cp *= (1.0e3);
    }

    return cp;
}

etr EVAP::get_vals(){
    /* ***********************************************************************
    Name:     EVAP.get_vals
    Input:    None
    Output:   etr struct, daily ET and radiation values
    Features: Returns the current set of daily ET and radiation values
    *********************************************************************** */
    d_etr.cond = cn;
    d_etr.eet = eet_d;
    d_etr.pet = pet_d;
    d_etr.aet = aet_d;
    return d_etr;
}

void EVAP::display(){
    /* ***********************************************************************
    Name:     EVAP.display
    Input:    None
    Output:   None
    Features: Prints the current set of daily STASH variables
    *********************************************************************** */
    solar.display();

    printf("EVAP variable list:\n");
    printf("  s: %0.6f Pa/K\n", s);
    printf("  lv: %0.6f MJ/kg\n", (1.0e-6)*lv);
    printf("  Patm: %0.6f bar\n", (1.0e-5)*patm);
    printf("  pw: %0.6f kg/m^3\n", pw);
    printf("  gamma: %0.6f Pa/K\n", g);
    printf("  Econ: %0.6f mm^3/J\n", (1.0e9)*econ);
    printf("  Cn: %0.6f mm\n", cn);
    printf("  rx: %0.6f\n", rx);
    printf("  hi: %0.6f degrees\n", hi);
    printf("  EET: %0.6f mm\n", eet_d);
    printf("  PET: %0.6f mm\n", pet_d);
    printf("  AET: %0.6f mm\n", aet_d);
}
