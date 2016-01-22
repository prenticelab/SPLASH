#include <cmath>
#include <iostream>
#include <vector>

#include "global.h"
#include "SOLAR.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SOLAR.cpp
 *
 * 2016-01-22 -- created
 * 2016-01-22 -- last updated
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
 * This script includes the definitions for the SOLAR class, which calculates
 * the daily quantities of solar radiation for the SPLASH model.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 *
 * //////////////////////////////////////////////////////////////////////// */
/*
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Class Constructors:
/////////////////////////////////////////////////////////////////////
*/
SOLAR::SOLAR()
     : lat(0), elv(0)
{
}

SOLAR::SOLAR(double a, double b)
     : lat(a), elv(b)
{
}

/*
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Class Function Definitions
/////////////////////////////////////////////////////////////////////
*/
void SOLAR::calculate_daily_fluxes(int n, int y, double sf, double tc){
    /* ***********************************************************************
    Name:     SOLAR.calculate_daily_fluxes
    Input:    - int, day of year (n)
              - int, year (y)
              - double, fraction of sunshine hours (sf)
              - double, mean daily air temperature, C (tc)
    Output:
    Features: Calculates the daily solar radiation fluxes
    *********************************************************************** */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 0. Save day of year
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // @TODO: check for erroneous values
    day = n;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate number of days in year (kN), days
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (y == 0){
        kN = 365;
    } else {
        kN = julian_day((y + 1), 1, 1) - julian_day(y, 1, 1);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate heliocentric longitudes (nu and lambda), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vector<double> my_vals = berger_tls(n);
    nu = my_vals[0];
    lam = my_vals[1];

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate distance factor (dr), unitless
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double kee = pow(Global::e, 2.0);
    rho = (1.0 - kee)/(1.0 + dcos(nu)*Global::e);
    dr = pow(rho, -1.0);
    dr = pow(dr, 2.0);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate declination angle (delta), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delta = dsin(lam)*dsin(Global::eps);
    delta = asin(delta);
    delta /= Global::pir;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Calculate variable substitutes (u and v), unitless
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ru = dsin(delta)*dsin(lat);
    rv = dcos(delta)*dcos(lat);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6. Calculate the sunset hour angle (hs), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((ru/rv) >= 1.0) {
        // Polar day (no sunset)
        hs = 180.0;
    } else if ((ru/rv) <= -1.0) {
        // Polar night (no sunrise)
        hs = 0.0;
    } else {
        hs = -1.0*(ru/rv);
        hs = acos(hs);
        hs /= Global::pir;
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ra_d = (86400.0/Global::PI)*dr*Global::Gsc;
    ra_d *= (ru*hs*Global::pir + rv*dsin(hs));

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 8. Calculate transmittivity (tau), unitless
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau_o = (Global::c + sf*Global::d);
    tau = tau_o*(1.0 + (2.67e-5)*elv);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 9. Calculate daily PPFD (ppfd_d), mol/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ppfd_d = (1.0e-6)*Global::fFEC*(1.0 - Global::alb_vis)*tau*ra_d;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 10. Estimate net longwave radiation (rnl), W/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnl = (Global::b + (1.0 - Global::b)*sf)*(Global::A - tc);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 11. Calculate variable substitute (rw), W/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw = (1.0 - Global::alb_sw)*tau*dr*Global::Gsc;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 12. Calculate net radiation cross-over hour angle (hn), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((rnl - rw*ru)/(rw*rv) >= 1.0) {
        // Net radiation negative all day
        hn = 0;
    } else if ((rnl - rw*ru)/(rw*rv) <= -1.0) {
        // Net radiation positive all day
        hn = 180.0;
    } else {
        hn = acos((rnl - rw*ru)/(rw*rv));
        hn /= Global::pir;
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 13. Calculate daytime net radiation (rn_d), J/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rn_d = Global::pir*hn*(rw*ru - rnl) + rw*rv*dsin(hn);
    rn_d *= (86400.0/Global::PI);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 14. Calculate nighttime net radiation (rnn_d), J/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnn_d = rw*ru*(hs - hn)*Global::pir;
    rnn_d += rw*rv*(dsin(hs) - dsin(hn));
    rnn_d += rnl*(Global::PI - 2.0*hs*Global::pir + hn*Global::pir);
    rnn_d *= (86400.0/Global::PI);
}

double SOLAR::dcos(double x){
    /* ***********************************************************************
    Name:     SOLAR.dcos
    Input:    double, angle, degrees (x)
    Output:   double, cos(x*pi/180)
    Features: Calculates the cosine of an angle given in degrees
    *********************************************************************** */
    double v = cos(x*Global::pir);
    return v;
 }

double SOLAR::dsin(double x){
    /* ***********************************************************************
    Name:     SOLAR.dsin
    Input:    double, angle, degrees (x)
    Output:   double, sin(x*pi/180)
    Features: Calculates the sine of an angle given in degrees
    *********************************************************************** */
    double v = sin(x*Global::pir);
    return v;
 }

vector<double> SOLAR::berger_tls(int n){
    /* ***********************************************************************
    Name:     SOLAR.berger_tls
    Input:    int, day of year
    Output:   vector:
              - double true anomaly, degrees
              - double true longitude, degrees
    Features: Returns true anomaly and true longitude for a given day
    Depends:  - ke
              - komega
              - kN
    Ref:      Berger, A. L. (1978), Long term variations of daily insolation
              and quaternary climatic changes, J. Atmos. Sci., 35, 2362-
              2367.
    *********************************************************************** */
    // Variable substitutes:
    vector<double> rvals;
    double xee = pow(Global::e, 2.0);
    double xec = pow(Global::e, 3.0);
    double xse = sqrt(1.0 - xee);

    // Mean longitude for vernal equinox:
    double xlam = (Global::e/2.0 + xec/8.0)*(1.0 + xse)*dsin(Global::omega);
    xlam -= xee/4.0*(0.5 + xse)*dsin(2.0*Global::omega);
    xlam += xec/8.0*(1.0/3.0 + xse)*dsin(3.0*Global::omega);
    xlam *= 2.0;
    xlam /= Global::pir;

    // Mean longitude for day of year:
    double dlamm = xlam + (n - 80.0)*(360.0/kN);

    // Mean anomaly:
    double anm = (dlamm - Global::omega);
    double ranm = anm*Global::pir;

    // True anomaly:
    double ranv = ranm;
    ranv += (2.0*Global::e - xec/4.0)*sin(ranm);
    ranv += 5.0/4.0*xee*sin(2.0*ranm);
    ranv += 13.0/12.0*xec*sin(3.0*ranm);
    double anv = ranv/Global::pir;

    // True longitude:
    double my_tls = (anv + Global::omega);
    if (my_tls < 0) {
        my_tls += 360.0;
    } else if (my_tls > 360) {
        my_tls -= 360.0;
    }

    // True anomaly:
    double my_nu = (my_tls - Global::omega);
    if (my_nu < 0) {
        my_nu += 360.0;
    }

    rvals.push_back(my_nu);
    rvals.push_back(my_tls);
    return rvals;
 }

int SOLAR::julian_day(int y, int m, int i){
    /* ***********************************************************************
    Name:     SOLAR.julian_day
    Input:    - int, year (y)
              - int, month (m)
              - int, day of month (i)
    Output:   int, Julian Ephemeris Day
    Features: Converts Gregorian date (year, month, day) to Julian
              Ephemeris Day
    Ref:      Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," Astronomical
              Algorithms
    *********************************************************************** */
    if (m <= 2.0){
        y -= 1.0;
        m += 12.0;
    }
    int a = int(y/100);
    int b = 2 - a + int(a/4);
    float jd = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5;
    int jde = int(jd);
    return jde;
 }

srad SOLAR::get_vals(){
    /* ***********************************************************************
    Name:     SOLAR.get_vals
    Input:    None
    Output:   srad struct, daily radiation values
    Features: Returns the current set of daily radiation values
    *********************************************************************** */
    dsr.ru = ru;        // variable substitute
    dsr.rv = rv;        // variable substitute
    dsr.rw = rw;        // variable substitute
    dsr.rnl = rnl;      // net longwave radiation, W/m^2
    dsr.hn = hn;        // cross-over hour angle, degrees
    dsr.rn_d = rn_d;    // daytime net radiation, J/m^2
    dsr.rnn_d = rnn_d;  // nighttime net radiation, J/m^2

    return dsr;
}

void SOLAR::display(){
    /* ***********************************************************************
    Name:     SOLAR.display
    Input:    None
    Output:   None
    Features: Prints the current set of daily SOLAR variables
    *********************************************************************** */
    cout<<"SOLAR variable list:"<<endl;
    cout<<"  kN: " << kN << endl;
    cout<<"  nu: " << nu << " degrees" << endl;
    cout<<"  lambda: " << lam << " degrees" << endl;
    cout<<"  rho: " << rho << endl;
    cout<<"  dr: " << dr << endl;
    cout<<"  delta: " << delta << " degrees" << endl;
    cout<<"  ru: " << ru << endl;
    cout<<"  rv: " << rv << endl;
    cout<<"  rw: " << rw << endl;
    cout<<"  hs: " << hs << " degrees" << endl;
    cout<<"  hn: " << hn << " degrees" << endl;
    cout<<"  tau_o: " << tau_o << endl;
    cout<<"  tau: " << tau << endl;
    cout<<"  PPFD: " << ppfd_d << " mol/m^2" << endl;
    cout<<"  Rnl: " << rnl << " W/m^2" << endl;
    cout<<"  Ho: " << (1.0e-6)*ra_d << " MJ/m^2" << endl;
    cout<<"  Hn: " << (1.0e-6)*rn_d << " MJ/m^2" << endl;
    cout<<"  Hnn: " << (1.0e-6)*rnn_d << " MJ/m^2" << endl;
}
