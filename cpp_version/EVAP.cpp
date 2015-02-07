#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * EVAP.cpp
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-06 -- created
 * 2015-02-07 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This script includes the definitions for the EVAP class, which calculates 
 * the daily quantities of radiation, evaporation, and condensation for the 
 * STASH model.
 * 
 * ----------
 * changelog:
 * ----------
 * 01. finished class function definitions and debugged [15.02.07]
 * 
 * //////////////////////////////////////////////////////////////////////// */


class EVAP {
    private:
        // Constants:
        double kA;       // constant for Rnl
        double kalb_sw;  // shortwave albedo
        double kalb_vis; // visible light albedo
        double kb;       // constant for Rnl
        double kc;       // cloudy transmittivity
        double kCw;      // supply constant, mm/hr 
        double kd;       // angular coefficient of transmittivity
        double ke;       // eccentricity
        double keps;     // obliquity, degrees 
        double kfFEC;    // from flux to energy conversion, umol/J 
        double kG;       // gravitational acceleration, m/s^2 
        double kGsc;     // solar constant, W/m^2 
        double kL;       // temperature lapse rate, K/m 
        double kMa;      // molecular weight of dry air, kg/mol 
        double kMv;      // molecular weight of water vapor, kg/mol 
        double kPo;      // standard atmosphere, Pa 
        double kR;       // universal gas constant, J/mol/K 
        double kTo;      // base temperature, K 
        double kWm;      // soil moisture capacity, mm 
        double kw;       // entrainment factor
        double komega;   // longitude of perihelion, degrees 
        double kPI;      // pi
        
        // Variables:
        int kN;                    // days in year
        double my_nu;              // true anomaly, degrees
        double my_lambda;          // true longitude, degrees
        double my_rho;             // rel. earth-sun distance
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
        double s;                  // slope sat. vap. press. temp. curve, Pa/K
        double lv;                 // enthalpy of vaporization, J/kg
        double pw;                 // density of water, kg/m^3
        double g;                  // psychrometric constant, Pa/K
        double econ;               // water-to-energy conversion, m^3/J
        double cn;                 // daily condensation, mm
        double eet_d, pet_d, aet_d; // daily ET, mm
        double cos_hi;             // cosine of hour angle, hi
        
        // Functions:
        double dcos(double x);
        double dsin(double x);
        vector<double> berger_tls(int n);
        int julian_day(int y, int m, int i);
        double sat_slope(double tc);
        double enthalpy_vap(double tc);
        double elv2pres(double z);
        double density_h2o(double tc, double p);
        double psychro(double tc, double p);
    
    public:
        // Constructors
        //EVAP(double lat, int n, double elv);  // default constructor
        //EVAP(double lat, int n, double elv, int y);  // default w/ year
        EVAP(double lat, int n, double elv, int y, float sf, float tc, 
             float sw);
        
        // Get Variable Functions:
        double get_ho();
        double get_ppfd();
        double get_aet();
        double get_cn();
        
        // Print variables to screen:
        void display();
};

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Constructors: 
// ////////////////////////////////////////////////////////////////////////
EVAP::EVAP(double lat, int n, double elv, int y, float sf, float tc, float sw){
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 0. Initialize constants:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    kA = 107.0;       // (Monteith & Unsworth, 1990)
    kalb_sw = 0.17;   // (Federer, 1968)
    kalb_vis = 0.03;  // (Sellers, 1985)
    kb = 0.20;        // (Linacre, 1968)
    kc = 0.25;        // (Linacre, 1968)
    kCw = 1.05;       // (Federer, 1982)
    kd = 0.50;        // (Linacre, 1968)
    ke = 0.0167;      // (Berger, 1978)
    keps = 23.44;     // (Berger, 1978)
    kfFEC = 2.04;     // (Meek et al., 1984)
    kG = 9.80665;     // (Allen, 1973)
    kGsc = 1360.8;    // (Kopp & Lean, 2011)
    kL = 0.0065;      // (Allen, 1973)
    kMa = 0.028963;   // (Tsilingiris, 2008)
    kMv = 0.01802;    // (Tsilingiris, 2008)
    kPo = 101325;     // (Allen, 1973)
    kR = 8.3143;      // (Allen, 1973)
    kTo = 298.15;     // (Prentice, unpublished)
    kWm = 150.0;      // (Cramer & Prentice, 1988)
    kw = 0.26;        // (Lhomme, 1997; Priestley & Taylor, 1972)
    komega = 283.0;   // (Berger, 1978)
    kPI = 3.14159265; // pi
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate number of days in year (kN), days
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (y == 0){
        kN = 365;
    } else {
        kN = julian_day((y+1), 1, 1) - julian_day(y, 1, 1);
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate heliocentric longitudes (nu and lambda), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vector<double> my_vals = berger_tls(n);
    my_nu = my_vals[0];
    my_lambda = my_vals[1];
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate distance factor (dr), unitless
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double kee = pow(ke, 2.0);
    my_rho = (1.0 - kee)/(1.0 + ke*dcos(my_nu));
    dr = pow(my_rho, -1.0);
    dr = pow(dr, 2.0);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate declination angle (delta), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double pir = kPI/180.0;
    delta = dsin(my_lambda)*dsin(keps);
    delta = asin(delta);
    delta /= pir;
    
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
        hs /= pir;
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ra_d = (86400.0/kPI)*kGsc*dr*(ru*pir*hs + rv*dsin(hs));
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 8. Calculate transmittivity (tau), unitless
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau_o = (kc + kd*sf);
    tau = tau_o*(1.0 + (2.67e-5)*elv);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 9. Calculate daily PPFD (ppfd_d), mol/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra_d;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 10. Estimate net longwave radiation (rnl), W/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnl = (kb + (1.0 - kb)*sf)*(kA - tc);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 11. Calculate variable substitute (rw), W/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw = (1.0 - kalb_sw)*tau*kGsc*dr;
    
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
        hn /= pir;
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 13. Calculate daytime net radiation (rn_d), J/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rn_d = (86400.0/kPI)*(hn*pir*(rw*ru - rnl) + rw*rv*dsin(hn));
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 14. Calculate nighttime net radiation (rnn_d), J/m^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnn_d = rw*ru*(hs - hn)*pir;
    rnn_d += rw*rv*(dsin(hs) - dsin(hn)); 
    rnn_d += rnl*(kPI - 2.0*hs*pir + hn*pir);
    rnn_d *= (86400.0/kPI);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 15. Calculate water-to-energy conversion (econ), m^3/J
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double patm = elv2pres(elv);
    s = sat_slope(tc);
    lv = enthalpy_vap(tc);
    pw = density_h2o(tc, patm);
    g = psychro(tc, patm);
    econ = s/(lv*pw*(s + g));
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 16. Calculate daily condensation (wc), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cn = (1.0e3)*econ*abs(rnn_d);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 17. Estimate daily EET (eet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eet_d = (1.0e3)*econ*rn_d;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 18. Estimate daily PET (pet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pet_d = (1.0 + kw)*eet_d;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx = (3.6e6)*(1.0 + kw)*econ;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 20. Calculate the intersection hour angle (hi), degrees
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
        hi /= pir;
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 21. Estimate daily AET (aet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aet_d = sw*hi*pir;
    aet_d += rx*rw*rv*(dsin(hn) - dsin(hi)); 
    aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*pir;
    aet_d *= (24.0/kPI);
}

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Function Definitions
// ////////////////////////////////////////////////////////////////////////
double EVAP::dcos(double x){
    /* ***********************************************************************
    Name:     EVAP.dcos
    Input:    double, angle, degrees (x)
    Output:   double, cos(x*pi/180)
    Features: Calculates the cosine of an angle given in degrees
    *********************************************************************** */
    double v = kPI/180.0;
    v *= x;
    v = cos(v);
    return v;
}

double EVAP::dsin(double x){
    /* ***********************************************************************
    Name:     EVAP.dsin
    Input:    double, angle, degrees (x)
    Output:   double, sin(x*pi/180)
    Features: Calculates the sine of an angle given in degrees
    *********************************************************************** */
    double v = kPI/180.0;
    v *= x;
    v = sin(v);
    return v;
}

vector<double> EVAP::berger_tls(int n){
    /* ***********************************************************************
    Name:     EVAP.berger_tls
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
    double xee = pow(ke, 2.0); 
    double xec = pow(ke, 3.0);
    double xse = sqrt(1.0 - xee);
    double pir = (kPI/180.0);
    
    // Mean longitude for vernal equinox:
    double xlam = (ke/2.0 + xec/8.0)*(1.0 + xse)*dsin(komega); 
    xlam -= xee/4.0*(0.5 + xse)*dsin(2.0*komega); 
    xlam += xec/8.0*(1.0/3.0 + xse)*dsin(3.0*komega);
    xlam *= 2.0;
    xlam /= pir;
    
    // Mean longitude for day of year:
    double dlamm = xlam + (n - 80.0)*(360.0/kN);
    
    // Mean anomaly:
    double anm = (dlamm - komega);
    double ranm = anm*pir;
    
    // True anomaly:
    double ranv = ranm;
    ranv += (2.0*ke - xec/4.0)*sin(ranm);
    ranv += 5.0/4.0*xee*sin(2.0*ranm);
    ranv += 13.0/12.0*xec*sin(3.0*ranm);
    double anv = ranv/pir;
    
    // True longitude:
    double my_tls = (anv + komega);
    if (my_tls < 0) {
        my_tls += 360.0;
    } else if (my_tls > 360) {
        my_tls -= 360.0;
    }
     
    // True anomaly:
    double my_nu = (my_tls - komega);
    if (my_nu < 0) {
        my_nu += 360.0;
    }
    
    rvals.push_back(my_nu);
    rvals.push_back(my_tls);
    return rvals;
}

int EVAP::julian_day(int y, int m, int i){
    /* ***********************************************************************
    Name:     EVAP.julian_day
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
              - kPo
              - kTo
              - kL
              - kMa
              - kG
              - kR
    Ref:      Allen et al. (1998)
    *********************************************************************** */
    double ep = (kG*kMa)/(kR*kL);
    double pa = (1.0 - kL*z/kTo);
    pa = pow(pa, ep);
    pa *= kPo;
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
    //   Eq. 47, Tsilingiris (2008)
    double cp = 1.0045714270;
    cp += (2.050632750e-3)*tc;
    cp += -(1.631537093e-4)*tc*tc;
    cp += (6.212300300e-6)*tc*tc*tc;
    cp += -(8.830478888e-8)*tc*tc*tc*tc;
    cp += (5.071307038e-10)*tc*tc*tc*tc*tc;
    cp *= (1.0e3);
    
    // Calculate latent heat of vaporization, J/kg
    double lv = enthalpy_vap(tc);
    
    // Calculate psychrometric constant, Pa/K
    //   Eq. 8, Allen et al. (1998)
    double ps = (cp*kMa*p)/(kMv*lv);
    
    return ps;
}

double EVAP::get_ho(){
    return ra_d;
}

double EVAP::get_ppfd(){
    return ppfd_d;
}

double EVAP::get_aet(){
    return aet_d;
}

double EVAP::get_cn(){
    return cn;
}

void EVAP::display(){
    cout<<"Variable list:"<<endl;
    cout<<"  nu: " << my_nu << " degrees" << endl;
    cout<<"  lambda: " << my_lambda << " degrees" << endl;
    cout<<"  delta: " << delta << " degrees" << endl;
    cout<<"  dr: " << dr << endl;
    cout<<"  kN: " << kN << endl;
    cout<<"  hs: " << hs << " degrees" << endl;
    cout<<"  hn: " << hn << " degrees" << endl;
    cout<<"  hi: " << hi << " degrees" << endl;
    cout<<"  Ho: " << (1.0e-6)*ra_d << " MJ/m^2" << endl;
    cout<<"  HN: " << (1.0e-6)*rn_d << " MJ/m^2" << endl;
    cout<<"  PAR: " << ppfd_d << " mol/m^2" << endl;
    cout<<"  tau: " << tau << endl;
    cout<<"  s: " << s << " Pa/K" << endl;
    cout<<"  g: " << g << " Pa/K" << endl;
    cout<<"  lv: " << (1.0e-6)*lv << " MJ/kg" << endl;
    cout<<"  pw: " << pw << " kg/m^3" << endl;
    cout<<"  Econ: " << econ << " m^3/J" << endl;
    cout<<"  Cn: " << cn << " mm" << endl;
    cout<<"  EET: " << eet_d << " mm" << endl;
    cout<<"  PET: " << pet_d << " mm" << endl;
    cout<<"  AET: " << aet_d << " mm" << endl;
}
