#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class EVAP {
    private:
        // Constants:
        const double kA = 107.0     // (Monteith & Unsworth, 1990)
        const double kalb_sw = 0.17 // (Federer, 1968)
        const double kalb_vis = 0.03 // (Sellers, 1985)
        const double kb = 0.20      // (Linacre, 1968)
        const double kc = 0.25      // (Linacre, 1968)
        const double kCw = 1.05     // mm/hr (Federer, 1982)
        const double kd = 0.50      // (Linacre, 1968)
        const double ke = 0.0167    // (Berger, 1978)
        const double keps = 23.44   // degrees (Berger, 1978)
        const double kfFEC = 2.04   // umol/J (Meek et al., 1984)
        const double kG = 9.80665   // m/s^2 (Allen, 1973)
        const double kGsc = 1360.8  // W/m^2 (Kopp & Lean, 2011)
        const double kL = 0.0065    // K/m (Allen, 1973)
        const double kMa = 0.028963 // kg/mol (Tsilingiris, 2008)
        const double kMv = 0.01802  // kg/mol (Tsilingiris, 2008)
        const double kPo = 101325   // Pa (Allen, 1973)
        const double kR = 8.3143    // J/mol/K (Allen, 1973)
        const double kTo = 298.15   // K (Prentice, unpublished)
        const double kWm = 150.0    // mm (Cramer & Prentice, 1988)
        const double kw = 0.26      // (Lhomme, 1997; Priestley & Taylor, 1972)
        const double komega = 283.0 // degrees (Berger, 1978)
        const double kPI = 3.14159265 // pi
        
        // Variables:
        double kN;                 // days in year
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
        
        // Get Variable Functions
        double get_ho();
        double get_ppfd();
        double get_aet();
        double get_cn();
};

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Constructors: 
// ////////////////////////////////////////////////////////////////////////
EVAP::EVAP(double lat, int n, double elv, int y, float sf, float tc, float sw){
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate number of days in year (kN), days
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (y == 0){
        kN = 365
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
    double kee = pow(ke, 2);
    my_rho = (1.0 - kee)/(1.0 + ke*dcos(my_nu));
    
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
    rw = (1.0 - kalb_sw)*tau*kGsc*dr
    
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
    rnn_d += rnl*(kPI - 2.0*hs*pir + hn*pir));
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
    rx = (3.6e6)*(1.0 + kw)*econ
    
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
    return cos(v);
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
    return sin(v);
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
    Ref:      Berger, A. L. (1978), Long term variations of daily insolation
              and quaternary climatic changes, J. Atmos. Sci., 35, 2362-
              2367.
    *********************************************************************** */
    // Variable substitutes:
    vector<double> rvals;
    double xee = pow(ke, 2.0); 
    double exec = pow(ke, 3.0);
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
    ranv += 13.0/12.0*xec*sin(3.0*ranm));
    double anv = ranv/pir;
    
    // True longitude:
    double my_tls = (anv + komega);
    if (my_tls < 0) {
        my_tls += 360.0
    } else if (my_tls > 360) {
        my_tls -= 360.0
    }
     
    // True anomaly:
    double my_nu = (my_tls - komega);
    if (my_nu < 0) {
        my_nu += 360.0
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
    int a = int(y/100)
    int b = 2 - a + int(a/4)
    int jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
    return jde
}
