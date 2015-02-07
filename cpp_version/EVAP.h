#include <vector>

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * EVAP.h
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-07 -- created
 * 2015-02-07 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This is the C++ header file for the EVAP class.
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
