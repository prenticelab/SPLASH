/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * STASH.h
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-17 -- created
 * 2015-02-17 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This is the header file for the C++ STASH class.
 * 
 * //////////////////////////////////////////////////////////////////////// */
class STASH {
    private:
        // Constants:
        double kCw;      // supply constant, mm/hr 
        double kWm;      // soil moisture capacity, mm 
        
        // Static variables:
        double lat;                       // latitude, degrees
        double elv;                       // elevation, meters
        
        // Daily status variables:
        double ho;                        // daily solar irradiation, J/m2
        double hn;                        // daily net radiation, J/m2
        double ppfd;                      // daily PPFD, mol/m2
        double cond;                      // daily condensation water, mm
        double wn;                        // daily soil moisture, mm
        double precip;                    // daily precipitation, mm
        double ro;                        // daily runoff, mm
        double eet;                       // daily equilibrium ET, mm
        double pet;                       // daily potential ET, mm
        double aet;                       // daily actual ET, mm
        
    public:
        // Constructors:
        STASH(double latitude, double elevation);
        
        // Functions:
        void run_one_day(int n, int y, double wn, double sf, double tc, 
                         double pn);
        double get_lat();
        double get_elv();
};
