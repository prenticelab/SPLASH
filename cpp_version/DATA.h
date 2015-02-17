#include <string>
#include <vector>

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * DATA.h
 * 
 * written by Tyler W. Davis
 * Imperial College London
 * 
 * 2015-02-08 -- created
 * 2015-02-08 -- last updated
 * 
 * ------------
 * description:
 * ------------
 * This is the header file for the DATA C++ class.
 * 
 * //////////////////////////////////////////////////////////////////////// */

class DATA {
    private:
        // Variables:
        int num_lines;                  // number of lines read from file
        int year;                       // year for data from file
        string file_name;
        vector<string> file_names;       // all file names
        vector<double> sf_vec;          // sun hours fraction
        vector<double> tair_vec;        // air temperature, deg C
        vector<double> pn_vec;          // precipitation, mm
        
    public:
        // Constructors:
        DATA();
        
        // Functions:
        void read_csv(string fname, int y=-1);
        void read_txt(string fname, string var, int y=-1);
        
        vector<double> get_all_sf();
        vector<double> get_all_tair();
        vector<double> get_all_pn();
        double get_one_sf(int n);
        double get_one_tair(int n);
        double get_one_pn(int n);
        int nlines();
        int get_year();
};
