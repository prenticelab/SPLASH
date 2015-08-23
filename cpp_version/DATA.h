#include <string>
#include <vector>

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * DATA.h
 * 
 * 2015-02-08 -- created
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
 * This is the header file for the DATA C++ class.
 * 
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. Added iostream to module list [15.02.17]
 * 02. Created read_csv and read_txt class member functions [15.02.17]
 * 03. Added year to class member variables [15.02.17]
 * 04. Created header guard [15.02.19]
 * 05. Moved all necessary includes to header file [15.02.19]
 *     ---> removed from cpp file
 * 
 * //////////////////////////////////////////////////////////////////////// */
#ifndef DATA_H
#define DATA_H
class DATA {
    private:
        // Variables:
        int num_lines;                         // number of lines read from file
        int year;                              // year for data from file
        std::string file_name;
        std::vector<std::string> file_names;   // all file names
        std::vector<double> sf_vec;            // sun hours fraction
        std::vector<double> tair_vec;          // air temperature, deg C
        std::vector<double> pn_vec;            // precipitation, mm
        
    public:
        // Constructors:
        DATA();
        
        // Functions:
        void read_csv(std::string fname, int y=-1);
        void read_txt(std::string fname, std::string var, int y=-1);
        
        std::vector<double> get_all_sf();
        std::vector<double> get_all_tair();
        std::vector<double> get_all_pn();
        double get_one_sf(int n);
        double get_one_tair(int n);
        double get_one_pn(int n);
        int nlines();
        int get_year();
};
#endif
