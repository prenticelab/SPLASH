#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * DATA.cpp
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
 * This script reads a CSV file with input for the STASH model.
 * 
 * ----------
 * changelog:
 * ----------
 * 
 * 
 * //////////////////////////////////////////////////////////////////////// */

class DATA {
    private:
        // Variables:
        string file_name;
        int num_lines;                  // number of lines read from file
        vector<double> sf_vec;          // sun hours fraction
        vector<double> tair_vec;        // air temperature, deg C
        vector<double> pn_vec;          // precipitation, mm
        
        // Functions:
        int count_lines(string fname);
        bool inValidDouble(string mystring);
    
    public:
        // Constructors:
        DATA(string fname);
        
        // Functions:
        vector<double> get_all_sf();
        vector<double> get_all_tair();
        vector<double> get_all_pn();
        double get_one_sf(int n);
        double get_one_tair(int n);
        double get_one_pn(int n);
        int nlines();
};

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Constructors: 
// ////////////////////////////////////////////////////////////////////////
DATA::DATA(string fname)
    : file_name(fname)
{
    string line;
    string sf_str, tair_str, pn_str;
    double sf_val, tair_val, pn_val;
    size_t pos;
    vector<int> posit;
    int commas;
    
    // Open stream to file:
    ifstream my_file( fname.c_str() );
    if (my_file.is_open()){
        int i = 0;
        while (!my_file.eof()){
            getline(my_file, line);
            
            // Count the number of commas in the line:
            commas = 0;
            pos = line.find_first_of(",");
            if(pos!=string::npos){   
                commas++;  
            }
            while (pos != string::npos){
                pos = line.find_first_of(",",pos+1);
                if(pos != string::npos){ 
                    commas++; 
                }
            }
            
            if (commas == 2){
                // Save the positions for each comma:
                posit.resize(3, 0);
                pos = line.find_first_of(",");
                int j = 0;
                while (pos != string::npos)
                {
                    posit[j]=pos;
                    pos = line.find_first_of(",",pos+1);
                    j++;
                }
                
                // Skip headerline:
                if (i != 0) {
                    // Extract substrings based on comma locations:
                    sf_str = line.substr(0, posit[0]);
                    tair_str = line.substr((posit[0]+1),(posit[1]-posit[0]-1));
                    pn_str = line.substr((posit[1]+1),string::npos);
                    
                    // Convert strings to doubles:
                    sf_val = atof(sf_str.c_str());
                    tair_val = atof(tair_str.c_str());
                    pn_val = atof(pn_str.c_str());
                    
                    // Append doubles to vectors:
                    sf_vec.push_back(sf_val);
                    tair_vec.push_back(tair_val);
                    pn_vec.push_back(pn_val);
                }
            } // end if commas==2
            
            i++;
        } // end while file != eof
        num_lines = (i - 1);  // subtract one for last iteration
        num_lines--;          // subtract one for the headerline
    } // end if file is open
}

// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Class Function Definitions
// ////////////////////////////////////////////////////////////////////////
bool DATA::inValidDouble(string mystring){
    /* ***********************************************************************
    Name:     DATA.inValidDouble
    Input:    string (mystring)
    Output:   bool
    Features: This function returns 'false' if string is a valid double. The 
              checks include non-digits, minus signs and decimal places.
    Ref:      Based partially from code found at http://programmersheaven.com/
    *********************************************************************** */
    int n = 0;             // number of decimals
    int m = 0;             // number of minus signs
    bool invstr = false;  // invalid string boolean
    
    for (int i=0; i < (mystring.size()); i++){
        // Check to see if character is a digit //
        if (mystring[i] < '0' || mystring[i] > '9'){
            if (mystring[i] != '.'){
                // Check to see if non-digit is a period or minus
                if (mystring[i] != '-'){
                    // Non-digit is not a period, check for minus sign
                    invstr = true;
                } else if (mystring.size() == 1){
                    // Non-digit was a minus sign; is it a number?
                    invstr = true;
                } else {
                    // add 1 to number of minus signs found
                    m += 1;
                }
            } else if (mystring.size() == 1) {
                // Non-digit was a period, is it a number?
                invstr = true;
            } else {
                // add 1 to number of periods found
                n += 1;
            }
        }
        if (n > 1 || m > 1) {
            // There was more than 1 period or minus  present
            invstr = true;
        } else if ( (n == 1) && (m == 1) && (mystring.size() == 2) ) {
            // String consists of only a minus and a period //
            invstr = true;
        }
    }
    return invstr;
}

vector<double> DATA::get_all_sf() {
    /* ***********************************************************************
    Name:     DATA.get_all_sf
    Input:    None
    Output:   vector<double>
    Features: This function returns the vector of sun hour fractions
    *********************************************************************** */
    return sf_vec;
}

vector<double> DATA::get_all_tair() {
    /* ***********************************************************************
    Name:     DATA.get_all_tair
    Input:    None
    Output:   vector<double>
    Features: This function returns the vector of air temperature
    *********************************************************************** */
    return tair_vec;
}

vector<double> DATA::get_all_pn() {
    /* ***********************************************************************
    Name:     DATA.get_all_pn
    Input:    None
    Output:   vector<double>
    Features: This function returns the vector of precipitation
    *********************************************************************** */
    return pn_vec;
}

double DATA::get_one_sf(int n){
    /* ***********************************************************************
    Name:     DATA.get_one_sf
    Input:    int (n)
    Output:   double
    Features: This function returns the sun hour fraction for a given index.
    *********************************************************************** */
    return sf_vec[n];
}

double DATA::get_one_tair(int n){
    /* ***********************************************************************
    Name:     DATA.get_one_tair
    Input:    int (n)
    Output:   double
    Features: This function returns the air temperature for a given index.
    *********************************************************************** */
    return tair_vec[n];
}

double DATA::get_one_pn(int n){
    /* ***********************************************************************
    Name:     DATA.get_one_pn
    Input:    int (n)
    Output:   double
    Features: This function returns the precipitation for a given index.
    *********************************************************************** */
    return pn_vec[n];
}

int DATA::nlines(){
    /* ***********************************************************************
    Name:     DATA.nlines
    Input:    None
    Output:   int
    Features: This function returns the number of lines in the input file.
    *********************************************************************** */
    return num_lines;
}
