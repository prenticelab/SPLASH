import numpy as np
from ctypes import *

#Python wrapper class for the SPLASH model written in fortran
#   Loads the fortran shared library and sets some of its function input and outputs
#   Gives easy wrapper names for several of the splash functions
#   Note that the shared compiler puts '__' at the begining of every function, this is a reserved keyword for python classes and as such we use the 'get_function' get around this
class SPLASH_fortran:
    def __init__(self):
        #Load in the library
        self.nextgen_splash = CDLL('../f90_version/libsplash.so')

        #Set how the inputs and outputs for the julian day function work
        self.get_function('__splash_MOD_get_julian_day').restype = c_int 
        self.get_function('__splash_MOD_get_julian_day').argtypes = [ 
                POINTER(c_int),
                POINTER(c_int),
                POINTER(c_int)
                ]

        #Create a type yearlong, note we choose the maximum number of days in a year for this
        self.yearlong = c_double*366

        #Set how the inputs and outputs for the spin up function work
        self.get_function('__splash_MOD_spin_up_sm').restype = c_double 
        self.get_function('__splash_MOD_spin_up_sm').argtypes = [ 
                POINTER(c_int),
                POINTER(c_double),
                POINTER(c_double),
                POINTER(self.yearlong),
                POINTER(self.yearlong),
                POINTER(self.yearlong),
                POINTER(c_int),
                POINTER(c_double),
                POINTER(c_double),
                POINTER(c_double) 
                ]

    #Since python mangles any functions that start with '__' we create a different way of calling them to stop the mangling
    def get_function(self, name):
        return getattr(self.nextgen_splash, name)

    def get_julian_day(self, year, month, day):
        return self.get_function('__splash_MOD_get_julian_day')(
                byref(c_int(year)),
                byref(c_int(month)),
                byref(c_int(day)),
                )

    def write_output():
        self.get_function('__splash_MOD_write_to_file')()

    def spin_up(self, lat, elv, precip, temp, sunshine_f, year, eccentricity=0.0167, obliquity=23.44, long_of_peri=283.0):
        #TODO: Find out how to set initial soil moisture values
        #waterbal%sm = 0.0
        #waterbal%ro = 0.0

        #Find how many days in this particular year
        days_in_year = self.get_julian_day( year+1, 1, 1 ) - self.get_julian_day( year, 1, 1 )
    
        #Finally spin up the splash model with given parameters and modern focings
        getattr(self.nextgen_splash, '__splash_MOD_spin_up_sm')(
                byref(c_int(year)),
                byref(c_double(lat)),
                byref(c_double(elv)), 
                byref(self.yearlong(*precip)), 
                byref(self.yearlong(*temp)),
                byref(self.yearlong(*sunshine_f)),
                byref(c_int(days_in_year)), 
                byref(c_double(eccentricity)),  #Eccentricity
                byref(c_double(obliquity)),   #Obliquity 
                byref(c_double(long_of_peri))    #Longitude of perihelion
                )

    #Run throught all the daily output variables and return them
    def get_result(self):
        daily_output_variable_names = {
                "daet"  : "__splash_MOD_outdaet",
                "dcn"   : "__splash_MOD_outdcn",
                "deet"  : "__splash_MOD_outdeet",
                "dpet"  : "__splash_MOD_outdpet",
                "dppfd" : "__splash_MOD_outdppfd",
                "dra"   : "__splash_MOD_outdra",
                "drn"   : "__splash_MOD_outdrn",
                "dro"   : "__splash_MOD_outdro",
                "dsm"   : "__splash_MOD_outdsm"}
        return { key : np.array(self.yearlong.in_dll(self.nextgen_splash, daily_output_variable_names[key])) for key in daily_output_variable_names }

#Only run if this is called from the command line
if __name__ == "__main__":
    #Import the python data module
    import sys
    sys.path.insert(0, '../py_version')
    from data import DATA

    #Load in the input data 
    my_data = DATA()
    my_file = '../../data/example_data.csv'
    my_data.read_csv(my_file)

    # Consistency Test #4: Spin-Up
    my_lat = 37.7
    my_elv = 142.

    sp = SPLASH_fortran() #Initalise a model
    sp.spin_up(my_lat, my_elv, my_data.pn_vec, my_data.tair_vec, my_data.sf_vec, my_data.year) #Spin up the model
    print(sp.get_result()['dsm']) #Read a sample of the span up model
