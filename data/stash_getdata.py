#!/usr/bin/python
#
# stash_getdata.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2015-01-22 -- created
# 2015-01-26 -- last updated
#
# ------------
# description:
# ------------
# The script creates a CSV of the necessary daily meteorological variables for
# running the STASH code, including: fraction of bright sunshine, Sf (unitless), 
# air temperature, Tair (deg. C), and precipitation, Pn (mm).
#
# Supported input sources currently include:
#   Sf:
#    * CRU TS3.2 Cld (monthly cloudiness fraction)
#   Tair:
#    * CRU TS3.2 Tmn (monthly mean daily air temperature)
#    * WATCH daily Tair (daily mean air temperature)
#   Pn
#    * CRU TS3.2 Pre (monthly total precipitation)
#    * WATCH daily Rainf (daily mean mass flow rate of rainfall)
#
# ----------
# changelog:
# ----------
# 01. created STASH_DATA class [15.01.26]
#
###############################################################################
## IMPORT MODULES
###############################################################################
import datetime
import glob
import numpy
import os.path
from scipy.io import netcdf

###############################################################################
## CLASSES
###############################################################################
class STASH_DATA:
    """
    Name:     STASH_DATA
    Features: Processes daily data for STASH model
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, user_lat, user_lon, user_elv = 0):
        """
        Name:     STASH_DATA.__init__
        Input:    - float, latitude, degrees north (user_lat)
                  - float, longitude, degrees east (user_lon)
                  - float, elevation, meters AMSL (user_elv)
        Features: Initialize class
        """
        # Initialize variable 'is set' booleans:
        self.date_is_set = False         # for timestamp
        self.cru_sf_is_set = False       # for CRU-based Sf
        self.cru_pn_is_set = False       # for CRU-based Pn
        self.cru_tair_is_set = False     # for CRU-based Tair
        #
        # Initialize data directory and source preference dictionaries:
        self.data_dir = {'cld' : '',
                         'pre' : '',
                         'tmn' : '',
                         'Rainf' : '',
                         'Tair' : ''}
        #
        self.data_src = {'sf' : '',
                         'tair' : '',
                         'pn' : ''}
        #
        # Check user inputs:
        self.elv = user_elv 
        self.elv2pres()
        if user_lat > 90.0 or user_lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.lat = user_lat
        if user_lon > 180.0 or user_lon < -180.0:
            print "Longitude outside range of validity (-180 to 180)!"
            exit(1)
        else:
            self.lon = user_lon
        #
        # Find the 0.5 degree pixel associated with user coordinates: 
        (self.px_lon, self.px_lat) = self.grid_centroid()
        (self.px_x, self.px_y) = self.get_xy_nc()
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def set_cru_cld_dir(self, d):
        """
        Name:     STASH_DATA.set_cru_cld_dir
        Input:    str, directory path (d)
        Features: Define the directory to CRU monthly cloudiness netCDF file
        """
        self.data_dir['cld'] = d
    #
    def set_cru_pre_dir(self, d):
        """
        Name:     STASH_DATA.set_pre_cld_dir
        Input:    str, directory path (d)
        Features: Define the directory to CRU monthly precipitation netCDF file
        """
        self.data_dir['pre'] = d
    #
    def set_cru_tmn_dir(self, d):
        """
        Name:     STASH_DATA.set_tmn_cld_dir
        Input:    str, directory path (d)
        Features: Define the directory to CRU monthly mean air temperature  
                  netCDF file
        """
        self.data_dir['tmn'] = d
    #
    def set_watch_rainf_dir(self, d):
        """
        Name:     STASH_DATA.set_watch_rainf_dir
        Input:    str, directory path (d)
        Features: Define the directory to WATCH daily rainfall netCDF file
        """
        self.data_dir['Rainf'] = d
    #
    def set_watch_tair_dir(self, d):
        """
        Name:     STASH_DATA.set_watch_tair_dir
        Input:    str, directory path (d)
        Features: Define the directory to WATCH daily air temperature netCDF 
                  file
        """
        self.data_dir['Tair'] = d
    #
    def set_output_dir(self, d):
        """
        Name:     STASH_DATA.set_output_dir
        Input:    str, directory path (d)
        Features: Define the directory where the output data will be written
                  and create the output file
        Depends:  writeout
        """
        self.output_file = d + 'stash_data_out.csv'
        header = 'date,sf,tair,pn\n'
        if os.path.isfile(self.output_file):
            user_resp = raw_input("File already exists, overwrite (y/n)? ")
            if user_resp.lower() == 'y':
                self.writeout(self.output_file, header)
        else:
            self.writeout(self.output_file, header)
    #
    def set_sf_source(self, n):
        """
        Name:     STASH_DATA.set_sf_source
        Input:    str, source preference (n)
        Features: Define the source preference for fractional sunshine hours
                  (e.g., 'cru' or 'watch')
        """
        self.data_src['sf'] = n
    #
    def set_pn_source(self, n):
        """
        Name:     STASH_DATA.set_pn_source
        Input:    str, source preference (n)
        Features: Define the source preference for precipitation (e.g., 'cru' 
                  or 'watch')
        """
        self.data_src['pn'] = n
    #
    def set_tair_source(self, n):
        """
        Name:     STASH_DATA.set_tair_source
        Input:    str, source preference (n)
        Features: Define the source preference for air temperature (e.g., 'cru' 
                  or 'watch')
        """
        self.data_src['tair'] = n
    #
    def add_one_month(self, dt0):
        """
        Name:     STASH_DATA.add_one_month
        Input:    datetime date
        Output:   datetime date
        Features: Adds one month to datetime
        Ref:      A. Balogh (2010), ActiveState Code
                  http://code.activestate.com/recipes/577274-subtract-or-add-a-
                  month-to-a-datetimedate-or-datet/
        """
        dt1 = dt0.replace(day=1)
        dt2 = dt1 + datetime.timedelta(days=32) 
        dt3 = dt2.replace(day=1)
        return dt3
    #
    def density_h2o(self, tc, p):
        """
        Name:     STASH_DATA.density_h2o
        Input:    - float, air temperature (tc), degrees C
                  - float, atmospheric pressure (p), Pa
        Output:   float, density of water, kg/m^3
        Features: Calculates density of water at a given temperature and 
                  pressure
        Ref:      Chen et al. (1977)
        """
        # Calculate density at 1 atm:
        po = (
            0.99983952 + 
            (6.788260e-5)*tc + 
            -(9.08659e-6)*tc*tc +
            (1.022130e-7)*tc*tc*tc + 
            -(1.35439e-9)*tc*tc*tc*tc +
            (1.471150e-11)*tc*tc*tc*tc*tc +
            -(1.11663e-13)*tc*tc*tc*tc*tc*tc + 
            (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc + 
            -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
        )
        #
        # Calculate bulk modulus at 1 atm:
        ko = (
            19652.17 +
            148.1830*tc + 
            -2.29995*tc*tc + 
            0.01281*tc*tc*tc + 
            -(4.91564e-5)*tc*tc*tc*tc + 
            (1.035530e-7)*tc*tc*tc*tc*tc
        )
        #
        # Calculate temperature dependent coefficients:
        ca = (
            3.26138 + 
            (5.223e-4)*tc + 
            (1.324e-4)*tc*tc + 
            -(7.655e-7)*tc*tc*tc + 
            (8.584e-10)*tc*tc*tc*tc
        )
        cb = (
            (7.2061e-5) +
            -(5.8948e-6)*tc + 
            (8.69900e-8)*tc*tc + 
            -(1.0100e-9)*tc*tc*tc + 
            (4.3220e-12)*tc*tc*tc*tc
        )
        #
        # Convert atmospheric pressure to bar (1 bar = 100000 Pa)
        pbar = (1.0e-5)*p
        #
        pw = (1e3)*po*(ko + ca*pbar + cb*pbar**2)
        pw /= (ko + ca*pbar + cb*pbar**2 - pbar)
        return pw
    #
    def elv2pres(self):
        """
        Name:     STASH_DATA.elv2pres
        Input:    None.
        Features: Sets the atmospheric pressure (pascals) based on class elv
        Ref:      Allen et al. (1998)
        """
        # Define constants:
        kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
        kTo = 298.15   # base temperature, K (Prentice, unpublished)
        kL = 0.0065    # temperature lapse rate, K/m (Cavcar, 2000)
        kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
        kG = 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
        kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
        #
        p = kPo*(1.0 - kL*self.elv/kTo)**(kG*kMa/(kR*kL))
        self.patm = p
    #
    def grid_centroid(self):
        """
        Name:     STASH_DATA.grid_centroid
        Input:    None.
        Output:   tuple, longitude latitude pair (my_centroid)
        Features: Returns the nearest 0.5 deg. grid centroid for the class's
                  coordinates based on the Euclidean distance to each of the 
                  four surrounding grids; if any distances are equivalent, the 
                  pixel north and east is selected by default
        """
        # Create lists of regular latitude and longitude:
        grid_res = 0.5
        my_lon = self.lon 
        my_lat = self.lat
        lat_min = -90 + 0.5*grid_res
        lon_min = -180 + 0.5*grid_res
        lat_dim = 360
        lon_dim = 720
        lats = [lat_min + y * grid_res for y in xrange(lat_dim)]
        lons = [lon_min + x * grid_res for x in xrange(lon_dim)]
        #
        # Find bounding longitude:
        centroid_lon = None
        if my_lon in lons:
            centroid_lon = my_lon
        else:
            lons.append(my_lon)
            lons.sort()
            lon_index = lons.index(my_lon)
            bb_lon_min = lons[lon_index-1]
            try:
                bb_lon_max = lons[lon_index+1]
            except IndexError:
                bb_lon_max = lons[-1] + grid_res
            #
        # Find bounding latitude:
        centroid_lat = None
        if my_lat in lats:
            centroid_lat = my_lat
        else:
            lats.append(my_lat)
            lats.sort()
            lat_index = lats.index(my_lat)
            bb_lat_min = lats[lat_index-1]
            try:
                bb_lat_max = lats[lat_index+1]
            except IndexError:
                bb_lat_max = lats[-1] + grid_res
            #
        # Determine nearest centroid:
        # NOTE: if dist_A equals dist_B, then centroid defaults positively 
        #       i.e., north / east
        if centroid_lon and centroid_lat:
            my_centroid = (centroid_lon, centroid_lat)
        elif centroid_lon and not centroid_lat:
            # Calculate the distances between lat and bounding box:
            dist_A = bb_lat_max - my_lat
            dist_B = my_lat - bb_lat_min
            if dist_A > dist_B:
                centroid_lat = bb_lat_min
            else:
                centroid_lat = bb_lat_max
            my_centroid = (centroid_lon, centroid_lat)
        elif centroid_lat and not centroid_lon:
            # Calculate the distances between lon and bounding box:
            dist_A = bb_lon_max - my_lon
            dist_B = my_lon - bb_lon_min
            if dist_A > dist_B:
                centroid_lon = bb_lon_min
            else:
                centroid_lon = bb_lon_max
            my_centroid = (centroid_lon, centroid_lat)
        else:
            # Calculate distances between lat:lon and bounding box:
            # NOTE: if all distances are equal, defaults to NE grid
            dist_A = numpy.sqrt(
                (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0
                )
            dist_B = numpy.sqrt(
                (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0
                )
            dist_C = numpy.sqrt(
                (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0
                )
            dist_D = numpy.sqrt(
                (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0
                )
            min_dist = min([dist_A, dist_B, dist_C, dist_D])
            #
            # Determine centroid based on min distance:
            if dist_A == min_dist:
                my_centroid = (bb_lon_max, bb_lat_max)
            elif dist_B == min_dist:
                my_centroid = (bb_lon_max, bb_lat_min)
            elif dist_C == min_dist:
                my_centroid = (bb_lon_min, bb_lat_max)
            elif dist_D == min_dist:
                my_centroid = (bb_lon_min, bb_lat_min)
                #
        # Return nearest centroid:
        return my_centroid
    #
    def get_daily_status(self, d, write_out=True):
        """
        Name:     STASH_DATA.get_daily_status
        Input:    datetime.date, current date (d)
        Features: Process the status of daily variables
        Depends:  - process_sf
                  - process_tair
                  - process_pn
                  - get_month_days
                  - save_to_file
        """
        if self.date_is_set:
            # Check to see if we're in the same month:
            if self.timestamp.replace(day=1) == d.replace(day=1):
                self.is_same_month = True
            else:
                self.is_same_month = False
                self.nm = self.get_month_days(d)
            #
            # Update timestamp:
            self.timestamp = d
            #
            # Process data:
            sf = self.process_sf(d)
            tair = self.process_tair(d)
            pn = self.process_pn(d, tair)
            #
        else:
            self.timestamp = d
            self.nm = self.get_month_days(d)
            self.date_is_set = True
            self.is_same_month = False
            #
            # Process data:
            sf = self.process_sf(d)
            tair = self.process_tair(d)
            pn = self.process_pn(d, tair)
        #
        # Save variables to class:
        self.sf = sf
        self.tair = tair
        self.pn = pn
        #
        # Write to file:
        if write_out:
            self.save_to_file(d, sf, tair, pn)
    #
    def get_daily_watch(self, v, ct):
        """
        Name:     STASH_DATA.get_daily_watch
        Input:    - str, variable of interest (v)
                  - datetime.date, current time (ct)
        Output:   numpy nd.array
        Features: Returns 360x720 monthly WATCH dataset for a given variable 
                  of interest (e.g., Tair, Rainf)
        """
        # Save class variable to local variable:
        d = self.data_dir[v]
        #
        # Search directory for netCDF file:
        my_path = '%s%s*%d%02d.nc' % (d, v, ct.year, ct.month)
        try:
            my_file = glob.glob(my_path)[0]
        except IndexError:
            print "No WATCH file was found for variable: ", v
            print "and month: ", ct.month 
        else:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(my_file, "r")
            #
            #   VARIABLES:
            #    * day (int), 
            #    * lat (grid box center, degrees north)
            #    * lon (grid box center, degrees east)
            #    * timestp (days since beginning of month)
            #
            #   DATA:
            #    * 'Rainf' units = kg m^-2 s^-1
            #    * 'Tair' units = Kelvin  
            #    *  Missing value = 1.00e+20
            #
            # Find time index:
            f_time = f.variables['day'].data
            ti = numpy.where(ct.day == f_time)[0][0]
            #
            # Get the spatial data for current time:
            f_data = f.variables[v].data[ti]
            #
            f.close()
            return f_data
    #
    def get_month_days(self, ts):
        """
        Name:     STASH_DATA.get_month_days
        Input:    datetime.date (ts)
        Output:   int
        Depends:  add_one_month
        Features: Returns the number of days in the month
        """
        ts1 = ts.replace(day=1)
        ts2 = self.add_one_month(ts1)
        dts = (ts2 - ts1).days
        return dts
    #
    def get_monthly_cru(self, v, ct):
        """
        Name:     STASH_DATA.get_monthly_cru
        Input:    - str, variable of interest (v)
                  - datetime.date, current time (ct)
        Output:   numpy nd.array
        Depends:  get_time_index
        Features: Returns 360x720 monthly CRU TS dataset for a given variable 
                  of interest (e.g., cld, pre, tmp)
        """
        # Save class variables to local variables
        d = self.data_dir[v]
        #
        # Search directory for netCDF file:
        try:
            my_file = glob.glob(d + "*" + v + ".dat.nc")[0]
        except IndexError:
            print "No CRU file was found for variable: ", v
        else:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(my_file, "r")
            #
            # Save data for variables of interest:
            # NOTE: for CRU TS 3.21: 
            #       variables: 'lat', 'lon', 'time', v
            #       where v is 'tmp', 'pre', 'cld', 'vap'
            # LAT:  -89.75 -- 89.75
            # LON:  -179.75 -- 179.75
            # TIME:
            #       units: days since 1900-1-1
            #       shape: (1344,)
            #       values: mid-day of each month (e.g., 15th or 16th)
            # DATA:
            #       'cld' units = %
            #       'pre' units = mm
            #       'tmp' units = deg. C
            #       Missing value = 9.96e+36
            # Save the base time stamp:
            bt = datetime.date(1900,1,1)
            #
            # Read the time data as array:
            f_time = f.variables['time'].data
            #
            # Find the time index for the current date:
            ti = self.get_time_index(bt, ct, f_time)
            #
            # Get the spatial data for current time:
            f_data = f.variables[v].data[ti]
            f.close()
            return f_data
    #
    def get_time_index(self, bt, ct, aot):
        """
        Name:     STASH_DATA.get_time_index
        Input:    - datetime date, base timestamp (bt)
                  - datetime date, current timestamp to be found (ct)
                  - numpy nd.array, days since base timestamp (aot)
        Output:   int
        Features: Finds the index in an array of CRU TS days for a given 
                  timestamp 
        """
        # For CRU TS 3.21, the aot is indexed for mid-month days, e.g. 15--16th
        # therefore, to make certain that ct index preceeds the index for the
        # correct month in aot, make the day of the current month less than
        # the 15th or 16th (i.e., replace day with '1'):
        ct = ct.replace(day=1)
        #
        # Calculate the time difference between ct and bt:
        dt = (ct - bt).days
        #
        # Append dt to the aot array:
        aot = numpy.append(aot, [dt,])
        #
        # Find the first index of dt in the sorted array:
        idx = numpy.where(numpy.sort(aot)==dt)[0][0]
        return idx
    #
    def get_xy_nc(self):
        """
        Name:     STASH_DATA.get_xy_nc
        Input:    None.
        Output:   tuple, x-y indices
        Features: Returns array indices for the class's pixel location at 0.5 
                  degree resolution for a netCDF file (i.e., inverted raster)
        """
        lon = self.px_lon
        lat = self.px_lat
        r = 0.5
        # Solve x and y indices:
        x = (lon + 180.0)/r - 0.5
        y = (lat + 90.0)/r - 0.5
        #
        return (int(x), int(y))
    #
    def process_pn(self, ct, tc):
        """
        Name:     STASH_DATA.process_pn
        Input:    - datetime.date, current time (ct)
                  - float, air temperature, deg C (tc)
        Output:   float, precipitation, mm/day
        Depends:  - get_monthly_cru
                  - density_h2o
                  - get_daily_watch
        Features: Processes the daily precipitation
        """
        if self.data_src['pn'] == 'cru':
            # Check to see if you already have this month's data processed:
            if not self.is_same_month or not self.cru_pn_is_set:
                # Either in new month or no CRU-based Pn data is set:
                voi = 'pre'
                d = self.get_monthly_cru(voi, ct)
                pn = d[self.px_y, self.px_x]   # mm/month
                pn /= self.nm                  # mm/day
                self.cru_pn_is_set = True
            else:
                pn = self.pn
        elif self.data_src['pn'] == 'watch':
            # Calculate water density, kg/m^3:
            pw = self.density_h2o(tc, self.patm)
            #
            # Look up WATCH rainfall data:
            voi = 'Rainf'
            d = self.get_daily_watch(voi, ct)
            pn = d[self.px_y, self.px_x]       # kg m^-2 s^-1
            pn /= pw                           # m/s
            pn *= 8.64e7                       # mm/day
        else:
            print "Method currently unavailable for Pn!"
            pn = numpy.nan
        #
        return pn
    #
    def process_sf(self, ct):
        """
        Name:     STASH_DATA.process_sf
        Input:    datetime.date, current time (ct)
        Output:   float, sunshine fraction, unitless
        Depends:  get_monthly_cru
        Features: Processes the daily fractional sunshine hours and save as 
                  class variable
        """
        if self.data_src['sf'] == 'cru':
            # Check to see if you already have this month's data processed:
            if not self.is_same_month or not self.cru_sf_is_set:
                # Either in new month or no CRU-based Sf data is set:
                voi = 'cld'
                d = self.get_monthly_cru(voi, ct)
                sf = d[self.px_y, self.px_x]  # %
                sf /= 100.0                   # unitless
                sf = 1.0 - sf                 # complement of cloudiness
                self.cru_sf_is_set = True
            else:
                sf = self.sf
        else:
            print "Method currently unavailable for Sf!"
            sf = numpy.nan
        #
        return sf
    #
    def process_tair(self, ct):
        """
        Name:     STASH_DATA.process_tair
        Input:    datetime.date, current time (ct)
        Output:   float, air temperature, deg C
        Depends:  - get_monthly_cru
                  - get_daily_watch
        Features: Processes the daily air temperature and save as class variable
        """
        #
        if self.data_src['tair'] == 'cru':
            # Check to see if you already have this month's data processed:
            if not self.is_same_month or not self.cru_tair_is_set:
                # Either in new month or no CRU-based Tair data is set:
                voi = 'tmn'
                d = self.get_monthly_cru(voi, ct)
                tair = d[self.px_y, self.px_x]
                self.cru_tair_is_set = True
            else:
                tair = self.tair
        elif self.data_src['tair'] == 'watch':
            # Look up WATCH air temperature data:
            voi = 'Tair'
            d = self.get_daily_watch(voi, ct)
            tair = d[self.px_y, self.px_x]  # Kelvin
            tair -= 273.15                  # deg. C
        else:
            print "Method currently unavailable for Tair!"
            tair = numpy.nan
        #
        return tair
    #
    def save_to_file(self, ts, sf, tair, pn):
        """
        Name:     STASH_DATA.save_to_file
        Input:    - datetime.date, timestamp (ts)
                  - float, sunshine fraction (sf)
                  - float, air temperature (tair)
                  - float, precipitation (pn)
        Output:   None.
        Features: Writes daily variables (i.e., timestamp, sf, tair and pn) to 
                  the output file in CSV format
        """
        output_line = "%s,%f,%f,%f\n" % (ts, sf, tair, pn)
        try:
            OUT = open(self.output_file, 'a')
            OUT.write(output_line)
        except IOError:
            print "Error: cannot write to file: " 
            print self.output_file
        else:
            OUT.close()
    #
    def writeout(self, f, d):
        """
        Name:     STASH_DATA.writeout
        Input:    - str, file name with path (t)
                  - str, data to be written to file (d)
        Features: Writes new/overwrites existing file with data string
        """
        try:
            OUT = open(f, 'w')
            OUT.write(d)
        except IOError:
            print "Error: cannot write to file: ", f
        else:
            OUT.close()

###############################################################################
## FUNCTIONS
###############################################################################
def add_one_day(dt0):
    """
    Name:     add_one_day
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt1)
    Features: Adds one day to datetime
    """
    dt1 = dt0 + datetime.timedelta(days=1)
    return dt1

###############################################################################
## USER VARIABLES
###############################################################################
# Create a class instance with longitude, latitude and elevation based on your 
# location of interest:
user_lat = 37.7    # degrees north 
user_lon = -122.4  # degrees east  
user_elv = 142.0   # meters AMSV (use 0 if unknown) 
my_class = STASH_DATA(user_lat, user_lon, user_elv)

# Set the data input/output directories for your machine:
my_class.set_cru_cld_dir('/usr/local/share/database/cru/')
my_class.set_cru_pre_dir('/usr/local/share/database/cru/')
my_class.set_cru_tmn_dir('/usr/local/share/database/cru/')
my_class.set_watch_rainf_dir('/usr/local/share/database/watch/netcdf/rainf/')
my_class.set_watch_tair_dir('/usr/local/share/database/watch/netcdf/tair/')
my_class.set_output_dir('/home/user/Desktop/out/')

# Define the sources you want to use for processing the data (cru or watch):
my_class.set_sf_source('cru')
my_class.set_pn_source('watch')
my_class.set_tair_source('watch')

###############################################################################
## MAIN
###############################################################################
##
## Set the date range you want to process
##
start_date = datetime.date(2000, 1, 1)
end_date = datetime.date(2001, 1, 1)

##
## Iterate through time, saving daily data to file
##
cur_date = start_date
while cur_date < end_date:
    my_class.get_daily_status(cur_date, write_out=True)
    cur_date = add_one_day(cur_date)
