#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# data_grid.py
#
# LAST UPDATED: 2016-02-13
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
# led algorithms for simulating habitats (SPLASH): Robust indices of radiation,
# evapotranspiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)

###############################################################################
# IMPORT MODULES:
###############################################################################
import datetime
import glob
import logging
import os.path

import numpy
from scipy.io import netcdf

from utilities import get_x_y


###############################################################################
# CLASSES:
###############################################################################
class DATA_G:
    """
    Name:     DATA_G
    Features: This class handles the file IO for reading gridded CRU TS data.
    History   Version 1.0.0-dev
              - added longitude and latitude class variables [16.02.13]
              - created get lon lat function [16.02.13]
    """
    error_val = numpy.inf

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     DATA_G.__init__
        Inputs:   None.
        Features: Initialize class variables.
        """
        # Initialize data to None:
        self.date = None
        self.elv_file = None
        self.cld_file = None
        self.pre_file = None
        self.tmp_file = None
        self.tair = None
        self.pre = None
        self.sf = None
        self.elevation = None
        self.latitude = None
        self.longitude = None

        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("DATA_G class called")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def add_one_month(self, dt0):
        """
        Name:     DATA_G.add_one_month
        Input:    datetime.date
        Output:   datetime.date
        Features: Adds one month to datetime
        Ref:      A. Balogh (2010), ActiveState Code
                  http://code.activestate.com/recipes/577274-subtract-or-add-a-
                  month-to-a-datetimedate-or-datet/
        """
        dt1 = dt0.replace(day=1)
        dt2 = dt1 + datetime.timedelta(days=32)
        dt3 = dt2.replace(day=1)
        return dt3

    def add_one_year(self, dt0):
        """
        Name:     DATA_G.add_one_year
        Input:    datetime.date
        Output:   datetime.date
        Features: Adds one year to the datetime preserving calendar date, if it
                  exists, otherwise uses the following day (i.e., February 29
                  becomes March 1)
        Ref:      G. Rees (2013) Stack Overflow
                  http://stackoverflow.com/questions/15741618/add-one-year-in-
                  current-date-python
        """
        try:
            return dt0.replace(year=dt0.year + 1)
        except ValueError:
            return dt0 + (datetime.date(dt0.year + 1, 1, 1) -
                          datetime.date(dt0.year, 1, 1))

    def find_cru_files(self, path):
        """
        Features: Searches for the four CRU TS variable files (i.e., elv, cld,
                  pre, tmp) within a single directory
        Depends:  get_cru_file
        """
        cld_file = self.get_cru_file(path, 'cld')
        if os.path.isfile(cld_file):
            self.logger.info("found CRU TS cloudiness file %s", cld_file)
            self.cld_file = cld_file
        else:
            self.logger.warning("failed to load CRU TS cloudiness file")
            self.cld_file = None

        elv_file = self.get_cru_file(path, 'elv')
        if os.path.isfile(elv_file):
            self.logger.info("found CRU TS elevation file %s", elv_file)
            self.elv_file = elv_file
        else:
            self.logger.warning("failed to load CRU TS elevation file")
            self.elv_file = None

        pre_file = self.get_cru_file(path, 'pre')
        if os.path.isfile(pre_file):
            self.logger.info("found CRU TS precipitation file %s", pre_file)
            self.pre_file = pre_file
        else:
            self.logger.warning("failed to load CRU TS precipitation file")
            self.pre_file = None

        tmp_file = self.get_cru_file(path, 'tmp')
        if os.path.isfile(tmp_file):
            self.logger.info("found CRU TS temperature file %s", tmp_file)
            self.tmp_file = tmp_file
        else:
            self.logger.warning("failed to load CRU TS temperature file")
            self.tmp_file = None

    def get_cru_file(self, path, voi):
        """
        Name:     DATA_G.get_cru_file
        Input:    - str, directory path for CRU data files (path)
                  - str, variable of interest (voi)
        Output:   str OR list of file names
        Features: Returns the CRU TS file for given variable of interest
        """
        # Read through all files within the paths for voi:
        my_file = None
        my_pattern = os.path.join(path, "*%s*.*" % (voi))
        my_files = glob.glob(my_pattern)

        if my_files:
            if len(my_files) > 1:
                self.logger.warning("Found %d files!", len(my_files))
            else:
                my_file = my_files[0]
                self.logger.info("found file %s", my_file)
        else:
            self.logger.warning("Found 0 files!")

        return my_file

    def get_month_days(self, ts):
        """
        Name:     DATA_G.get_month_days
        Input:    datetime date
        Output:   int, days in the month
        Features: Returns the number of days in the month
        Depends:  add_one_month
        """
        ts1 = ts.replace(day=1)
        ts2 = self.add_one_month(ts1)
        dts = (ts2 - ts1).days
        self.logger.debug("month has %d days", dts)

        return dts

    def get_monthly_cru(self, ct, v):
        """
        Name:     DATA_G.get_monthly_cru
        Input:    - datetime date, current month datetime object (ct)
                  - str, variable of interest (v)
        Output:   numpy nd.array
        Features: Returns 360x720 monthly CRU TS dataset for a given month and
                  variable of interest (e.g., cld, pre, tmp)
                  NODATA = error_val
        Depends:  - get_cru_file
                  - get_time_index
        """
        self.logger.debug("reading %s for month %s" % (v, ct))

        if v == 'tmp':
            my_file = self.tmp_file
        elif v == 'pre':
            my_file = self.pre_file
        elif v == 'cld':
            my_file = self.cld_file

        if my_file:
            # Open netCDF file for reading:
            self.logger.debug("opening NetCDF file %s", my_file)
            f = netcdf.NetCDFFile(my_file, "r")

            # Save data for variables of interest:
            # NOTE: for CRU TS 3.2:
            #       variables: 'lat', 'lon', 'time', v
            #       where v is 'tmp', 'pre', 'cld'
            # LAT:  -89.75 -- 89.75
            # LON:  -179.75 -- 179.75
            # TIME: units: days since 1900-1-1
            #       shape: (1344,)
            #       values: mid-day of each month (e.g., 15th or 16th)
            # DATA: 'cld' units = %
            #       'pre' units = mm
            #       'tmp' units = deg. C
            #       Missing value = 9.96e+36

            # Save the base time stamp:
            bt = datetime.date(1900, 1, 1)

            # Read the time data as array:
            f_time = f.variables['time'].data.copy()

            # Find the time index for the current date:
            ti = self.get_time_index(bt, ct, f_time)

            # Get the spatial data for current time:
            f_var = f.variables[v]
            f_noval = f_var.missing_value
            f_temp = f_var.data[ti]
            f_data = numpy.copy(f_temp)
            f_var = None
            f.close()

            self.logger.debug("setting error value")
            noval_idx = numpy.where(f_data == f_noval)
            f_data[noval_idx] *= 0.0
            f_data[noval_idx] += self.error_val

            self.logger.debug("finished reading %s for month %s" % (v, ct))
            return f_data

    def get_time_index(self, bt, ct, aot):
        """
        Name:     get_time_index
        Input:    - datetime.date, base timestamp (bt)
                  - datetime.date, current timestamp (ct)
                  - numpy.ndarray, array of days since base timestamp (aot)
        Output:   int, time index for given month
        Features: Finds the index in an array of days for a given timestamp
        """
        # For CRU TS 3.2, the aot is indexed for mid-month days, e.g. 15--16th
        # therefore, to make certain that ct index preceeds the index for the
        # correct month in aot, make the day of the current month less than
        # the 15th or 16th (i.e., replace day with '1'):
        ct = ct.replace(day=1)

        # Calculate the time difference between ct and bt:
        dt = (ct - bt).days
        if dt < 0:
            self.logger.warning("Month %s preceeds CRU base date!", ct)
            idx = None
        else:
            try:
                idx = numpy.where(aot > dt)[0][0]
            except IndexError:
                self.logger.warning("Month %s is out of bounds!", ct)
                idx = None
            else:
                self.logger.debug("Found index %d for month %s" % (idx, ct))
                return idx

    def get_year_days(self, ts):
        """
        Name:     DATA_G.get_year_days
        Input:    datetime date
        Output:   int
        Features: Returns the total number of days in the year
        Depends:  add_one_year
        """
        ts1 = datetime.date(ts.year, 1, 1)
        ts2 = self.add_one_year(ts1)
        ryd = (ts2 - ts1).days
        self.logger.debug("year has %d days", ryd)

        return ryd

    def print_vals(self, lon, lat):
        """
        Name:     DATA_G.print_vals
        Inputs:   - float, longitude, degrees (lon)
                  - float, latitude, degrees (lat)
        Outputs:  None.
        Features: Prints the four variables (i.e., elv, pre, tair, and sf) for
                  a given location
        Depends:  get_x_y
        """
        x, y = get_x_y(lon, lat)
        print("Date: %s" % (self.date))
        print("Longitude: %0.4f degrees (%d)" % (lon, x))
        print("Latitude: %0.4f degrees (%d)" % (lat, y))
        print("Elevation: %0.4f m" % (self.elevation[y, x]))
        print("Mean air temperature: %0.6f deg. C" % (self.tair[y, x]))
        print("Precipitation: %0.6f mm/d" % (self.pre[y, x]))
        print("Sunshine fraction: %0.6f" % (self.sf[y, x]))

    def read_elv(self):
        """
        Name:     DATA_G.read.elv
        Inputs:   None.
        Outputs:  None.
        Features: Reads elevation data from file and sets good and no value
                  indexes
        """
        if self.elv_file:
            try:
                self.logger.debug("reading elevation data")
                f = numpy.loadtxt(self.elv_file)
            except:
                self.logger.exception("failed to read elevation data")
                self.elevation = None
                self.good_idx = (numpy.array([]), numpy.array([]))
                self.noval_idx = (numpy.array([]), numpy.array([]))
            else:
                self.logger.debug("read %d values", f.size)
                self.good_idx = numpy.where(f != -999.0)
                self.noval_idx = numpy.where(f == -999.0)
                f[self.noval_idx] *= 0.0
                f[self.noval_idx] += self.error_val
                self.elevation = f
        else:
            self.logger.warnig("no elevation file found!")
            self.elevation = None
            self.good_idx = (numpy.array([]), numpy.array([]))
            self.noval_idx = (numpy.array([]), numpy.array([]))

    def read_lon_lat(self):
        """
        Name:     DATA_G.read_lon_lat
        Inputs:   None.
        Outputs:  None.
        Features: Retrieves latitude array from CRU temperature netCDF file
        """
        self.logger.debug("retrieving lon and lat arrays from CRU file")
        my_file = self.tmp_file
        if my_file:
            f = netcdf.NetCDFFile(my_file, "r")
            f_lat = f.variables['lat'].data.copy()
            f_lon = f.variables['lon'].data.copy()
            f.close()
        else:
            self.logger.debug("Please set CRU files before loading lon/lat")
            f_lat = numpy.array([])
            f_lon = numpy.array([])
        self.logger.debug("Read %d latitude values", f_lat.size)
        self.logger.debug("Read %d longitude values", f_lon.size)
        self.latitude = f_lat
        self.longitude = f_lon

    def read_monthly_clim(self, m):
        """
        Name:     DATA_G.read_monthly_clim
        Inputs:   datetime.date, month of interest (m)
        Features: Reads monthly climatology (i.e., temperature, precipitation,
                  and cloudiness), converts cloudiness to sunshine fraction and
                  monthly precipitation to daily fraction
        Depends:  - get_monthly_cru
                  - set_date
        """
        # Check to see if new monthly data needs to be processed:
        to_process = False
        self.logger.debug("checking whether to process month... ")
        if self.date:
            # Check to see if this is a new month:
            if not self.year == m.year or not self.month == m.month:
                to_process = True
        else:
            to_process = True
        self.logger.debug("... %s", to_process)

        if to_process:
            # Set the date:
            self.set_date(m)

            # Reset the data arrays:
            self.logger.debug("initializing climate arrays")
            self.tair = numpy.zeros(shape=(360, 720))
            self.pre = numpy.zeros(shape=(360, 720))
            self.sf = numpy.zeros(shape=(360, 720))

            # Read monthly data:
            self.logger.debug("reading monthly climatology")
            tmp = self.get_monthly_cru(m, 'tmp')
            pre = self.get_monthly_cru(m, 'pre')
            cld = self.get_monthly_cru(m, 'cld')

            # Update good and noval indexes:
            self.logger.debug("updating good and no-value indexes")
            self.noval_idx = numpy.where(
                (self.elevation == self.error_val) | (tmp == self.error_val) |
                (pre == self.error_val) | (cld == self.error_val))
            self.good_idx = numpy.where(
                (self.elevation != self.error_val) & (tmp != self.error_val) &
                (pre != self.error_val) & (cld != self.error_val))

            # Convert cloudiness to fractional sunshine, sf
            self.logger.debug("converting cloudiness to sunshine fraction")
            sf = numpy.copy(cld)
            sf[self.good_idx] *= 1e-2
            sf[self.good_idx] *= -1.0
            sf[self.good_idx] += 1.0

            # Convert monthly precip to daily fraction & clip missing to zero:
            self.logger.debug("converting monthly to daily precipitation")
            pre[self.good_idx] /= float(self.nm)

            # Save data:
            self.tair = tmp
            self.pre = pre
            self.sf = sf

    def set_date(self, date):
        """
        Name:     DATA_G.set_date
        Input:    datetime.date, date of interest (date)
        Output:   None
        Features: Sets the days and years for a given date
        Depends:  - get_month_days
                  - get_year_days
        """
        self.logger.debug("setting date to %s", date)
        self.date = date
        self.year = date.year
        self.month = date.month
        self.n = date.timetuple().tm_yday
        self.nm = self.get_month_days(date)
        self.ny = self.get_year_days(date)

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Instantiating logging handler and record format:
    root_handler = logging.StreamHandler()
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    cru_dir = "/usr/local/share/data/cru"
    my_class = DATA_G()
    my_class.find_cru_files(cru_dir)
    my_class.read_elv()
    my_class.read_lon_lat()

    my_day = 172
    my_year = 2000
    my_date = datetime.date(my_year, 1, 1)
    my_date += datetime.timedelta(days=(my_day-1))

    my_class.read_monthly_clim(my_date)
    root_logger.info("read %d precipitation", my_class.pre.size)
    root_logger.info("read %d sunshine fraction", my_class.sf.size)
    root_logger.info("read %d air temperature", my_class.tair.size)
