#!/usr/bin/python
#
# data.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2016-07-26
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
#
# This file is part of the SPLASH model.
#
# SPLASH is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# SPLASH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
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
import logging
import os

import numpy
from scipy.io import netcdf


###############################################################################
# FUNCTIONS:
###############################################################################
def get_cru_lat(fname):
    """
    Name:     get_cru_lat
    Input:    string, CRU netcdf file (fname)
    Output:   numpy nd.array
    Features: Returns NetCDF data array of latitudes
    """
    lat = numpy.array([])
    if os.path.isfile(fname):
        try:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(fname, "r")
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
        else:
            # Read the latitude data as array:
            lat = f.variables['lat'].data.copy()
        finally:
            f.close()
    else:
        logging.warning("Failed to find file, %s, for reading", fname)
    return lat


def get_cru_lon(fname):
    """
    Name:     get_cru_lon
    Input:    string, CRU netcdf file (fname)
    Output:   numpy nd.array
    Features: Returns NetCDF data array of longitudes
    """
    lon = numpy.array([])
    if os.path.isfile(fname):
        try:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(fname, "r")
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
        else:
            # Read the latitude data as array:
            lon = f.variables['lon'].data.copy()
        finally:
            f.close()
    else:
        logging.warning("Failed to find file, %s, for reading", fname)
    return lon


def get_monthly_cru(fname, ct, v):
    """
    Name:     get_monthly_cru
    Input:    - string, CRU netcdf file (fname)
              - datetime.date, current month datetime object (ct)
              - string, variable of interest (v)
    Output:   numpy nd.array
    Depends:  get_time_index
    Features: Returns 360x720 monthly CRU TS dataset for a given month and
              variable of interest (e.g., cld, pre, tmp)
    """
    if os.path.isfile(fname):
        try:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(fname, "r")
            #
            # Save data for variables of interest:
            # NOTE: for CRU TS 3.21:
            #       variables: 'lat', 'lon', 'time', v
            #       where v is 'tmp', 'pre', 'cld'
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
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
        else:
            # Save the base time stamp:
            bt = datetime.date(1900, 1, 1)
            #
            # Read the time data as array:
            f_time = f.variables['time'].data
            #
            # Find the time index for the current date:
            ti = get_time_index(bt, ct, f_time)
            #
            # Get the spatial data for current time:
            f_data = f.variables[v].data[ti]
            f.close()
            return f_data


def process_cru(v, cru_dir, my_dir):
    """
    Name:     process_cru
    Input:    - string, variable name, e.g., tmp, pre, cld (v)
              - string, directory name for CRU TS data file (cru_dir)
              - string, directory for output files (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF file by month into variable list and data
              set table output files
    Depends:  - writeout
              - add_one_month
              - get_monthly_cru
    """
    # Define the start and end dates you want to process (2002-2006):
    start_date = datetime.date(2002, 1, 1)
    end_date = datetime.date(2007, 1, 1)
    #
    # Set flag for varlist:
    #
    # Prepare var output file:
    #
    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # # Open and read netcdf files in the file directory:
        my_data = get_monthly_cru(cru_dir, cur_date, v)
        (sh_lat, sh_lon) = my_data.shape
        #
        # Prepare data set output file:
        #
        # Iterate through each lat:lon pair
        # * row-major ordering from bottom left
        #  x (longitude): 0...719
        #  y (latitude): 0...359
        for y in xrange(sh_lat):
            for x in xrange(sh_lon):
                # Calc station ID:
                st_id = 720*y + x
                station_parts = ('HDG', st_id)
                print("%s, %s" % station_parts)
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Turn off varlist flag after processing first month:
        #
        # Increment cur_date:


def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime.date, base timestamp
              - datetime.date, current timestamp
              - numpy nd.array, days since base timestamp
    Output:   int
    Features: Finds the index in an array of CRU TS days for a given timestamp
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
    aot = numpy.append(aot, [dt, ])
    #
    # Find the first index of dt in the sorted array:
    idx = numpy.where(numpy.sort(aot) == dt)[0][0]
    return idx


if __name__ == '__main__':
    pass
