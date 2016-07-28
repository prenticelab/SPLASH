#!/usr/bin/python
#
# data.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2016-07-28
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
def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3


def cld_to_sf(cld):
    """
    Name:     cld_to_sf
    Inputs:   float or numpy.ndarray, CRU TS cloudiness (cld)
    Outputs:  float, estimate of fractional sunshine hours
    Features: Converts CRU TS cloudiness to fractional sunshine hours, Sf,
              assuming that Sf is the complement of cloudiness
    """
    cld /= 100.0    # unitless
    sf = 1.0 - cld  # complement of cloudiness
    return sf


def get_cru_elv(fname):
    """
    Name:     get_cru_elv
    Input:    string, input/output file directory (d)
    Output:   None.
    Features: Return array of CRU TS 3.00 elevation integers from the .dat file

              NOTE: data is read into an array with shape (360, 720)
                    'lat' goes from -89.75 to 89.75 (360 points south to north)
                    'lon' goes from -179.75 to 179.75 (720 points east to west)
                    'nodata' value == -999.0
                    'max' value is only 5734, so array can be stored as int16
    """
    elv = numpy.zeros((360, 720), numpy.int16)

    # Open and read dat file:
    if fname is not None and os.path.isfile(fname):
        try:
            elv = numpy.loadtxt(fname, numpy.int16)
        except:
            logging.exception("Failed to read CRU elevation file!")
    else:
        logging.exception("CRU elevation file does not exist!")

    return elv


def get_cru_lat(fname):
    """
    Name:     get_cru_lat
    Input:    string, CRU netcdf file (fname)
    Output:   numpy nd.array
    Features: Returns NetCDF data array of latitudes
    """
    if os.path.isfile(fname):
        try:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(fname, "r")
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
            raise
        else:
            try:
                # Read the latitude data as array:
                lat = f.variables['lat'].data.copy()
            except:
                logging.exception(
                    "Failed to read 'lat' variable from %s", fname)
                raise
        finally:
            f.close()
    else:
        logging.error("Failed to find file, %s, for reading", fname)
        raise OSError("Could not find file!")

    return lat


def get_cru_lon(fname):
    """
    Name:     get_cru_lon
    Input:    string, CRU netcdf file (fname)
    Output:   numpy nd.array
    Features: Returns NetCDF data array of longitudes
    """
    if os.path.isfile(fname):
        try:
            # Open netCDF file for reading:
            f = netcdf.NetCDFFile(fname, "r")
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
            raise
        else:
            try:
                # Read the latitude data as array:
                lon = f.variables['lon'].data.copy()
            except:
                logging.exception(
                    "Failed to read 'lon' variable from %s", fname)
                raise
        finally:
            f.close()
    else:
        logging.error("Failed to find file, %s, for reading", fname)
        raise OSError("Could not find file!")

    return lon


def get_cru_time(fname):
    """
    Name:     get_cru_time
    Input:    string, CRU netcdf file (fname)
    Output:   numpy nd.array
    Features: Returns NetCDF data array of time differences
    """
    if os.path.isfile(fname):
        try:
            f = netcdf.NetCDFFile(fname, "r")
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
            raise
        else:
            try:
                tm = f.variables['time'].data.copy()
            except:
                logging.exception(
                    "Failed to read 'time' variable from %s", fname)
                raise
        finally:
            f.close()
    else:
        logging.error("Failed to find file, %s, for reading", fname)
        raise OSError("Could not find file!")

    return tm


def get_time_index(ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime.date, current timestamp
              - numpy.ndarray, array of time differentials (aot)
    Output:   int
    Features: Finds the index in an array of CRU TS days for a given timestamp
    """
    # CRU TS 'time' variable units are given as "days since 1900-1-1"
    bt = datetime.date(1900, 1, 1)

    # For CRU TS 3.21, the aot is indexed for mid-month days, e.g. 15--16th
    # therefore, to make certain that ct index preceeds the index for the
    # correct month in aot, make the day of the current month less than
    # the 15th or 16th (i.e., replace day with '1'):
    ct = ct.replace(day=1)
    dt = (ct - bt).days                  # datetime difference in days
    idx = numpy.searchsorted(aot, dt)    # numpy array index
    return idx


def get_monthly_cru(fname, ct, lat, lon, v):
    """
    Name:     get_monthly_cru
    Input:    - string, CRU netcdf file (fname)
              - datetime.date, current datetime (ct)
              - float, latitude (lat)
              - float, longitude (lon)
              - string, variable of interest (v)
    Output:   float
    Depends:  get_time_index
    Features: Returns monthly CRU TS value for a given month and variable of
              interest (e.g., cld, pre, tmp)
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
            #       shape = ('time', 'lat', 'lon')
            #       'cld' units = %
            #       'pre' units = mm
            #       'tmp' units = deg. C
            #       Missing value = 9.96e+36
        except:
            logging.exception("Failed to open file, %s, for reading", fname)
            raise
        else:
            try:
                # Read the time data as array:
                f_time = f.variables['time'].data
            except:
                logging.exception(
                    "Failed to read 'time' variable from %s", fname)
                raise
            else:
                # Find the time index for the current date:
                ti = get_time_index(ct, f_time)

                # Find lon/lat indexes:
                try:
                    f_lon = f.variables['lon'].data
                    f_lat = f.variables['lat'].data
                except:
                    logging.exception(
                        "Failed to read lon/lat variables from %s", fname)
                    raise
                else:
                    if lat in f_lat:
                        ilat = numpy.where(f_lat == lat)[0]
                        logging.debug("lat index = %d", ilat)
                    else:
                        logging.error("Could not index latitude %s", lat)
                        raise ValueError(
                            "Latitude %s not found in CRU TS file!" % lat)
                    if lon in f_lon:
                        ilon = numpy.where(f_lon == lon)[0]
                        logging.debug("lon index = %d", ilon)
                    else:
                        logging.error("Could not index longitude %s", lon)
                        raise ValueError(
                            "Longitude %s not found in CRU TS file!" % lon)

                    # Get data value:
                    try:
                        # Get the spatial data for current time:
                        f_var = f.variables[v]
                    except:
                        logging.exception(
                            "Failed to read '%s' variable from %s" % (
                                v, fname))
                        raise
                    else:
                        f_data = f_var.data[ti].copy()
                        dp = f_data[ilat, ilon]
                        try:
                            f_missing = f_var.missing_value
                        except:
                            logging.warning("Failed to retrieve missing value "
                                            "for variable %s", v)
                            f_missing = None
                        else:
                            if dp == f_missing:
                                dp = numpy.nan
                        f_time = None
                        f_var = None
                        f_lon = None
                        f_lat = None
                        f.close()
                        return dp
    else:
        logging.error("Failed to open file, %s, for reading", fname)
        raise OSError("File not found!")


def pre_to_pn(pre, cm):
    """
    Name:     pre_to_pn
    Inputs:   - float or numpy.ndarray, monthly precipitation (pre)
              - datetime.date, current month (cm)
    Outputs:  float/numpy.ndarray
    Features: Converts CRU TS precipitation (mm/mo) to daily precipitation
              (mm/d)
    """
    # Calculate the number of days in the month:
    tcur = cm.replace(day=1)
    tend = add_one_month(tcur)
    nm = (tend - datetime.timedelta(days=1)).day

    # Convert mm/mo to mm/d
    pn = pre / nm
    return pn
