#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# data.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2016-07-27
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

# @TODO: maybe move the netcdf handling to SPLASH_DATA class and make it an
#        attribute to DATA

###############################################################################
# IMPORT MODULES:
###############################################################################
import datetime
import logging
import os

import numpy

from crutsutil import add_one_month
from crutsutil import get_cru_elv
from crutsutil import get_cru_lat
from crutsutil import get_cru_lon
from crutsutil import get_monthly_cru


###############################################################################
# CLASSES
###############################################################################
class DATA(object):
    """
    Name:     DATA
    Features: This class handles the file IO for reading and writing data.
    History:  Version 1.1.0-dev
              - added logging statements [16.02.05]
              - added netCDF file handing [16.07.26]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, mode="txt"):
        """
        Name:     DATA.__init__
        Input:    str, operation mode (mode)
                  > 'txt' for plain text files (point data)
                  > 'latlon' for NetCDF files (gridded data)
        Features: Initialize empty class variables
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("DATA class called; using mode %s", mode)

        self._year = 0
        self.mode = mode
        self.file_name = ""
        self.num_lines = 0.
        self.sf_vec = numpy.array([])
        self.tair_vec = numpy.array([])
        self.pn_vec = numpy.array([])

        # For netCDF handling:
        self._cld_file = None      # CRU TS cloudiness netCDF file
        self._pre_file = None      # CRU TS precipitation netCDF file
        self._tmp_file = None      # CRU TS mean air temperature netCDF file
        self._elv_file = None      # CRUT TS elevation data file
        self._latitude = None      # array of latitudes
        self._longitude = None     # array of longitudes
        self._lat = None           # current latitude
        self._lon = None           # current longitude
        self._elv = None           # current elevation

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Properties
    # ////////////////////////////////////////////////////////////////////////
    @property
    def cru_elv_file(self):
        """CRU TS 3.00 elevation data file"""
        return self._elv_file

    @cru_elv_file.setter
    def cru_elv_file(self, val):
        if os.path.isfile(val):
            if os.path.splitext(val)[1] == ".dat":
                self.logger.info("Setting CRU ELV grid file to %s", val)
                self._elv_file = val
            else:
                self.logger.error("CRU elevation must be a data file!")
                raise TypeError("CRU elevation must be a data file!")
        else:
            self.logger.error("CRU elevation file not found!")
            raise OSError("File not found!")

    @property
    def cru_pre_file(self):
        """NetCDF CRU Monthly Precipitation File"""
        return self._pre_file

    @cru_pre_file.setter
    def cru_pre_file(self, val):
        if os.path.isfile(val):
            if os.path.splitext(val)[1] == ".nc":
                self.logger.info("Setting CRU PRE file to %s", val)
                self._pre_file = val
            else:
                self.logger.error("CRU precipitation file must be NetCDF!")
                raise TypeError("CRU precipitation file must be NetCDF.")
        else:
            self.logger.error("CRU precipitation file not found!")
            raise OSError("File not found!")

    @property
    def cru_tmp_file(self):
        """NetCDF CRU Monthly Air Temperature File"""
        return self._tmp_file

    @cru_tmp_file.setter
    def cru_tmp_file(self, val):
        if os.path.isfile(val):
            if os.path.splitext(val)[1] == ".nc":
                self.logger.info("Setting CRU TMP file to %s", val)
                self._tmp_file = val
            else:
                self.logger.error("CRU air temperature file must be NetCDF!")
                raise TypeError("CRU air temperature file must be NetCDF.")
        else:
            self.logger.error("CRU air temperature file not found!")
            raise OSError("File not found!")

    @property
    def cru_cld_file(self):
        """NetCDF CRU Monthly Cloudiness File"""
        return self._cld_file

    @cru_cld_file.setter
    def cru_cld_file(self, val):
        if os.path.isfile(val):
            if os.path.splitext(val)[1] == ".nc":
                self.logger.info("Setting CRU CLD file to %s", val)
                self._cld_file = val
            else:
                self.logger.error("CRU cloudiness file must be NetCDF!")
                raise TypeError("CRU cloudiness file must be NetCDF.")
        else:
            self.logger.error("CRU cloudiness file not found!")
            raise OSError("File not found!")

    @property
    def latitude(self):
        """Array of latitude values"""
        if self._latitude is None:
            if self.cru_pre_file is not None:
                self._latitude = get_cru_lat(self.cru_pre_file)
                return self._latitude
            elif self.cru_tmp_file is not None:
                self._latitude = get_cru_lat(self.cru_tmp_file)
                return self._latitude
            elif self.cru_cld_file is not None:
                self._latitude = get_cru_lat(self.cru_cld_file)
                return self._latitude
            else:
                self.logger.error("No CRU files found for loading lat!")
                raise OSError("No CRU files set! "
                              "Please first assign them to access latitude.")
        else:
            return self._latitude

    @latitude.setter
    def latitude(self, val):
        if isinstance(val, numpy.ndarray):
            self._latitude = val
        elif isinstance(val, list):
            self._latitude = numpy.array(val)
        else:
            self.logger.error("Expected a list or array of values")
            raise TypeError("Expected a list or array of values")

    @property
    def longitude(self):
        """Array of longitude values"""
        if self._longitude is None:
            if self.cru_pre_file is not None:
                self._longitude = get_cru_lon(self.cru_pre_file)
                return self._longitude
            elif self.cru_tmp_file is not None:
                self._longitude = get_cru_lon(self.cru_tmp_file)
                return self._longitude
            elif self.cru_cld_file is not None:
                self._longitude = get_cru_lon(self.cru_cld_file)
                return self._longitude
            else:
                self.logger.error("No CRU files found for loading lon!")
                raise OSError("No CRU files set! "
                              "Please first assign them to access longitude.")
        else:
            return self._longitude

    @longitude.setter
    def longitude(self, val):
        if isinstance(val, numpy.ndarray):
            self._longitude = val
        elif isinstance(val, list):
            self._longitude = numpy.array(val)
        else:
            self.logger.error("Expected a list or array of values")
            raise TypeError("Expected a list or array of values")

    @property
    def year(self):
        """The process year"""
        return self._year

    @year.setter
    def year(self, val):
        if isinstance(val, int):
            self._year = val
        elif isinstance(val, str):
            try:
                tmp = int(val)
            except:
                raise
            else:
                self._year = tmp
        else:
            self.logger.error("Expected an integer for year.")
            raise TypeError("Year is not an integer!")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def check_cru_files(self):
        """Checks that all CRU files have been assigned"""
        if self._cld_file is None:
            self.logger.error("No CRU TS cloudiness file found!")
            raise AttributeError("CRU TS cloudiness file not assigned!")
        if self._pre_file is None:
            self.logger.error("No CRU TS precipitation file found!")
            raise AttributeError("CRU TS precipitation file not assigned!")
        if self._tmp_file is None:
            self.logger.error("No CRU TS air temperature file found!")
            raise AttributeError("CRU TS air temperature file not assigned!")
        if self._elv_file is None:
            self.logger.warning("No CRU TS elevation file found!")
            print("Warning: CRU TS elevation is not set!")

    def get_annual_cru(self, year, lat, lon):
        """
        Name:     DATA.get_annual_cru
        Input:    - int, process year (year)
                  - float, latitude (lat)
                  - float, longitude (lon)
        Output:   None
        Features: Reads all three input variables (sf, tair, and pn) for
                  a single year from their respective netCDF files.
        """
        # Make certain CRU TS files have been assigned and set elevation:
        self.check_cru_files()
        self.set_elevation(lat, lon)

        # Find the time index for each CRU TS file
        # * should be the same if they are all from the same CRU TS release
        yc = datetime.date(year, 1, 1)       # current year
        ye = datetime.date(year + 1, 1, 1)   # end year
        day_cld = numpy.array([])
        day_pre = numpy.array([])
        day_tmp = numpy.array([])
        while yc < ye:
            try:
                cld_mo = get_monthly_cru(
                    self.cru_cld_file, yc, lat, lon, 'cld')
            except:
                self.logger.exception("Failed to read from cloudiness file")
                raise
            else:
                self.logger.debug("cloudiness value = %s", cld_mo)
                cld_day = self.monthly_to_daily(cld_mo, yc)
                day_cld = numpy.append(day_cld, [cld_day, ])

            try:
                pre_mo = get_monthly_cru(
                    self.cru_pre_file, yc, lat, lon, 'pre')
            except:
                self.logger.exception(
                    "Failed to read from precipitation file")
                raise
            else:
                self.logger.debug("precipitation value = %s", pre_mo)
                pre_day = self.monthly_to_daily(pre_mo, yc)
                day_pre = numpy.append(day_pre, [pre_day, ])

            try:
                tmp_mo = get_monthly_cru(
                    self.cru_tmp_file, yc, lat, lon, 'tmp')
            except:
                self.logger.exception(
                    "Failed to read from air temperature file")
                raise
            else:
                self.logger.debug("air temperature value = %s", tmp_mo)
                tmp_day = self.monthly_to_daily(tmp_mo, yc)
                day_tmp = numpy.append(day_tmp, [tmp_day, ])

            yc = add_one_month(yc)
        self.logger.debug("ready with %d cloudiness points", len(day_cld))
        self.logger.debug("ready with %d precipitation points", len(day_pre))
        self.logger.debug("ready with %d temperature points", len(day_tmp))

    def monthly_to_daily(self, mo_ds, mo_ts, method="const"):
        """
        Name:     DATA.monthly_to_daily
        Inputs:   - float, monthly dataset (mo_ds)
                  - datetime.date, monthly timestamp (mo_ts)
                  - [optional] str, quasi-daily method (method)
        Outputs:  numpy.ndarray
        Features: Returns quasi-daily array based on the given monthly value
        """
        tcur = mo_ts.replace(day=1)
        tend = add_one_month(tcur)
        mo_days = (tend - datetime.timedelta(days=1)).day
        if method == "const":
            day_ds = numpy.ones((mo_days, ), dtype=numpy.float)
            day_ds *= mo_ds
        else:
            self.logger.error("No method called %s", method)
            raise AttributeError("Method %s is unavailable" % method)
        return day_ds

    def read_csv(self, fname, y=-1):
        """
        Name:     DATA.read_csv
        Input:    - str, input CSV filename (fname)
                  - int, year (y)
        Output:   None
        Features: Reads all three daily input variables (sf, tair, and pn) for
                  a single year from a CSV file that includes a headerline.
        """
        self.file_name = fname

        try:
            data = numpy.loadtxt(fname,
                                 dtype={'names': ('sf', 'tair', 'pn'),
                                        'formats': ('f4', 'f4', 'f4')},
                                 delimiter=',',
                                 skiprows=1)
        except IOError:
            self.logger.exception("could not read input file %s", fname)
            raise
        else:
            self.sf_vec = data['sf']
            self.tair_vec = data['tair']
            self.pn_vec = data['pn']
            self.num_lines = data.shape[0]

            if y == -1:
                if data.shape[0] == 366:
                    self.year = 2000
                elif data.shape[0] == 365:
                    self.year = 2001
            else:
                self.year = y

    def read_txt(self, fname, var, y=-1):
        """
        Name:     DATA.read_txt
        Input:    - str, input text file (fname)
                  - str, variable name (i.e., 'pn', 'sf', 'tair')
                  - int, year (y)
        Output:   None.
        Features: Reads plain text file (no header) into one of daily input
                  arrays.
        """
        # Add filename to list:
        if not isinstance(self.file_name, list):
            self.file_name = []
        self.file_name.append(fname)

        try:
            data = numpy.loadtxt(fname, dtype='f4')
        except IOError:
            self.logger.exception("could not read input file %s", fname)
            raise
        else:
            if var == 'sf':
                self.sf_vec = data
            elif var == 'pn':
                self.pn_vec = data
            elif var == 'tair':
                self.tair_vec = data
            else:
                self.logger.error("variable %s undefined!", var)
                raise ValueError("Unrecognized variable in read_txt")

            # Add line numbers to list:
            if not isinstance(self.num_lines, list):
                self.num_lines = []
            self.num_lines.append(data.shape[0])

            if y == -1:
                if data.shape[0] == 366:
                    self.year = 2000
                elif data.shape[0] == 365:
                    self.year = 2001
            else:
                self.year = y

    def set_elevation(self, lat, lon):
        """
        Name:     DATA.set_elevation
        Input:    - float, latitude (lat)
                  - float, longitude (lon)
        Output:   None.
        Features: Sets elevation based on latitude and longitude.
        """
        self._elv = None
        if self._elv_file is not None:
            if self._latitude is not None:
                if self._longitude is not None:
                    # Find the array indexes for the longitude and latitude:
                    if lat in self.latitude:
                        lat_idx = numpy.where(self.latitude == lat)[0]
                        self.logger.debug("lat index = %d", lat_idx)
                    else:
                        self.logger.error("Could not index latitude %s", lat)
                        raise ValueError(
                            "Latitude %s not found in CRU TS file!" % lat)

                    if lon in self.longitude:
                        lon_idx = numpy.where(self.longitude == lon)[0]
                        self.logger.debug("lon index = %d", lon_idx)
                    else:
                        self.logger.error("Could not index longitude %s", lon)
                        raise ValueError(
                            "Longitude %s not found in CRU TS file!" % lon)

                    # Extract gridded elevation based on lon and lat indexes
                    elevation = get_cru_elv(self.cru_elv_file)
                    elv = elevation[lat_idx, lon_idx]
                    self._elv = elv
                    self.logger.debug("elv = %f", elv)
                else:
                    self.logger.warning("Longitude is not set!")
            else:
                self.logger.warning("Latitude is not set!")
        else:
            self.logger.warning("Elevation file is not set!")


###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("temp.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # Set CRU TS file paths:
    my_data = DATA(mode="latlon")
    my_data.cru_cld_file = os.path.join(
        os.path.expanduser("~"), "Data", "cru_ts",
        "cru_ts3.22.1901.2013.cld.dat.nc")
    my_data.cru_pre_file = os.path.join(
        os.path.expanduser("~"), "Data", "cru_ts",
        "cru_ts3.22.1901.2013.pre.dat.nc")
    my_data.cru_tmp_file = os.path.join(
        os.path.expanduser("~"), "Data", "cru_ts",
        "cru_ts3.22.1901.2013.tmp.dat.nc")
    my_data.cru_elv_file = os.path.join(
        os.path.expanduser("~"), "Data", "cru_ts",
        "cru_ts3.00_halfdeg.elv.grid.dat")

    my_lat = my_data.latitude[260]
    my_lon = my_data.longitude[200]
    my_data.get_annual_cru(2000, my_lat, my_lon)
