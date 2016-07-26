#!/usr/bin/python
# -*- coding: utf-8 -*-
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
import logging
import os

import numpy

from crutsutil import get_cru_lat
from crutsutil import get_cru_lon


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

        self.mode = mode
        self.file_name = ""
        self.year = 0
        self.num_lines = 0.
        self.sf_vec = numpy.array([])
        self.tair_vec = numpy.array([])
        self.pn_vec = numpy.array([])

        # For netCDF handling:
        self._cld_file = None
        self._pre_file = None
        self._tmp_file = None
        self._lat = None
        self._lon = None

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Properties
    # ////////////////////////////////////////////////////////////////////////
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
                self._pre_file = val
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
                self._pre_file = val
            else:
                self.logger.error("CRU cloudiness file must be NetCDF!")
                raise TypeError("CRU cloudiness file must be NetCDF.")
        else:
            self.logger.error("CRU cloudiness file not found!")
            raise OSError("File not found!")

    @property
    def lat(self):
        """Array of latitude values"""
        if self._lat is None:
            if self.cru_pre_file is not None:
                self._lat = get_cru_lat(self.cru_pre_file)
                return self._lat
            elif self.cru_tmp_file is not None:
                self._lat = get_cru_lat(self.cru_tmp_file)
                return self._lat
            elif self.cru_cld_file is not None:
                self._lat = get_cru_lat(self.cru_cld_file)
                return self._lat
            else:
                self.logger.error("No CRU files found for loading lat!")
                raise OSError("No CRU files set! "
                              "Please first assign them to access latitude.")
        else:
            return self._lat

    @lat.setter
    def lat(self, val):
        if isinstance(val, numpy.ndarray):
            self._lat = val
        elif isinstance(val, list):
            self._lat = numpy.array(val)
        else:
            self.logger.error("Expected a list or array of values")
            raise TypeError("Expected a list or array of values")

    @property
    def lon(self):
        """Array of longitude values"""
        if self._lon is None:
            if self.cru_pre_file is not None:
                self._lon = get_cru_lon(self.cru_pre_file)
                return self._lon
            elif self.cru_tmp_file is not None:
                self._lon = get_cru_lon(self.cru_tmp_file)
                return self._lon
            elif self.cru_cld_file is not None:
                self._lon = get_cru_lon(self.cru_cld_file)
                return self._lon
            else:
                self.logger.error("No CRU files found for loading lon!")
                raise OSError("No CRU files set! "
                              "Please first assign them to access longitude.")
        else:
            return self._lon

    @lon.setter
    def lon(self, val):
        if isinstance(val, numpy.ndarray):
            self._lon = val
        elif isinstance(val, list):
            self._lon = numpy.array(val)
        else:
            self.logger.error("Expected a list or array of values")
            raise TypeError("Expected a list or array of values")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
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

    def read_netcdf(self, fname, type="cru"):
        """
        """
        self.file_name = fname


###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("data.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

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
    for lat in my_data.lat:
        for lon in my_data.lon:
            print("(%s, %s)" % (lat, lon))
