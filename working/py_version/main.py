#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# main.py
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
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 00. v.1.0 release [16.02.05]
# 01.
#
# -----
# todo:
# -----
# 1. create a function in DATA class to save daily, monthly, and
#    annual quantities & write out
# 2. upate plot commands

###############################################################################
# IMPORT MODULES
###############################################################################
import os
import logging
import sys

from data import DATA
from splash import SPLASH

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.WARNING)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("main.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    method = "point"

    if method == "point":
        my_data = DATA(mode='txt')
        my_file = '../../data/example_data.csv'
        my_data.read_csv(my_file)
        output_dir = os.path.join(os.path.expanduser("~"), "Data", "splash")
        output_file = "splash_out_2000sf.txt"
        output_path = os.path.join(output_dir, output_file)
        try:
            f = open(output_path, 'w')
            tmp = SPLASH(0, 0)
            f.write(tmp._header)
        except:
            raise
        else:
            f.close()
            tmp = None
        my_class = SPLASH(37.7, 142.0)
        my_class.spin_up(my_data)

        # Loop through a year:
        for i in range(my_data.npoints):
            # Get preceding soil moisture status:
            if i == 0:
                wn = my_class.wn_vec[-1]
            else:
                wn = my_class.wn_vec[i-1]

            # Calculate soil moisture and runoff:
            my_class.run_one_day(n=i+1,
                                 y=my_data.year,
                                 wn=wn,
                                 sf=my_data.sf_vec[i],
                                 tc=my_data.tair_vec[i],
                                 pn=my_data.pn_vec[i])
            my_class.wn_vec[i] = my_class.wn

            # Save daily results to file:
            if os.path.isfile(output_path):
                try:
                    f = open(output_path, 'a')
                    f.write(my_class.data_str)
                except:
                    logging.exception("Write failed!")
                else:
                    f.close()

    if method == "lonlat":
        # Testing the lon-lat version of SPLASH:
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

        # Set the output directory for saving lon-lat gridded data files:
        output_dir = os.path.join(os.path.expanduser("~"), "Data", "splash")

        old_z = -1

        # Iterate through each latitude:
        for y in range(len(my_data.latitude)):
            lat = my_data.latitude[y]

            # Print progress:
            z = 100.0 * (y - 1.0) / (len(my_data.latitude) - 1.0)
            if int(z) % 5 == 0 and int(z) != old_z:
                msg = "[{}{}] {}% \r".format(
                    '#'*int(z / 5), '\u2013'*(20 - int(z / 5)), int(z))
                sys.stdout.write(msg)
                old_z = int(z)

            # Iterate through each longitude:
            for lon in my_data.longitude:
                my_data.get_annual_cru(2000, lat, lon)

                # Prepare the output file:
                if os.path.isdir(output_dir):
                    output_file = os.path.join(
                        output_dir,
                        "splash_cru_%0.2f_%0.2f.out" % (lat, lon))
                    try:
                        f = open(output_file, 'w')
                        tmp = SPLASH(0, 0)
                        f.write(tmp._header)
                    except:
                        raise
                    else:
                        f.close()
                        tmp = None

                my_class = SPLASH(lat, my_data.elv)
                my_class.lon = lon

                if my_data.is_okay:
                    my_class.spin_up(my_data)

                    # Loop through a year:
                    for i in range(my_data.npoints):
                        # Get preceding soil moisture status:
                        if i == 0:
                            wn = my_class.wn_vec[-1]
                        else:
                            wn = my_class.wn_vec[i-1]

                        # Calculate soil moisture and runoff:
                        my_class.run_one_day(n=i+1,
                                             y=my_data.year,
                                             wn=wn,
                                             sf=my_data.sf_vec[i],
                                             tc=my_data.tair_vec[i],
                                             pn=my_data.pn_vec[i])
                        my_class.wn_vec[i] = my_class.wn

                        # Save daily results to file:
                        if os.path.isfile(output_file):
                            try:
                                f = open(output_file, 'a')
                                f.write(my_class.data_str)
                            except:
                                logging.exception("Write failed!")
                            else:
                                f.close()
                else:
                    if os.path.isfile(output_file):
                        my_class.year = my_data.year
                        try:
                            f = open(output_file, 'a')
                            f.write(my_class.data_str)
                        except:
                            logging.exception("Write failed!")
                        else:
                            f.close()
        # Print the 100%
        z = 100.0
        msg = "[{}{}] {}%\n".format(
            '#'*int(z / 5), '\u2013'*(20 - int(z / 5)), int(z))
        sys.stdout.write(msg)
