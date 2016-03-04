#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# grid.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2016-02-19
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
###############################################################################
# IMPORT MODULES:
###############################################################################
import datetime
import logging

from data_grid import DATA_G
from splash import SPLASH


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("grid.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # Instantiate DATA_G class:
    cru_dir = "/usr/local/share/data/cru"
    data = DATA_G()
    data.find_cru_files(cru_dir)
    data.read_elv()
    data.read_lon_lat()

    day = 172
    year = 2000
    my_date = datetime.date(year, 1, 1)
    my_date += datetime.timedelta(days=(day-1))

    data.read_monthly_clim(my_date)
    root_logger.info("read %d precipitation", data.pre.size)
    root_logger.info("read %d sunshine fraction", data.sf.size)
    root_logger.info("read %d air temperature", data.tair.size)

    sm = 75
    for i in range(len(data.latitude)):
        for j in range(len(data.longitude)):
            lat = data.latitude[i]
            lon = data.longitude[j]
            elv = data.elevation[i, j]
            sf = data.sf[i, j]
            tc = data.tair[i, j]
            pn = data.pre[i, j]
            if elv != data.error_val:
                root_logger.info("lat: %f, lon: %f, elv: %f" % (lat, lon, elv))
                my_class = SPLASH(lat, elv)
                my_class.run_one_day(day, year, sm, sf, tc, pn)
