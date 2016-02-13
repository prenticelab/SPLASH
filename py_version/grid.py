#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# grid.py
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
#
###############################################################################
# IMPORT MODULES:
###############################################################################
import datetime
import logging

from data_grid import DATA_G


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

    my_day = 172
    my_year = 2000
    my_date = datetime.date(my_year, 1, 1)
    my_date += datetime.timedelta(days=(my_day-1))

    data.read_monthly_clim(my_date)
    root_logger.info("read %d precipitation", data.pre.size)
    root_logger.info("read %d sunshine fraction", data.sf.size)
    root_logger.info("read %d air temperature", data.tair.size)
