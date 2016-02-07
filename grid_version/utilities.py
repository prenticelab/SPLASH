#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# utilities.py
#
# LAST UPDATED: 2016-01-29
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


def get_x_y(lon, lat, r=0.5):
    """
    Name:     get_x_y
    Input:    - float, longitude, degrees (lon)
              - float, latitude, degrees (lat)
              - [optional] float, resolution, degrees (r)
    Output:   tuple, x-y indices
    Features: Returns x and y indices for a given lon-lat pair and pixel
              resolution. Y is latitude (0--359). X is longitude (0--719)
    """
    if lon > 179.75 or lon < -179.75:
        logging.warning("longitude outside range of validity")
        return (None, None)
    elif lat > 89.75 or lat < -89.75:
        logging.warning("latitude outside range of validity")
        return (None, None)
    else:
        # Solve x and y indices:
        x = (lon + 180.0)/r - 0.5
        y = (lat + 90.0)/r - 0.5
        return (int(x), int(y))
