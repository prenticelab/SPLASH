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

import numpy

from const import pir


###############################################################################
# FUNCTIONS:
###############################################################################
def calculate_latitude(y, r):
    """
    Name:     calculate_latitude
    Input:    - nd.array, latitude index (y)
              - float, pixel resolution (r)
    Output:   nd.array, latitude, degrees
    Features: Returns latitude for an index array (numbered from the
              bottom-left corner) and pixel resolution
    """
    # Offset lat to pixel centroid and calucate based on index:
    lat = -90.0 + (0.5*r)
    lat += (y*r)
    if isinstance(y, int):
        logging.debug(
            "latitude at index %d is %f at resolution %0.3f" % (y, lat, r))
    else:
        logging.debug("calculating latitude at %0.3f degrees", r)
    return lat


def dcos(x):
    """
    Name:     dcos
    Input:    float/nd.array, angle, degrees (x)
    Output:   float/nd.array, cos(x*pi/180)
    Features: Calculates the cosine of an array of angles in degrees
    """
    if isinstance(x, float):
        logging.debug("calculating cosine of %f degrees", x)
    elif isinstance(x, numpy.ndarray):
        logging.debug("calculating cosine of numpy array of length %d", x.size)
    return numpy.cos(x*pir)


def dsin(x):
    """
    Name:     dsin
    Input:    float/nd.array, angle, degrees (x)
    Output:   float/nd.array, sin(x*pi/180)
    Features: Calculates the sine of an array of angles in degrees
    """
    if isinstance(x, float):
        logging.debug("calculating sine of %f degrees", x)
    elif isinstance(x, numpy.ndarray):
        logging.debug("calculating cosine of numpy array of length %d", x.size)
    return numpy.sin(x*pir)


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
