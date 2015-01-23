#!/usr/bin/python
#
# stash_getdata.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2015-01-22 -- created
# 2015-01-23 -- last updated
#
# ------------
# description:
# ------------
# The script creates a CSV of the necessary daily meteorological variables for
# running the STASH code, including: fraction of bright sunshine, Sf (unitless), 
# air temperature, Tair (deg. C), and precipitation, Pn (mm).
#
# Supported input sources currently include:
#   Sf:
#    * CRU TS3.2 Cld (cloudiness fraction)
#   Tair:
#    * CRU TS3.2 Tmn (monthly mean daily air temperature)
#    * WATCH daily Tair
#   Pn
#    * CRU TS3.2 Pre (monthly total precipitation)
#    * WATCH daily Rainf
#
# -----
# todo:
# -----
#  * loop through days of interest
#  * write data to CSV file
#
###############################################################################
## IMPORT MODULES
###############################################################################
import datetime
import glob
import numpy
from scipy.io import netcdf

###############################################################################
## FUNCTIONS
###############################################################################
def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime date
    Output:   datetime date
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32) 
    dt3 = dt2.replace(day=1)
    return dt3

def density_h2o(tc, p):
    """
    Name:     density_h2o
    Input:    - float, air temperature (tc), degrees C
              - float, atmospheric pressure (p), Pa
    Output:   float, density of water, kg/m^3
    Features: Calculates density of water at a given temperature and 
              pressure
    Ref:      Chen et al. (1977)
    """
    # Calculate density at 1 atm:
    po = (
        0.99983952 + 
        (6.788260e-5)*tc + 
        -(9.08659e-6)*tc*tc +
        (1.022130e-7)*tc*tc*tc + 
        -(1.35439e-9)*tc*tc*tc*tc +
        (1.471150e-11)*tc*tc*tc*tc*tc +
        -(1.11663e-13)*tc*tc*tc*tc*tc*tc + 
        (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc + 
        -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
    )
    #
    # Calculate bulk modulus at 1 atm:
    ko = (
        19652.17 +
        148.1830*tc + 
        -2.29995*tc*tc + 
        0.01281*tc*tc*tc + 
        -(4.91564e-5)*tc*tc*tc*tc + 
        (1.035530e-7)*tc*tc*tc*tc*tc
    )
    #
    # Calculate temperature dependent coefficients:
    ca = (
        3.26138 + 
        (5.223e-4)*tc + 
        (1.324e-4)*tc*tc + 
        -(7.655e-7)*tc*tc*tc + 
        (8.584e-10)*tc*tc*tc*tc
    )
    cb = (
        (7.2061e-5) +
        -(5.8948e-6)*tc + 
        (8.69900e-8)*tc*tc + 
        -(1.0100e-9)*tc*tc*tc + 
        (4.3220e-12)*tc*tc*tc*tc
    )
    #
    # Convert atmospheric pressure to bar (1 bar = 100000 Pa)
    pbar = (1.0e-5)*p
    #
    pw = (
        (1e3)*po*(ko + ca*pbar + cb*pbar**2.0)/(
            ko + ca*pbar + cb*pbar**2.0 - pbar
        )
    )
    return pw

def get_daily_watch(d, ct, v):
    """
    Name:     get_daily_watch
    Input:    - str, directory to WATCH netcdf file (d)
              - datetime date, current datetime.date object (ct)
              - str, variable of interest (v)
    Output:   numpy nd.array
    Features: Returns 360x720 monthly WATCH dataset for a given date and 
              variable of interest (e.g., Tair, Rainf)
    """
    # Search directory for netCDF file:
    my_path = '%s%s*%d%02d.nc' % (d, v, ct.year, ct.month)
    my_file = glob.glob(my_path)[0]
    #
    if my_file:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")
        #
        #   VARIABLES:
        #    * day (int), 
        #    * lat (grid box center, degrees north)
        #    * lon (grid box center, degrees east)
        #    * timestp (days since beginning of month)
        #
        #   DATA:
        #    * 'Rainf' units = kg m^-2 s^-1
        #    * 'Tair' units = Kelvin  
        #    *  Missing value = 1.00e+20
        #
        # Find time index:
        f_time = f.variables['day'].data
        ti = numpy.where(ct.day == f_time)[0][0]
        #
        # Get the spatial data for current time:
        f_data = f.variables[v].data[ti]
        #
        f.close()
        return f_data

def get_month_days(ts):
    """
    Name:     get_month_days
    Input:    datetime date
    Output:   int
    Depends:  add_one_month
    Features: Returns the number of days in the month
    """
    ts1 = ts.replace(day=1)
    ts2 = add_one_month(ts1)
    dts = (ts2 - ts1).days
    return dts

def get_monthly_cru(d, ct, v):
    """
    Name:     get_monthly_cru
    Input:    - str, directory to CRU netcdf file (d)
              - datetime date, current month datetime object (ct)
              - str, variable of interest (v)
    Output:   numpy nd.array
    Depends:  get_time_index
    Features: Returns 360x720 monthly CRU TS dataset for a given month and 
              variable of interest (e.g., cld, pre, tmp)
    """
    # Search directory for netCDF file:
    my_file = glob.glob(d + "*" + v + ".dat.nc")[0]
    #
    if my_file:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")
        #
        # Save data for variables of interest:
        # NOTE: for CRU TS 3.21: 
        #       variables: 'lat', 'lon', 'time', v
        #       where v is 'tmp', 'pre', 'cld', 'vap'
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
        #       'vap' units = hPa
        #       Missing value = 9.96e+36
        # Save the base time stamp:
        bt = datetime.date(1900,1,1)
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

def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime date, base timestamp (bt)
              - datetime date, current timestamp to be found (ct)
              - numpy nd.array, days since base timestamp (aot)
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
    aot = numpy.append(aot, [dt,])
    #
    # Find the first index of dt in the sorted array:
    idx = numpy.where(numpy.sort(aot)==dt)[0][0]
    return idx

def get_xy_nc(lon, lat, r):
    """
    Name:     get_xy_nc
    Input:    - float, longitude (lon)
              - float, latitude (lat)
              - float, resolution (r)
    Output:   tuple, x-y indices
    Features: Returns array indices for a given lon-lat pair and pixel 
              resolution for a netCDF file (i.e., inverted raster)
    """
    # Solve x and y indices:
    x = (lon + 180.0)/r - 0.5
    y = (lat + 90.0)/r - 0.5
    #
    return (int(x), int(y))

def grid_centroid(my_lon, my_lat):
    """
    Name:     grid_centroid
    Input:    - float, longitude (my_lon)
              - float, latitude (my_lat)
    Output:   tuple, longitude latitude pair
    Features: Returns the nearest 0.5 deg. grid centroid per given coordinates
              based on the Euclidean distance to each of the four surrounding 
              grids; if any distances are equivalent, the pixel north and east
              is selected by default
    """
    # Create lists of regular latitude and longitude:
    grid_res = 0.5
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = 360
    lon_dim = 720
    lats = [lat_min + y*grid_res for y in xrange(lat_dim)]
    lons = [lon_min + x*grid_res for x in xrange(lon_dim)]
    #
    # Find bounding longitude:
    centroid_lon = None
    if my_lon in lons:
        centroid_lon = my_lon
    else:
        lons.append(my_lon)
        lons.sort()
        lon_index = lons.index(my_lon)
        bb_lon_min = lons[lon_index-1]
        try:
            bb_lon_max = lons[lon_index+1]
        except IndexError:
            bb_lon_max = lons[-1] + grid_res
        #
    # Find bounding latitude:
    centroid_lat = None
    if my_lat in lats:
        centroid_lat = my_lat
    else:
        lats.append(my_lat)
        lats.sort()
        lat_index = lats.index(my_lat)
        bb_lat_min = lats[lat_index-1]
        try:
            bb_lat_max = lats[lat_index+1]
        except IndexError:
            bb_lat_max = lats[-1] + grid_res
        #
    # Determine nearest centroid:
    # NOTE: if dist_A equals dist_B, then centroid defaults positively 
    #       i.e., north / east
    if centroid_lon and centroid_lat:
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lon and not centroid_lat:
        # Calculate the distances between lat and bounding box:
        dist_A = bb_lat_max - my_lat
        dist_B = my_lat - bb_lat_min
        if dist_A > dist_B:
            centroid_lat = bb_lat_min
        else:
            centroid_lat = bb_lat_max
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lat and not centroid_lon:
        # Calculate the distances between lon and bounding box:
        dist_A = bb_lon_max - my_lon
        dist_B = my_lon - bb_lon_min
        if dist_A > dist_B:
            centroid_lon = bb_lon_min
        else:
            centroid_lon = bb_lon_max
        my_centroid = (centroid_lon, centroid_lat)
    else:
        # Calculate distances between lat:lon and bounding box:
        # NOTE: if all distances are equal, defaults to NE grid
        dist_A = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_B = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        dist_C = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_D = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        min_dist = min([dist_A, dist_B, dist_C, dist_D])
        #
        # Determine centroid based on min distance:
        if dist_A == min_dist:
            my_centroid = (bb_lon_max, bb_lat_max)
        elif dist_B == min_dist:
            my_centroid = (bb_lon_max, bb_lat_min)
        elif dist_C == min_dist:
            my_centroid = (bb_lon_min, bb_lat_max)
        elif dist_D == min_dist:
            my_centroid = (bb_lon_min, bb_lat_min)
            #
    # Return nearest centroid:
    return my_centroid

###############################################################################
## USER VARIABLES
###############################################################################
mac = 0
if mac:
    # Directories on your Mac machine:
    cru_dir = {'cld' : '',
               'pre' : '',
               'tmn' : ''}
    watch_dir = {'Rainf' : '',
                 'Tair' : ''}
else:
    # Directories on you Linux machine:
    cru_dir = {'cld' : '/usr/local/share/database/cru/',
               'pre' : '/usr/local/share/database/cru/',
               'tmn' : '/usr/local/share/database/cru/'}
    watch_dir = {'Rainf' : '/usr/local/share/database/watch/netcdf/rainf/',
                 'Tair' : '/usr/local/share/database/watch/netcdf/tair/'}

# Define where you want to receive data from:
data_src = {'sf' : 'cru',
            'tair' : 'watch',
            'pn' : 'watch'}

# Longitude and latitude of your location of interest:
user_lon = -122.4
user_lat = 37.7

###############################################################################
## MAIN
###############################################################################
# CRU TS3.2 and WATCH WFDEI data are in 0.5x0.5 degree resolution.
##
## 1. Get 0.5 degree pixel associated with user coordinates: 
##
(px_lon, px_lat) = grid_centroid(user_lon, user_lat)
(px_x, px_y) = get_xy_nc(px_lon, px_lat, 0.5)

##
## 2. Set the date you want data for
##
cur_date = datetime.date(2000, 1, 4)

##
## 3. Get netCDF files for Sf, Tair, and Pn:
##
if data_src['sf'] == 'cru':
    # Look up CRU cloudiness data:
    voi = 'cld'
    d = get_monthly_cru(cru_dir[voi], cur_date, voi)
    my_sf = d[px_y, px_x]    # %
    my_sf /= 100.0           # unitless
    my_sf = 1.0 - my_sf      # complement of cloudiness

if data_src['tair'] == 'cru':
    # Look up CRU air temperature data:
    voi = 'tmn'
    d = get_monthly_cru(cru_dir[voi], cur_date, voi)
    my_tair = d[px_y, px_x]  # deg. C
elif data_src['tair'] == 'watch':
    # Look up WATCH air temperature data:
    voi = 'Tair'
    d = get_daily_watch(watch_dir[voi], cur_date, voi)
    my_tair = d[px_y, px_x]  # Kelvin
    my_tair -= 273.15        # deg. C

if data_src['pn'] == 'cru':
    # Get number of days in current month:
    nm = get_month_days(cur_date)
    #
    # Look up CRU precipitation data:
    voi = 'pre'
    d = get_monthly_cru(cru_dir[voi], cur_date, voi)
    my_pre = d[px_y, px_x]   # mm/month
    my_pre /= nm             # mm/day
elif data_src['pn'] == 'watch':
    # Calculate water density (or use standard value):
    if my_tair:
        pw = density_h2o(my_tair, 101325)
    else:
        pw = 1e3
    #
    # Look up WATCH rainfall data:
    voi = 'Rainf'
    d = get_daily_watch(watch_dir[voi], cur_date, voi)
    my_pre = d[px_y, px_x]   # kg m^-2 s^-1
    my_pre /= pw             # m/s
    my_pre *= 8.64e7         # mm/day
