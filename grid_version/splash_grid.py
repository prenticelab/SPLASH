#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# splash_grid.py
#
# 2014-01-30 -- created
# 2016-01-29 -- last updated
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
# led algorithms for simulating habitats (SPLASH): Robust indices of radiation,
# evapotranspiration and plant-available moisture, Geoscientific Model
# Development, 2015 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
#
# NOT WORKING
#
# This script calculates the monthly outputs at 0.5 deg resolution based on
# the SPLASH model
#
# IMPORTANT NOTE:
#   Global variables are defined inside a function at definition; therefore,
#   you must re-run function definitions if you change global variable values
#
# Key (monthly) outputs:
# 1. PPFD, mol/m^2
# 2. EET, mm
# 3. PET, mm
# 4. CP_alpha = SUM_mo(AETi) / SUM_mo(EETi)
#    where:
#     SUM_mo = sum over the days in a month (i)
#     AETi = daily actual evapotranspiration, mm
#     EETi = daily equilibrium evapotranspiration, mm
# 5. CWD = SUM_mo(PETi - AETi)
#    where:
#      PETi = daily potential evapotranspiration, mm
#
# The following input files are required:
# 1. Tc, CRU TS 3.21 Monthly Mean Daily Temperature, deg. C
#    "cru_ts3.21.1901.2012.tmp.dat.nc"
# 2. Pre, CRU TS 3.21 Monthly Precipitation, mm
#    "cru_ts3.21.1901.2012.pre.dat.nc"
# 3. Cld, CRU TS 3.21 Monthly Fractional Cloud Cover, %
#    "cru_ts3.21.1901.2012.cld.dat.nc"
# 4. z, CRU TS 3.00 Elevation Data, m
#    "halfdeg.elv.grid.dat"
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 00. created based on cramer_prentice.py [14.01.30]
# 01. file handling for cru ts 3.00 ('elv') [14.01.30]
# 02. file handling for cru ts 3.21 ('pre', 'cld', and 'tmp') [14.01.31]
# 03. file handling for wfdei ('SWdown') [14.01.31]
# 04. avoiding overflow error, set_missing on cru data [14.01.31]
# 05. reverse column ordering for raster output [14.01.31]
# 06. fixed precip in omega update [14.01.31]
# 07. parsed out functions from calc_demand [14.01.31]
# 08. created clip_field and error_field [14.02.01]
# 09. added conversion from mm/h to mm/d [14.02.01]
# 10. fixed the update for omega (daily aet, not commulative) [14.02.01]
# 11. fixed demand equation (1000/rho_water) [14.02.02]
# 12. started ALPHA class [14.02.17]
# --> removed ALPHA class (incomplete) [14.03.30]
# 13. implementing STASH 2.0 features [14.08.28]
# 14. EVAP_G class (based on EVAP from stash.py) [14.09.02]
# 15. Moved land_clip and error_field outside loop [14.09.04]
# 16. created kerror global constant [14.09.04]
# 17. created save_to_file using kerror [14.09.04]
# 18. corrected equation for CP_alpha (AET/EET) [14.09.04]
# 19. created spin up for soil moisture [14.09.05]
# 20. updated extrapolation equation [14.09.24]
# 21. made my_lambda & my_new public variables [14.09.24]
# 22. added out_no and an_dic for easier processing [14.09.24]
# 23. updated komega to float --- may influence cooper delta [14.09.29]
# --> no difference in results
# 24. added 'berger' lamm method [14.09.30]
# 25. added check to woolf's method for lambda > 360 [14.10.01]
# 26. added spencer method for delta [14.10.10]
# 27. replaced tau with Allen (1996); removed kZ [14.10.10]
# 28. distinguish shortwave from visible light albedo [14.10.16]
# 29. added some if 0 statements for annual PPFD [14.10.21]
# 30. created mean_monthly_w function [14.10.23]
# 31. updated value and reference for semi-major axis, a [14.10.31]
# 32. fixed Cooper's and Spencer's declination equations [14.11.25]
# 33. replaced simplified kepler with full kepler [14.11.25]
# 34. removed options for approximation methods not considering variable
#     orbital velocity (e.g. Spencer, Woolf, Klein, Cooper, and Circle
#     methods) [15.01.13]
# 35. updated get_monthly_cru to set nodata values equal to kerror [15.02.25]
# 36. updated get_elevation to set nodata values equal to kerror [15.02.25]
# 37. updated EVAP_G according to stash.py [15.02.25]
#     --> removed loops for hs, hn, and hi
# 38. updated econ variable functions with NODATA preservation [15.02.25]
# 39. created DATA_G & STASH_G classes [15.02.25]
# 40. continued work on DATA_G, STASH_G & main program [15.02.26]
#     --> added nc_lat, nc_lon, nc_history, nc_write, get_date functions
# 41. renaming (addresses issue #3) [15.08.23]
# 42. import global constants from const.py [15.08.23]
# 43. PEP8 style fixes [15.11.11]
#
# ~~~~~
# todo:
# ~~~~~
# ! Monthly AET divided by monthly EET goes greater than 1.26
# ! Daily AET/EET as high as 1.67 (Namibia, Africa)
#
# 0. Go through each part of EVAP_G to make sure it is working properly
# 1. complete DATA_G class
#    --> read temperature, precip and cloudiness data for given month
#    --> don't forget to scale precip!
# 2. complete STASH_G class
#    --> update the outline from STASH class
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import numpy
from scipy.io import netcdf

from const import (kA, kalb_sw, kalb_vis, kb, kc, kCw, kd, ke, keps, kfFEC,
                   kG, kGsc, kL, kMa, kMv, kPo, kR, kTo, kWm, kw, komega)
from data_grid import DATA_G

###############################################################################
## GLOBAL CONSTANTS:
###############################################################################
kerror = -9999.  # error value


###############################################################################
## CLASSES:
###############################################################################
class EVAP_G:
    """
    Name:     EVAP_G
    Features: This class calculates daily 360x720 gridded radiation and
              evapotranspiration quantities
              - daily PPFD (ppfd_d), mol/m^2
              - daily EET (eet_d), mm
              - daily PET (pet_d), mm
              - daily AET (aet_d), mm
              - daily condensation (wc), mm
    Refs:     Allen, R.G. (1996), Assessing integrity of weather data for
                reference evapotranspiration estimation, Journal of Irrigation
                and Drainage Engineering, vol. 122, pp. 97--106.
              Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998),
                'Meteorological data,' Crop evapotranspiration - Guidelines for
                computing crop water requirements - FAO Irrigation and drainage
                paper 56, Food and Agriculture Organization of the United
                Nations, online: http://www.fao.org/docrep/x0490e/x0490e07.htm
              Berger, A.L. (1978), Long-term variations of daily insolation and
                quarternary climatic changes, Journal of Atmospheric Sciences,
                vol. 35, pp. 2362--2367.
              Berger, A.L., M.F. Loutre, and C. Tricot (1993), Insolation and
                Earth's orbital periods, J. Geophys. Res., 98, 10341--10362.
              Duffie, J. A. and W. A. Beckman (1991). Solar engineering of
                thermal processes. 4th ed. New Jersey: John Wiley and Sons
              Federer (1982), Transpirational supply and demand: plant, soil,
                and atmospheric effects evaluated by simulation, Water
                Resources Research, vol. 18, no. 2, pp. 355--362.
              Ge, S., R.G. Smith, C.P. Jacovides, M.G. Kramer, R.I. Carruthers
                (2011), Dynamics of photosynthetic photon flux density (PPFD)
                and estimates in coastal northern California, Theoretical and
                Applied Climatology, vol. 105, pp. 107--118.
              Henderson-Sellers, B. (1984), A new formula for latent heat of
                vaporization of water as a function of temperature, Quarterly
                Journal of the Royal Meteorological Society 110, pp. 1186–1190
              Linacre (1968), Estimating the net-radiation flux, Agricultural
                Meteorology, vol. 5, pp. 49--63.
              Prentice, I.C., M.T. Sykes, W. Cramer (1993), A simulation model
                for the transient effects of climate change on forest
                landscapes, Ecological Modelling, vol. 65, pp. 51--70.
              Priestley, C.H.B. and R.J. Taylor (1972), On the assessment of
                surface heat flux and evaporation using large-scale parameters,
                Monthly Weather Review, vol. 100 (2), pp. 81--92.
              Spencer, J. W. (1971), Fourier series representation of the
                position of the sun, Search, vol. 2, p. 172.
              Stine, W. B. and M. Geyer (2001). “Power from the Sun”.
                online: http://www.powerfromthesun.net/Book/chapter03/chapter03
              Wetherald, R.T., S. Manabe (1972), Response to joint ocean-
                atmosphere model to the seasonal variation of the solar
                radiation, Monthly Weather Review, vol. 100 (1), pp. 42--59.
              Woolf, H. M. (1968). On the computation of solar evaluation
                angles and the determination of sunrise and sunset times.
                Tech. rep. NASA-TM-X-164. National Aeronautics and Space
                Administration (NASA).
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, n, elv, sf, tc, sw, y=0):
        """
        Name:     EVAP_G.__init__
        Input:    - int, day of the year (n)
                  - numpy nd.array, elevation, m (elv)
                  - numpy nd.array, fraction of sunshine hours (sf)
                  - numpy nd.array, mean daily air temperature, deg C (tc)
                  - numpy nd.array, evaporative supply rate, mm/hr (sw)
                  - int, year (y)
        """
        # Assign default public variables:
        self.user_year = y
        #self.user_elv = elv
        #self.user_sf = sf
        #self.user_tc = tc
        #self.user_sw = sw
        #
        # Error handle the day of the year:
        if n < 1 or n > 366:
            print "Day of year outside range of validity (1 to 366)!"
            exit(1)
        else:
            self.user_day = n
            #
        #
        # Create lon and lat arrays (degrees):
        my_x = numpy.array([j for j in xrange(720)])
        my_y = numpy.array([j for j in xrange(360)])
        (lon_array, lat_array) = self.get_lon_lat(my_x, my_y, 0.5)
        #
        # Convert lat array to grids (degrees)
        lat_grid = numpy.reshape(numpy.repeat(lat_array, 720), (360, 720), 'C')
        #
        # Define missing data indexes:
        #   sw doesn't have error values
        nodata_idx = numpy.where(
            (elv == kerror) | (tc == kerror) | (sf == kerror)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year, days
        #    kN, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            kN = 365
        else:
            kN = self.julian_day((y+1), 1, 1) - self.julian_day(y, 1, 1)
        self.kN = kN
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes, degrees
        #    my_nu, SCALAR
        #    my_lambda, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor, unitless
        #    dr, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        kee = ke**2
        my_rho = (1.0 - kee)/(1.0 + ke*self.dcos(my_nu))
        dr = (1.0/my_rho)**2
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle, degrees
        #    delta, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        pir = (numpy.pi/180.0)
        delta = numpy.arcsin(self.dsin(my_lambda)*self.dsin(keps))
        delta /= pir
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes, unitless
        #    ru, MATRIX (360x720)
        #    rv, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(delta)*self.dsin(lat_grid)
        rv = self.dcos(delta)*self.dcos(lat_grid)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the sunset hour angle, degrees
        #    hs, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed hour angle grid
        hs = numpy.zeros(shape=(360, 720))
        #
        # Indexes of pixels under polar night conditions:
        # Note: hs == 0 degrees (no further comp's)
        # hs_neg = numpy.where(ru/rv <= -1.0)
        #
        # Indexes of pixels under polar day conditions:
        # Note: hs == 180 degrees
        hs_pos = numpy.where(ru/rv >= 1.0)
        hs[hs_pos] += 180.0
        #
        # Indexes of pixels for regular conditions:
        # Note: hs = acos(-u/v)
        hs_reg = numpy.where((ru/rv < 1.0) & (ru/rv > -1.0))
        hs[hs_reg] -= 1.0
        hs[hs_reg] *= ru[hs_reg]
        hs[hs_reg] /= rv[hs_reg]
        hs[hs_reg] = numpy.arccos(hs[hs_reg])
        hs[hs_reg] /= pir
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate daily extraterrestrial solar radiation, J/m^2
        #    ra_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 40, Documentation
        # Note: ru = sin(delta)*sin(phi); rv = cos(delta)*cos(phi)
        ra_d = (86400.0/numpy.pi)*kGsc*dr*((ru*hs)*pir + (rv*self.dsin(hs)))
        self.ra_d = ra_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate transmittivity, unitless
        #    tau, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
        tau_o = kc + (kd*sf)
        tau = tau_o*(1.0 + (2.67e-5)*elv)
        #
        # Note: missing == kerror
        tau[nodata_idx] *= 0.0
        tau[nodata_idx] += kerror
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        #    ppfd_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd_d = ra_d
        ppfd_d *= tau
        ppfd_d *= (1e-6)*(kfFEC)*(1.0 - kalb_vis)
        #
        # Maintain error values:
        ppfd_d[nodata_idx] *= 0.0
        ppfd_d[nodata_idx] += kerror
        self.ppfd_d = ppfd_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation, W/m^2
        #     rnl, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        rnl = (kb + (1.0 - kb)*sf)*(kA - tc)
        #
        # Maintain error values:
        rnl[nodata_idx] *= 0.0
        rnl[nodata_idx] += kerror
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate variable substitute, W/m^2
        #     rw, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = tau
        rw *= (1.0 - kalb_sw)*kGsc*dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate net radiation cross-over hour angle, degrees
        #     hn, MATRIX (360,720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed cross-over angle:
        hn = numpy.zeros(shape=(360, 720))
        #
        # Define cosine of hn:
        cos_hn = numpy.copy(rnl)
        cos_hn -= (rw*ru)
        cos_hn /= rw
        cos_hn /= rv
        #
        # Indexes of pixels where Rnl is all-day positive:
        hn_pos = numpy.where(cos_hn <= -1.0)
        hn[hn_pos] += 180.0
        #
        # Indexes of pixels for regular Rnl:
        hn_reg = numpy.where((cos_hn < 1.0) & (cos_hn > -1.0))
        hn[hn_reg] = numpy.arccos(cos_hn[hn_reg])
        hn[hn_reg] /= pir
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate daytime net radiation, J/m^2
        #     rn_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rn_d = (86400.0/numpy.pi)*((rw*ru - rnl)*hn*pir + rw*rv*self.dsin(hn))
        #
        # Maintain error values:
        rn_d[nodata_idx] *= 0.0
        rn_d[nodata_idx] += kerror
        self.rn_d = rn_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate nighttime net radiation, J/m^2
        #     rnn_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rnn_d = rw*ru*(hs - hn)*pir
        rnn_d += rw*rv*(self.dsin(hs) - self.dsin(hn))
        rnn_d += rnl*(numpy.pi - 2.0*hs*pir + hn*pir)
        rnn_d *= (86400.0/numpy.pi)
        #
        # Maintain error values:
        rnn_d[nodata_idx] *= 0.0
        rnn_d[nodata_idx] += kerror
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 15. Calculate water-to-energy conversion, m^3/J
        #     econ, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        s = self.sat_slope(tc)                         # Pa/K
        lv = self.enthalpy_vap(tc)                     # J/kg
        pw = self.density_h2o(tc, self.elv2pres(elv))  # kg/m^3
        g = self.psychro(tc, self.elv2pres(elv))       # Pa/K
        econ = s/(lv*pw*(s + g))
        #
        # Maintain error values:
        econ[nodata_idx] *= 0.0
        econ[nodata_idx] += kerror
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 16. Calculate daily condensation, mm
        #     cn, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cn = (1e3)*econ
        cn *= numpy.abs(rnn_d)
        #
        # Maintain error values:
        cn[nodata_idx] *= 0.0
        cn[nodata_idx] += kerror
        self.cond = cn
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 17. Estimate daily EET, mm
        #     eet_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eet_d = (1e3)*econ
        eet_d *= rn_d
        #
        # Maintain error values:
        eet_d[nodata_idx] *= 0.0
        eet_d[nodata_idx] += kerror
        self.eet_d = eet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 18. Estimate daily PET, mm
        #     pet_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pet_d = (1.0 + kw)*eet_d
        #
        # Maintain error values:
        pet_d[nodata_idx] *= 0.0
        pet_d[nodata_idx] += kerror
        self.pet_d = pet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 19. Calculate variable substitute, (mm/hr)/(W/m^2)
        #     rx, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = (3.6e6)*(1.0 + kw)*econ
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 20. Calculate the intersection hour angle, degrees
        #     hi, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed intersection hour angle grid
        hi = numpy.zeros(shape=(360, 720))
        #
        # Compute cos(hi)
        cos_hi = sw
        cos_hi /= (rw*(rv*rx))
        cos_hi += rnl/(rw*rv)
        cos_hi -= ru/rv
        #
        # Indexes where supply limits demand everywhere
        hi_pos = numpy.where(cos_hi <= -1.0)
        hi[hi_pos] += 180.0
        #
        # Indexes for regular supply
        hi_reg = numpy.where((cos_hi < 1.0) & (cos_hi > -1.0))
        hi[hi_reg] = numpy.arccos(cos_hi[hi_reg])
        hi[hi_reg] /= pir
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 21. Estimate daily AET (aet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        aet_d = sw*hi*pir
        aet_d += rx*rw*rv*(self.dsin(hn) - self.dsin(hi))
        aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*pir
        aet_d *= (24.0/numpy.pi)
        #
        # Maintain error values:
        aet_d[nodata_idx] *= 0.0
        aet_d[nodata_idx] += kerror
        self.aet_d = aet_d

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def dcos(self, x):
        """
        Name:     EVAP_G.dcos
        Input:    float/nd.array, angle, degrees (x)
        Output:   float/nd.array, cos(x*pi/180)
        Features: Calculates the cosine of angle(s) given in degrees
        """
        return numpy.cos(x*numpy.pi/180.0)

    def dsin(self, x):
        """
        Name:     EVAP_G.dsin
        Input:    float/nd.array, angle, degrees (x)
        Output:   float/nd.array, sin(x*pi/180)
        Features: Calculates the sine of angle(s) given in degrees
        """
        return numpy.sin(x*numpy.pi/180.0)

    def berger_tls(self, n):
        """
        Name:     EVAP_G.berger_tls
        Input:    int, day of year
        Output:   tuple,
                  - true anomaly, degrees
                  - true longitude, degrees
        Features: Returns true anomaly and true longitude for a given day
        Depends:  - ke
                  - komega
                  - dsin
        Ref:      Berger, A. L. (1978), Long term variations of daily
                  insolation and quaternary climatic changes, J. Atmos. Sci.,
                  35, 2362-2367.
        """
        # Variable substitutes:
        xee = ke**2
        xec = ke**3
        xse = numpy.sqrt(1.0 - xee)
        pir = (numpy.pi/180.0)
        #
        # Mean longitude for vernal equinox:
        xlam = (ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(komega)
        xlam -= xee/4.0*(0.5 + xse)*self.dsin(2.0*komega)
        xlam += xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*komega)
        xlam *= 2.0
        xlam /= pir
        #
        # Mean longitude for day of year:
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        #
        # Mean anomaly:
        anm = (dlamm - komega)
        ranm = (anm*pir)
        #
        # True anomaly:
        ranv = ranm
        ranv += (2.0*ke - xec/4.0)*numpy.sin(ranm)
        ranv += 5.0/4.0*xee*numpy.sin(2.0*ranm)
        ranv += 13.0/12.0*xec*numpy.sin(3.0*ranm)
        anv = ranv/pir
        #
        # True longitude:
        my_tls = anv + komega
        if my_tls < 0:
            my_tls += 360.0
        elif my_tls > 360:
            my_tls -= 360.0
        #
        # True anomaly:
        my_nu = (my_tls - komega)
        if my_nu < 0:
            my_nu += 360.0
        #
        return(my_nu, my_tls)

    def get_lon_lat(self, x, y, r):
        """
        Name:     EVAP_G.get_lon_lat
        Input:    - int/nd.array, longitude index (x)
                  - int/nd.array, latitude index (y)
                  - float, pixel resolution (r)
        Output:   float/nd.array tuple, longitude(s) and latitude(s), degrees
        Features: Returns lon-lat pair for an x-y index pair (numbered from the
                  bottom-left corner) and pixel resolution
        """
        # Offset lat, lon to pixel centroid
        lon = -180.0 + (0.5*r)
        lat = -90.0 + (0.5*r)
        #
        # Offset lat, lon based on pixel index
        lon = lon + (x*r)
        lat = lat + (y*r)
        #
        return (lon, lat)

    def julian_day(self, y, m, i):
        """
        Name:     EVAP_G.julian_day
        Input:    - int, year (y)
                  - int, month (m)
                  - int, day of month (i)
        Output:   float, Julian Ephemeris Day
        Features: Converts Gregorian date (year, month, day) to Julian
                  Ephemeris Day
        Ref:      Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," Astronomical
                  Algorithms
        """
        if m <= 2.0:
            y -= 1.0
            m += 12.0
        #
        a = int(y/100)
        b = 2 - a + int(a/4)
        #
        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
        return jde

    def sat_slope(self, tc):
        """
        Name:     EVAP_G.sat_slope
        Input:    numpy nd.array, air temperatures (tc), degrees C
        Output:   numpy nd.array, slopes of sat vap press temp curve (s)
        Features: Calculates the slopes of the sat pressure temp curve, Pa/K
                  NODATA = kerror
        Ref:      Eq. 13, Allen et al. (1998)
        """
        # Get no data indexes:
        nodata_idx = numpy.where(tc == kerror)
        #
        s = 17.269*tc
        s /= (tc + 237.3)
        s = numpy.exp(s)
        s /= numpy.power((tc + 237.3), 2.0)
        s *= (17.269)*(237.3)*(610.78)
        #
        # Maintain error values:
        s[nodata_idx] *= 0.0
        s[nodata_idx] += kerror
        #
        return s

    def enthalpy_vap(self, tc):
        """
        Name:     EVAP_G.enthalpy_vap
        Input:    numpy nd.array, air temperatures (tc), degrees C
        Output:   numpy nd.array, latent heats of vaporization
        Features: Calculates the enthalpy of vaporization, J/kg
                  NODATA = kerror
        Ref:      Eq. 8, Henderson-Sellers (1984)
        """
        # Get no data indexes:
        nodata_idx = numpy.where(tc == kerror)
        #
        lv = (tc + 273.15)
        lv /= (tc + 273.15 - 33.91)
        lv = numpy.power(lv, 2.0)
        lv *= (1.91846e6)
        #
        # Maintain error values:
        lv[nodata_idx] *= 0.0
        lv[nodata_idx] += kerror
        #
        return lv

    def elv2pres(self, z):
        """
        Name:     EVAP_G.elv2pres
        Input:    numpy nd.array, elevations above sea level (z), m
        Output:   numpy nd.array, atmospheric pressures, Pa
        Features: Calculates atm. pressure for a given elevation
                  NODATA = kerror
        Depends:  Global constants:
                  - kPo     - kTo       - kL
                  - kMa     - kG        - kR
        Ref:      Allen et al. (1998)
        """
        # Find nodata points:
        nodata_idx = numpy.where(z == kerror)
        #
        p_exp = (kG*kMa/(kR*kL))
        p = -1.0*kL*z
        p /= kTo
        p += 1.0
        p = numpy.power(p, p_exp)
        p *= kPo
        #
        # Maintain error values:
        p[nodata_idx] *= 0.0
        p[nodata_idx] += kerror
        #
        return p

    def density_h2o(self, tair, p):
        """
        Name:     EVAP_G.density_h2o
        Input:    - numpy nd.array, air temperatures (tair), degrees C
                  - numpy nd.array, atmospheric pressures (p), Pa
        Output:   numpy nd.array, densities of water, kg/m^3
        Features: Calculates water density at a given temperature and pressure
                  NODATA = kerror
        Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
                  pure water and sea water, Tech. Rept., Marine Physical
                  Laboratory, San Diego, CA.
        """
        # Find nodata indexes:
        nodata_idx = numpy.where((tair == kerror) | (p == kerror))
        #
        # For computational efficiency, set error values equal to one:
        tc = numpy.copy(tair)
        tc[nodata_idx] *= 0.0
        #
        # Calculate lambda, (bar cm^3)/g:
        my_lambda = 21.55053*tc
        my_lambda += -0.4695911*tc*tc
        my_lambda += (3.096363e-3)*tc*tc*tc
        my_lambda += -(7.341182e-6)*tc*tc*tc*tc
        my_lambda += 1788.316
        my_lambda[nodata_idx] *= 0.0
        my_lambda[nodata_idx] += 1.0
        #
        # Calculate po, bar
        po = 58.05267*tc
        po += -1.1253317*tc*tc
        po += (6.6123869e-3)*tc*tc*tc
        po += -(1.4661625e-5)*tc*tc*tc*tc
        po += 5918.499
        po[nodata_idx] *= 0.0
        po[nodata_idx] += 0.5
        #
        # Calculate vinf, cm^3/g
        vinf = -(7.435626e-4)*tc
        vinf += (3.704258e-5)*tc*tc
        vinf += -(6.315724e-7)*tc*tc*tc
        vinf += (9.829576e-9)*tc*tc*tc*tc
        vinf += -(1.197269e-10)*tc*tc*tc*tc*tc
        vinf += (1.005461e-12)*tc*tc*tc*tc*tc*tc
        vinf += -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc
        vinf += (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc
        vinf += -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc
        vinf += 0.6980547
        vinf[nodata_idx] *= 0.0
        #
        # Convert pressure to bars (1 bar = 100000 Pa)
        pbar = (1e-5)*p
        pbar[nodata_idx] *= 0.0
        pbar[nodata_idx] += 0.5
        #
        # Calculate the specific volume (cm^3 g^-1):
        v = my_lambda
        v /= (po + pbar)
        v += vinf
        v[nodata_idx] *= 0.0
        v[nodata_idx] += 1.0
        #
        # Convert to density (g cm^-3) -> 1000 g/kg;
        # 1000000 cm^3/m^3 -> kg/m^3:
        rho = numpy.power(v, -1.)
        rho *= (1e3)
        #
        # Maintain error values:
        rho[nodata_idx] *= 0.0
        rho[nodata_idx] += kerror
        #
        return rho

    def psychro(self, tair, p):
        """
        Name:     EVAP_G.psychro
        Input:    - numpy nd.array, air temperatures (tair), degrees C
                  - numpy nd.array, atm. pressures (p), Pa
        Output:   numpy nd.array, psychrometric constant, Pa/K
        Features: Calculates the psychrometric constant for a given temperature
                  and pressure
                  NODATA = kerror
        Depends:  - kMa
                  - kMv
                  - enthalpy_vap
        Refs:     Allen et al. (1998); Tsilingiris (2008)
        """
        # Define nodata indexes:
        nodata_idx = numpy.where((tair == kerror) | (p == kerror))
        #
        # For compuational efficiency, set no data values to zero:
        tc = numpy.copy(tair)
        tc[nodata_idx] *= 0.0
        #
        # Calculate the specific heat capacity of water, J/kg/K
        # Eq. 47, Tsilingiris (2008)
        cp = (2.050632750e-3)*tc
        cp += -(1.631537093e-4)*tc*tc
        cp += (6.212300300e-6)*tc*tc*tc
        cp += -(8.830478888e-8)*tc*tc*tc*tc
        cp += (5.071307038e-10)*tc*tc*tc*tc*tc
        cp += 1.0045714270
        cp *= (1e3)
        cp[nodata_idx] *= 0.0
        #
        # Calculate latent heat of vaporization, J/kg
        lv = self.enthalpy_vap(tc)
        #
        # Multiply by the molecular weight of water vapor:
        lv *= kMv
        lv[nodata_idx] *= 0.0
        lv[nodata_idx] += 1.0
        #
        # Calculate psychrometric constant, Pa/K
        # Eq. 8, Allen et al. (1998)
        pc = kMa*p*cp
        pc /= lv
        #
        # Maintain error values:
        pc[nodata_idx] *= 0.0
        pc[nodata_idx] += kerror
        #
        return pc


class SPLASH_G:
    """
    Name:     SPLASH_G
    Features: This class updates daily gridded quantities of radiation,
              evapotranspiration, soil moisture and runoff based on the
              SPLASH methodology.

    @TODO: finish class
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, elv):
        """
        Name:     SPLASH.__init__
        Input:    numpy.ndarray, elevation, meters (elv)
        """
        # Assign public variables:
        self.elv = elv
        #
        # Initialize daily status variables:
        self.ho = numpy.zeros(shape=(360, 720))      # solar irradiation, J/m2
        self.hn = numpy.zeros(shape=(360, 720))      # net radiation, J/m2
        self.ppfd = numpy.zeros(shape=(360, 720))    # PPFD, mol/m2
        self.cond = numpy.zeros(shape=(360, 720))    # condensation water, mm
        self.wn = numpy.zeros(shape=(360, 720))      # soil moisture, mm
        self.precip = numpy.zeros(shape=(360, 720))  # precipitation, mm
        self.ro = numpy.zeros(shape=(360, 720))      # runoff, mm
        self.eet = numpy.zeros(shape=(360, 720))     # equilibrium ET, mm
        self.pet = numpy.zeros(shape=(360, 720))     # potential ET, mm
        self.aet = numpy.zeros(shape=(360, 720))     # actual ET, mm

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def add_one_day(self, dt0):
        """
        Name:     SPLASH_G.add_one_day
        Input:    datetime.date (dt0)
        Output:   datetime.date (dt1)
        Features: Adds one day to datetime
        """
        dt1 = dt0 + datetime.timedelta(days=1)
        #
        return dt1

    def spin_up(self, d, y):
        """
        Name:     SPLASH_G.spin
        Input:    - DATA class, (d)
                  - int, year (y)
        Output:   None.
        Features: Spins up the daily soil moisture.
        Depends:  - add_one_day
                  - quick_run
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Initialize dates & load data
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sd = datetime.date(y, 1, 1)
        ed = datetime.date(y+1, 1, 1)
        d.read_monthly_clim(sd)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Create a soil moisture array:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        wn_vec = numpy.zeros(shape=(d.ny, 360, 720))
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Run one year:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cur_date = sd
        while cur_date < ed:
            # Get preceding soil moisture status:
            i = (cur_date - sd).days
            if i == 0:
                wn = wn_vec[-1]
            else:
                wn = wn_vec[i-1]
            #
            # Update monthly climatology (if necessary):
            d.read_monthly_clim(cur_date)
            #
            # Calculate soil moisture and runoff:
            sm, ro = self.quick_run(n=i+1,
                                    y=y,
                                    wn=wn,
                                    sf=d.sf,
                                    tc=d.tair,
                                    pn=d.pre)
            #
            # Set no data values equal to zero & save to array:
            sm[d.noval_idx] *= 0.0
            wn_vec[i] = sm
            #
            # Increment date by one day:
            cur_date = self.add_one_day(cur_date)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate change in starting soil moisture:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        start_sm = wn_vec[0]
        d.read_monthly_clim(sd)
        end_sm, ro = self.quick_run(n=1,
                                    y=y,
                                    wn=wn_vec[-1],
                                    sf=d.sf,
                                    tc=d.tair,
                                    pn=d.pre)
        #
        # Set no data values equal to zero & calc difference:
        end_sm[d.noval_idx] *= 0.0
        diff_sm = numpy.abs(end_sm - start_sm)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Equilibrate:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        spin_count = 1
        while diff_sm.mean() > 1.0:
            cur_date = sd
            while cur_date < ed:
                # Get preceding soil moisture status:
                i = (cur_date - sd).days
                if i == 0:
                    wn = wn_vec[-1]
                else:
                    wn = wn_vec[i-1]
                #
                # Update monthly climatology (if necessary):
                d.read_monthly_clim(cur_date)
                #
                # Calculate soil moisture and runoff:
                sm, ro = self.quick_run(n=i+1,
                                        y=y,
                                        wn=wn,
                                        sf=d.sf,
                                        tc=d.tair,
                                        pn=d.pre)
                #
                # Set no data values equal to zero & save to array:
                sm[d.noval_idx] *= 0.0
                wn_vec[i] = sm
                #
                # Increment date by one day:
                cur_date = self.add_one_day(cur_date)
            #
            start_sm = wn_vec[0]
            d.read_monthly_clim(sd)
            end_sm, ro = self.quick_run(n=1,
                                        y=y,
                                        wn=wn_vec[-1],
                                        sf=d.sf,
                                        tc=d.tair,
                                        pn=d.pre)
            #
            # Set no data values equal to zero & calc difference:
            end_sm[d.noval_idx] *= 0.0
            diff_sm = numpy.abs(end_sm - start_sm)
            spin_count += 1
        #
        print "Spun", spin_count, "years"
        self.wn_vec = wn_vec
        self.wn = wn_vec[-1]
        #
        # Preserve error values:
        self.wn[d.noval_idx] *= 0.0
        self.wn[d.noval_idx] += kerror
        #

    def quick_run(self, n, y, wn, sf, tc, pn):
        """
        Name:     SPLASH_G.quick_run
        Inputs:   - int, day of year (n)
                  - int, year (y)
                  - numpy.ndarray, daily soil water content, mm (wn)
                  - numpy.ndarray, daily fraction of bright sunshine (sf)
                  - numpy.ndarray, daily air temperature, deg C (tc)
                  - numpy.ndarray, daily precipitation, mm (pn)
        Output:   - numpy.ndarray, soil moisture, mm (sm)
                  - numpy.ndarray, runoff, mm (ro)
        Features: Returns gridded daily soil moisture and runoff.
        Depends:  - kCw
                  - kWm
                  - EVAP_G
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate evaporative supply rate, mm/h
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sw = (kCw/kWm)*wn
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate radiation and evaporation quantities
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        my_evap = EVAP_G(n, self.elv, sf, tc, sw, y)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate today's soil moisture, mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sm = wn + pn + my_evap.cond - my_evap.aet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate runoff, mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ro = numpy.zeros(shape=(360, 720))
        #
        # Where bucket is too full, allocate to runoff and reset to max:
        full_idx = numpy.where(sm > kWm)
        ro[full_idx] = sm[full_idx] - kWm
        sm[full_idx] *= 0.0
        sm[full_idx] += kWm
        #
        # Where bucket is too empty, set soil moisture to min:
        empty_idx = numpy.where(sm < 0)
        sm[empty_idx] *= 0.0
        #
        return(sm, ro)

    def run_one_day(self, n, y, wn, sf, tc, pn):
        """
        Name:     SPLASH_G.run_one_day
        Inputs:   - int, day of year (n)
                  - int, year (y)
                  - numpy.ndarray, soil water content, mm (wn)
                  - numpy.ndarray, fraction of bright sunshine (sf)
                  - numpy.ndarray, air temperature, deg C (tc)
                  - numpy.ndarray, precipitation, mm (pn)
        Outputs:  None
        Features: Runs SPLASH model for one day.
                  model.
        Depends:  - kCw
                  - kWm
                  - EVAP_G
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Set meteorological variables:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.precip = pn    # daily precipitation, mm
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate evaporative supply rate (sw), mm/h
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sw = (kCw/kWm)*wn
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate radiation and evaporation quantities
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        my_evap = EVAP_G(n, self.elv, sf, tc, sw, y)
        self.ho = my_evap.ra_d      # daily solar irradiation, J/m2
        self.hn = my_evap.rn_d      # daily net radiation, J/m2
        self.ppfd = my_evap.ppfd_d  # daily PPFD, mol/m2
        self.cond = my_evap.cond    # daily condensation water, mm
        self.eet = my_evap.eet_d    # daily equilibrium ET, mm
        self.pet = my_evap.pet_d    # daily potential ET, mm
        self.aet = my_evap.aet_d    # daily actual ET, mm
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate today's soil moisture (sm), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        sm = wn + pn + my_evap.cond - my_evap.aet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate runoff (ro), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ro = numpy.zeros(shape=(360, 720))
        #
        # Where bucket is too full, allocate to runoff and reset to max:
        full_idx = numpy.where(sm > kWm)
        ro[full_idx] = sm[full_idx] - kWm
        sm[full_idx] *= 0.0
        sm[full_idx] += kWm
        #
        # Where bucket is too empty, reduce AET by discrepancy amount
        # & set soil moisture to min:
        empty_idx = numpy.where(sm < 0)
        self.aet[empty_idx] += sm[empty_idx]
        sm[empty_idx] *= 0.0
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Update soil moisture & runoff
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.wn = sm  # daily soil moisture, mm
        self.ro = ro  # daily runoff, mm


###############################################################################
## FUNCTIONS
###############################################################################
def get_days(ts):
    """
    Name:     get_days
    Input:    datetime date (ts)
    Output:   int (delta.days)
    Features: Returns the number of days since 1 Jan 1860 for a given timestamp
    """
    # Set base time:
    base_time = datetime.date(1860, 1, 1)
    #
    # Calculate days since base time:
    delta = (ts - base_time)
    #
    # Return the number of days:
    return delta.days


def nc_history():
    """
    Name:     nc_history
    Input:    None.
    Output:   string (my_str)
    Features: Returns a string for netCDF file history field based on the
              file's creation date
    """
    my_str = "created %s" % datetime.date.today()
    return my_str


def nc_lat():
    """
    Name:     nc_lat
    Input:    None.
    Output:   numpy.ndarray, latitudes, degrees (my_lats)
    Features: Returns an array of latitudes from -90 to 90 at 0.5 deg
              resolution
    """
    my_lats = [(-90 + 0.5*0.5) + 0.5*i for i in xrange(360)]
    #
    return numpy.array(my_lats)


def nc_lon():
    """
    Name:     nc_lon
    Input:    None.
    Output:   numpy.ndarray, longitudes, degrees (my_lons)
    Features: Returns an array of longitudes from -180 to 180 at 0.5 deg
              resolution
    """
    my_lons = [(-180 + 0.5*0.5) + 0.5*i for i in xrange(720)]
    #
    return numpy.array(my_lons)


def nc_write(my_array, nc_file, nc_title, var_sname, var_lname, var_units, d):
    """
    Name:     nc_write
    Inputs:   - numpy.ndarray, 360x720 data (my_array)
              - str, netCDF output file with path (nc_file)
              - str, netCDF file title (nc_title)
              - str, variable short name (var_sname)
              - str, variable long name (var_lname)
              - str, variable units (var_units)
              - datetime.date, output date (d)
    Output:   None
    Features: Writes 360x720 data to netCDF file
    Depends:  - nc_history
              - nc_lat
              - nc_lon
              - get_days
              - kerror
    """
    # Create a new netCDF file for saving the output:
    f = netcdf.netcdf_file(nc_file, 'w')
    #
    # Add meta data variables to file:
    f.contact1 = 'tyler.davis@imperial.ac.uk'
    f.contact2 = 'c.prentice@imperial.ac.uk'
    f.history = nc_history()
    f.institution = 'Imperial College London'
    f.note1 = (
        'Monthly data is processed based on the daily iterations of the '
        'SPLASH model, (Davis et al., in prep)'
        )
    f.note2 = (
        'SPLASH model results are based on soil moisture fields initialized '
        'over the year 2001 using CRU TS 3.22 climatology'
        )
    f.note3 = (
        'Pixels without data are set equal to ' + str(kerror)
        )
    f.title = nc_title
    #
    # Create latitude dimension & variable:
    f.createDimension('lat', 360)
    lat = f.createVariable('lat', 'd', ('lat',))
    lat[:] = nc_lat()
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lat.axis = 'Y'
    #
    # Create longitude dimension & variable:
    f.createDimension('lon', 720)
    lon = f.createVariable('lon', 'd', ('lon',))
    lon[:] = nc_lon()
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'
    lon.axis = 'X'
    #
    # Create time dimension & variable
    f.createDimension('time', 1)
    t = f.createVariable('time', 'd', ('time',))
    t.standard_name = 'time'
    t.units = 'days since 1860-01-01 00:00:00'
    t[0] = get_days(d)
    #
    # Create data variable:
    data = f.createVariable(var_sname, 'd', ('time', 'lat', 'lon'))
    data._FillValue = kerror
    data.missing_value = kerror
    data.long_name = var_lname
    data.units = var_units
    data[0] = my_array
    #
    # Close and save netCDF file:
    f.close()


def writeout(f, d):
    """
    Name:     writeout
    Input:    - string, file name with path (t)
              - string, data to be written to file (d)
    Output:   None
    Features: Writes new/overwrites existing file with data string
    """
    try:
        OUT = open(f, 'w')
        OUT.write(d)
    except IOError:
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

###############################################################################
## DEFINITIONS
###############################################################################
mac = False
if mac:
    # Mac:
    cld_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_22/"
    pre_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_22/"
    tmp_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_22/"
    elv_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_00/"
    output_dir = "/Users/twdavis/Projects/data/alpha/analysis/"
else:
    # Linux:
    cld_dir = "/usr/local/share/database/cru/"
    elv_dir = "/usr/local/share/database/cru/"
    pre_dir = "/usr/local/share/database/cru/"
    swd_dir = "/usr/local/share/database/watch/netcdf/"
    tmp_dir = "/usr/local/share/database/cru/"
    output_dir = "/home/user/Projects/gepisat/data/splash/out/"

###############################################################################
## INITIALIZATIONS
###############################################################################
# Initialize data class:
my_data = DATA_G(cld_dir, elv_dir, pre_dir, tmp_dir)

# Initialize SPLASH class:
my_class = SPLASH_G(my_data.elv)
my_class.spin_up(my_data, 2001)  # <10 mins to spin up

# Define start and ending dates:
start_date = datetime.date(2002, 1, 1)
end_date = datetime.date(2003, 1, 1)

###############################################################################
## MAIN PROGRAM
###############################################################################
# Initialize monthly variables & daily soil moisture:
wn = numpy.copy(my_class.wn)
eet_mo = numpy.zeros(shape=(360, 720))
aet_mo = numpy.zeros(shape=(360, 720))

# Initialize the current day & input data:
cur_date = start_date
out_date = cur_date
my_data.read_monthly_clim(cur_date)

# Iterate through each month:
while cur_date < my_class.add_one_day(end_date):
    # Month check.
    if cur_date.year == my_data.year and cur_date.month == my_data.month:
        # Still in the same month! Continue to process!
        my_class.run_one_day(n=cur_date.timetuple().tm_yday,
                             y=cur_date.year,
                             wn=wn,
                             sf=my_data.sf,
                             tc=my_data.tair,
                             pn=my_data.pre)
        wn = numpy.copy(my_class.wn)
        eet_mo[my_data.good_idx] += my_class.eet[my_data.good_idx]
        aet_mo[my_data.good_idx] += my_class.aet[my_data.good_idx]
        #
    else:
        # New month! Calculate Cramer-Prentice alpha:
        missing_eet = numpy.where(eet_mo == 0)
        good_eet = numpy.where(eet_mo > 0)
        cpa_mo = numpy.copy(aet_mo)
        cpa_mo[good_eet] /= eet_mo[good_eet]
        cpa_mo[missing_eet] *= 0.0
        cpa_mo[missing_eet] += kerror
        #
        # Maintain error values:
        cpa_mo[my_data.noval_idx] *= 0.0
        cpa_mo[my_data.noval_idx] += kerror
        #
        # Save data to file:
        nc_output_file = "%sSPLASH_%d-%02d_alpha.nc" % (output_dir,
                                                        out_date.year,
                                                        out_date.month)
        nc_write(cpa_mo, nc_output_file, 'Monthly Cramer-Prentice alpha',
                 'CPA', 'Cramer-Prentice bioclimatic moisture index', 'none',
                 out_date)
        #
        # Update output date:
        out_date = cur_date
        #
        # Reset monthly totals
        eet_mo *= 0.0
        aet_mo *= 0.0
        #
        # Load new month's data:
        my_data.read_monthly_clim(cur_date)
        #
        # Calc new month's starting values:
        my_class.run_one_day(n=cur_date.timetuple().tm_yday,
                             y=cur_date.year,
                             wn=wn,
                             sf=my_data.sf,
                             tc=my_data.tair,
                             pn=my_data.pre)
        wn = numpy.copy(my_class.wn)
        eet_mo[my_data.good_idx] += my_class.eet[my_data.good_idx]
        aet_mo[my_data.good_idx] += my_class.aet[my_data.good_idx]
    #
    # Increment day:
    cur_date = my_class.add_one_day(cur_date)

# TEST EVAP_G
my_date = datetime.date(2001, 6, 1)
my_data.read_monthly_clim(my_date)
my_evap = EVAP_G(n=174,
                 elv=my_data.elv,
                 sf=my_data.sf,
                 tc=my_data.tair,
                 sw=0.5*numpy.ones(shape=(360, 720)),
                 y=2001)
