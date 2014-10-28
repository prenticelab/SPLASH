#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# stash_grid.py
# * based on cramer_prentice.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-30 -- created
# 2014-10-23 -- last updated
#
# ------------
# description:
# ------------
# This script calculates the monthly outputs at 0.5 deg resolution based on 
# the STASH 2.0 model
#
# IMPORTANT NOTE: 
#   Global variables are defined inside a function at runtime; therefore, you
#   must re-run function definitions if you change a global variable's value
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
# ANALYSIS OF METHODS
# 1. Kepler method:
#    drm='loutre'; lamm='kepler'; delm='loutre'
# 2. Berger's method
#    drm='loutre'; lamm='berger'; delm='loutre'
# 3. Woolf's method
#    drm='loutre'; lamm='woolf'; delm='loutre'
# 4. Klein's method
#    drm='klein'; lamm='kepler'; delm='loutre'
# 5. Spencer's method:
#    drm='loutre'; lamm='kepler'; delm='spencer'
# 6. Cooper's method
#    drm='loutre'; lamm='kepler'; delm='cooper'
# 7. Circle method
#    drm='loutre'; lamm='kepler'; delm='circle'
#
# ----------
# changelog:
# ----------
# 00. created [14.01.30]
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
#
# -----
# todo:
# -----
# * check why correct_kepler throws divide by zero warning
#   --> actually thrown in simplified_kepler, but caught in arctan function
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import glob
import numpy
from scipy.io import netcdf

###############################################################################
## GLOBAL CONSTANTS:
###############################################################################
kA = 107       # constant for Rnl (Monteith & Unsworth, 1990)
ka = 1.496e8   # semi-major axis, km
kalb_sw = 0.17 # shortwave albedo (Federer, 1968)
kalb_vis = 0.03 # visible light albedo (Sellers, 1985)
kb = 0.20      # constant for Rnl (Linacre, 1968)
kc = 0.25      # cloudy transmittivity (Linacre, 1968)
kCw = 1.05     # supply constant, mm/hr (Federer, 1982)
kd = 0.50      # angular coefficient of transmittivity (Linacre, 1968)
ke = 0.0167    # eccentricity for 2000 CE (Berger, 1978)
keps = 23.44   # obliquity for 2000 CE, degrees (Berger, 1978)
kerror = -9999 # error value
kfFEC = 2.04   # from flux to energy conversion, umol/J (Meek et al., 1984)
kG = 9.80665   # gravitational acceleration, m/s^2
kGM = 1.32712e11 # standard gravity of the sun, km^3/s^2
kGsc = 1360.8  # solar constant, W/m^2 (Kopp & Lean, 2011)
kL = 0.0065    # temperature lapse rate, K/m (Cavcar, 2000)
kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kMv = 0.01802  # molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
kPo = 101325   # base pressure, Pa (Cavcar, 2000)
kR = 8.314     # universal gas constant, J/mol/K
ksb = 5.670373e-8 # Stefan-Boltzman constant, W/m^2/K^4
kTo = 298.15   # base temperature, K (Prentice, unpublished)
kWm = 150      # soil moisture capacity, mm (Cramer & Prentice, 1988)
kw = 0.26      # entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
komega = 283.0 # longitude of perihelion for 2000 CE, degrees (Berger, 1978)

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
              Cooper, P. I. (1969). “The absorption of radiation in solar 
                stills,” Solar Energy 12.3, pp. 333–346 
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
              Loutre, M.F. (2003) Ice ages (Milankovitch theory), Encyclopedia 
                of Atmospheric Sciences, edited by J.R. Holton, J.A. Curry, 
                and J.A. Pyle, pp. 995–1003, Elsevier Ltd.
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
    def __init__(
        self, n, elv, sf, tc, sw, y=0,
        drm='loutre', lamm = 'kepler', delm = 'loutre'):
        """
        Name:     EVAP_G.__init__
        Input:    - int, day of the year (n)
                  - numpy nd.array, elevation, m (elv)
                  - numpy nd.array, fraction of sunshine hours (sf)
                  - numpy nd.array, mean daily air temperature, C (tc)
                  - numpy nd.array, evaporative supply rate, mm/hr (sw)
                  - int, year (y)
                  - string, distance method (drm)
                  - string, lambda method (lamm)
                  - string, delta method (delm)
        
        """
        # Assign default public variables:
        self.user_elv = elv
        self.user_year = y
        self.user_sf = sf
        self.user_tc = tc
        self.user_sw = sw
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
        # Convert lon and lat arrays to grids (degrees)
        #lon_grid = numpy.reshape(numpy.repeat(lon_array, 360), (360,720), 'F')
        lat_grid = numpy.reshape(numpy.repeat(lat_array, 720), (360,720), 'C')
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year, days
        #    kN, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.user_year == 0:
            self.kN = 365
        else:
            self.kN = (
                self.julian_day((self.user_year + 1), 1, 1) - 
                self.julian_day(self.user_year, 1, 1)
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes, degrees
        #    my_nu, SCALAR
        #    my_lambda, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if lamm == 'kepler':
            # Kepler Method:
            self.my_nu, self.my_lambda = self.map_days(self.user_day, 
                                                       self.user_year)
        elif lamm == 'woolf':
            # Woolf Method:
            # Eq. 23, STASH 2.0 Documentation
            woolf_b = (self.user_day - 1.0)*(360.0/self.kN)
            self.my_lambda = (
                279.9348 + woolf_b + 1.914827*self.dsin(woolf_b) - 
                0.079525*self.dcos(woolf_b) + 0.019938*self.dsin(2.0*woolf_b) - 
                0.00162*self.dcos(2.0*woolf_b)
            )
            if self.my_lambda < 0:
                self.my_lambda += 360.0
            elif self.my_lambda > 360:
                self.my_lambda -= 360.0
            self.my_nu = (self.my_lambda - komega)
            if self.my_nu < 0:
                self.my_nu += 360.0
        elif lamm == 'berger':
            # Berger'78 Method:
            xee = ke**2 
            xec = ke**3
            xse = numpy.sqrt(1.0 - xee)
            xlam = ((ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(komega) - 
                xee/4.0*(0.5 + xse)*self.dsin(2.0*komega) + 
                xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*komega))
            xlam = numpy.degrees(2.0*xlam)
            dlamm = xlam + (n - 80.0)*(360.0/self.kN)
            anm = dlamm - komega
            ranm = numpy.radians(anm)
            ranv = (ranm + (2.0*ke - xec/4.0)*numpy.sin(ranm) + 
                5.0/4.0*xee*numpy.sin(2.0*ranm) + 
                13.0/12.0*xec*numpy.sin(3.0*ranm))
            anv = numpy.degrees(ranv)
            self.my_lambda = anv + komega
            if self.my_lambda < 0:
                self.my_lambda += 360.0
            elif self.my_lambda > 360:
                self.my_lambda -= 360.0
            self.my_nu = (self.my_lambda - komega)
            if self.my_nu < 0:
                self.my_nu += 360.0
        else:
            print "Heliocentric longitude method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor, unitless
        #    dr, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if drm == 'loutre':
            # Eq. 7--9, STASH 2.0 Documentation
            my_rho = (1.0 - ke**2.0)/(1.0 + ke*self.dcos(self.my_nu))
            self.dr = (1.0/my_rho)**2.0
        elif drm == 'klein':
            # Eq. 11, STASH 2.0 Documentation
            self.dr = (
                1.0 + 2.0*ke*numpy.cos(2.0*numpy.pi*self.user_day/self.kN)
            )
        else:
            print "Distance factor method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle, degrees
        #    delta, SCALAR
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if delm == 'loutre':
            # Eq. 14, STASH 2.0 Documentation (Loutre, 2002)
            self.delta = numpy.arcsin(self.dsin(self.my_lambda)*self.dsin(keps))
            self.delta *= (180.0/numpy.pi)
        elif delm == 'cooper':
            # Eq. 16, STASH 2.0 Documentation (Cooper, 1969)
            self.delta = keps*numpy.sin(
                2.0*numpy.pi*(komega + self.user_day)/self.kN
            )
        elif delm == 'circle':
            # Eq. 15, STASH 2.0 Documentation
            self.delta = -1.0*keps*numpy.cos(
                2.0*numpy.pi*(self.user_day + 10.0)/self.kN
            )
        elif delm == 'spencer':
            # Eq. 17, STASH 2.0 Documentation
            spencer_b = (self.user_day - 1.0)*(2.0*numpy.pi/self.kN)
            self.delta = (0.006918 - 
                          0.399912*numpy.cos(spencer_b) + 
                          0.070257*numpy.sin(spencer_b) - 
                          0.006758*numpy.cos(2.0*spencer_b) +
                          0.000907*numpy.sin(2.0*spencer_b) - 
                          0.0022697*numpy.cos(3.0*spencer_b) + 
                          0.00148*numpy.sin(3.0*spencer_b))
            self.delta *= (180.0/numpy.pi)
        else:
            print "Declination angle method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes, unitless
        #    ru, MATRIX (360x720)
        #    rv, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(self.delta)*self.dsin(lat_grid)
        rv = self.dcos(self.delta)*self.dcos(lat_grid)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the sunset hour angle, degrees
        #    hs, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed hour angle grid 
        hs = numpy.zeros(shape=(360,720))
        #
        # Indexes of pixels under polar night conditions:
        # Note: hs == 0 degrees (no further comp's)
        # hs_neg = numpy.where(ru/rv <= -1.0)
        #
        # Indexes of pixels under polar day conditions:
        # Note: hs == 180 degrees
        hs_pos = numpy.where(ru/rv >= 1.0)
        for j in xrange(len(hs_pos[0])):
            (a,b) = (hs_pos[0][j], hs_pos[1][j])
            hs[a,b] = 180.0
        #
        # Indexes of pixels for regular conditions:
        # Note: hs = acos(-u/v)
        hs_reg = numpy.where((ru/rv < 1.0) & (ru/rv > -1.0))
        for j in xrange(len(hs_reg[0])):
            (a,b) = (hs_reg[0][j], hs_reg[1][j])
            # Eq. 37, STASH 2.0 Documentation
            hs[a,b] = (180.0/numpy.pi)*numpy.arccos(-1.0*ru[a,b]/rv[a,b])
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate daily extraterrestrial solar radiation, J/m^2
        #    ra_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 40, STASH 2.0 Documentation
        # Note: ru = sin(delta)*sin(phi); rv = cos(delta)*cos(phi)
        self.ra_d = (86400.0/numpy.pi)*kGsc*self.dr*(
            (ru*hs)*(numpy.pi/180.0) + (rv*self.dsin(hs))
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate transmittivity, unitless
        #    tau, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 48, STASH 2.0 Documentation
        tau_o = kc + (kd*self.user_sf)
        #
        # Eq. 49, STASH 2.0 Documentation
        # Note: elv missing == -999
        tau = tau_o*(1.0 + (2.67e-5)*self.user_elv)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        #    ppfd_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 57, STASH 2.0 Documentation
        self.ppfd_d = (kfFEC*1.0e-6)*(1.0 - kalb_vis)*(tau*self.ra_d)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation, W/m^2
        #     rnl, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 46, STASH 2.0 Documentation
        # Note: missing user_tc == -9999
        rnl = (
            (kb + (1.0 - kb)*self.user_sf)*(kA - self.user_tc)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate variable substitute, W/m^2
        #     rw, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = (1.0 - kalb_sw)*tau*kGsc*self.dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate net radiation cross-over hour angle, degrees
        #     hn, MATRIX (360,720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed cross-over angle: 
        hn = numpy.zeros(shape=(360,720))
        #
        # Indexes of pixels where Rnl is all-day negative:
        # hn == 0 (no further processing)
        #hn_neg = numpy.where((rnl - rw*ru)/(rw*rv) >= 1.0)
        #
        # Indexes of pixels where Rnl is all-day positive:
        hn_pos = numpy.where((rnl - rw*ru)/(rw*rv) <= -1.0)
        for j in xrange(len(hn_pos[0])):
            (a,b) = (hn_pos[0][j], hn_pos[1][j])
            hn[a,b] = 180.0
        #
        # Indexes of pixels for regular Rnl:
        hn_reg = numpy.where(
            ((rnl - rw*ru)/(rw*rv) < 1.0) & 
            ((rnl - rw*ru)/(rw*rv) > -1.0)
        )
        for j in xrange(len(hn_reg[0])):
            (a,b) = (hn_reg[0][j], hn_reg[1][j])
            # Eq. 51, STASH 2.0 Documentation
            hn[a,b] = (180.0/numpy.pi)*numpy.arccos(
                (rnl[a,b] - rw[a,b]*ru[a,b])/(rw[a,b]*rv[a,b])
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate daytime net radiation, J/m^2
        #     rn_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 53, STASH 2.0 Documentation
        rn_d = (86400.0/numpy.pi)*(
            (rw*ru - rnl)*hn*(numpy.pi/180.0) + rw*rv*self.dsin(hn)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate nighttime net radiation, J/m^2
        #     rnn_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 56, STASH 2.0 Documentation
        rnn_d = (86400.0/numpy.pi)*(
            rw*ru*(hs - hn)*(numpy.pi/180.0) + 
            rw*rv*(self.dsin(hs) - self.dsin(hn)) + rnl*(
                numpy.pi - 2.0*hs*(numpy.pi/180.0) + hn*(numpy.pi/180.0)
            )
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 15. Calculate water-to-energy conversion, m^3/J
        #     econ, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Slope of saturation vap press temp curve, Pa/K
        s = self.sat_slope(self.user_tc)
        # Enthalpy of vaporization, J/kg
        lv = self.enthalpy_vap(self.user_tc)
        # Density of water, kg/m^3
        pw = self.density_h2o(self.user_tc, self.elv2pres(self.user_elv))
        # Psychrometric constant, Pa/K
        g = self.psychro(self.user_tc, self.elv2pres(self.user_elv))
        #
        # Eq. 58, STASH 2.0 Documentation
        econ = s/(lv*pw*(s + g))
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 16. Calculate daily condensation, mm
        #     wc, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 69, STASH 2.0 Documentation
        self.wc = 1000.0*econ*numpy.abs(rnn_d)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 17. Estimate daily EET, mm
        #     eet_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 71, STASH 2.0 Documentation
        self.eet_d = 1000.0*(econ*rn_d)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 18. Estimate daily PET, mm
        #     pet_d, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 73, STASH 2.0 Documentation
        self.pet_d = (1.0 + kw)*self.eet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 19. Calculate variable substitute, (mm/hr)/(W/m^2)
        #     rx, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = 1000.0*3600.0*(1.0 + kw)*econ
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 20. Calculate the intersection hour angle, degrees
        #     hi, MATRIX (360x720)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed intersection hour angle grid
        hi = numpy.zeros(shape=(360,720))
        #
        # Compute cos(hi)
        # Eq. 79, STASH 2.0 Documentation
        cos_hi = (
            self.user_sw/(rw*(rv*rx)) + rnl/(rw*rv) - ru/rv
        )
        #
        # Indexes where supply exceeds demand:
        # Note: hi == 0 (no further comp's)
        #hi_neg = numpy.where(cos_hi >= 1.0)
        #
        # Indexes where supply limits demand everywhere
        # Note: hi == 180
        hi_pos = numpy.where(cos_hi <= -1.0)
        for j in xrange(len(hi_pos[0])):
            (a,b) = (hi_pos[0][j], hi_pos[1][j])
            hi[a,b] = 180.0
        #
        # Indexes for regular supply
        hi_reg = numpy.where((cos_hi < 1.0) & (cos_hi > -1.0))
        for j in xrange(len(hi_reg[0])):
            (a,b) = (hi_reg[0][j], hi_reg[1][j])
            hi[a,b] = numpy.arccos(cos_hi[a,b])*(180.0/numpy.pi)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 21. Estimate daily AET (aet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 82, STASH 2.0 Documentation
        self.aet_d = (24.0/numpy.pi)*(
            self.user_sw*hi*(numpy.pi/180.0) + 
            rx*(rw*rv)*(self.dsin(hn) - self.dsin(hi)) + 
            (rx*(rw*ru) - rx*rnl)*(hn - hi)*(numpy.pi/180.0)
        )
    #
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
    #
    def dsin(self, x):
        """
        Name:     EVAP_G.dsin
        Input:    float/nd.array, angle, degrees (x)
        Output:   float/nd.array, sin(x*pi/180)
        Features: Calculates the sine of angle(s) given in degrees
        """
        return numpy.sin(x*numpy.pi/180.0)
    #
    def dtan(self, x):
        """
        Name:     EVAP_G.dtan
        Input:    float/nd.array, angle, degrees (x)
        Output:   float/nd.array, tan(x*pi/180)
        Features: Calculates the tangent of angle(s) given in degrees
        """
        return numpy.tan(x*numpy.pi/180.0)
    #
    def get_lon_lat(self,x,y,r):
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
    #
    def earth_velocity(self, lon):
        """
        Name:     EVAP_G.earth_velocity
        Input:    longitude(s) w.r.t. the perihelion (lon)
        Output:   angular velocity(ies), rad/day (w)
        Features: Calculates earth's angular velocity at a given longitude
                  relative to the perihelion
        Depends:  Global constants:
                  - ke
        Ref:      Kepler's Second Law of Planetary Motion
        """
        w = (
            2.0*numpy.pi/self.kN*(
                1.0 + ke*self.dcos(lon)
                )**2.0*(1.0 - ke**2)**(-1.5)
        )
        return w
    #
    def julian_day(self,y,m,i):
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
    #
    def equinox(self, year, opt=0):
        """
        Name:     EVAP_G.equinox
        Input:    - int, year (year)
                  - int, option (opt)
                    0: vernal equinox     1: summer solstice
                    2: autumnal equinox   3: winter solstice
        Output:   float, day of the year
        Features: Calculates the day of the year on which seasonal dates fall
        Depends:  julian_day
        Ref:      J. Meeus (1991), Ch.26 "Equinoxes and solstices," 
                  Astronomical Algorithms
        """
        # Table 26.C (Meeus, 1991)
        periodic_terms = numpy.array([
            # A    B          C
            485, 324.96,   1934.136,
            203, 337.23,  32964.467,
            199, 342.08,     20.186,
            182,  27.85, 445267.112,
            156,  73.14,  45036.886,
            136, 171.52,  22518.443,
            77, 222.54,  65928.934,
            74, 296.72,   3034.906,
            70, 243.58,   9037.513,
            58, 119.81,  33718.147,
            52, 297.17,    150.678,
            50,  21.02,   2281.226,
            45, 247.54,  29929.562,
            44, 325.15,  31555.956,
            29,  60.93,   4443.417,
            18, 155.12,  67555.328,
            17, 288.79,   4562.452,
            16, 198.04,  62894.029,
            14, 199.76,  31436.921,
            12,  95.39,  14577.848,
            12, 287.11,  31931.756,
            12, 320.81,  34777.259,
            9, 227.73,   1222.114,
            8,  15.45,  16859.074
            ])
        #
        # Table 26.A (Meeus, 1991)
        jde_table_a = numpy.array([
            numpy.array([1721139.29189, 365242.13740,  0.06134,  
                         0.00111, -0.00071]), # March equinox
            numpy.array([1721233.25401, 365241.72562, -0.05323,  
                         0.00907,  0.00025]), # June solstice
            numpy.array([1721325.70455, 365242.49558, -0.11677, 
                         -0.00297,  0.00074]), # September equinox
            numpy.array([1721414.39987, 365242.88257, -0.00769, 
                         -0.00933, -0.00006])  # December solstice
            ])
        #
        # Table 26.B (Meeus, 1991)
        jde_table_b = numpy.array([
            numpy.array([2451623.80984, 365242.37404,  0.05169, 
                         -0.00411, -0.00057]), # March equinox
            numpy.array([2451716.56767, 365241.62603,  0.00325,  
                         0.00888, -0.00030]), # June solstice
            numpy.array([2451810.21715, 365242.01767, -0.11575,  
                         0.00337,  0.00078]), # September equinox
            numpy.array([2451900.05952, 365242.74049, -0.06223, 
                         -0.00823,  0.00032])  # December solstice
            ])
        #
        if year < 1000:
            # Use Table 26.A for years -1000 to 1000
            jde_table = jde_table_a
            y = year/1000.0
        else:
            # Use Table 26.B for years 1000 to 3000
            jde_table = jde_table_b
            y = (year - 2000.0)/1000.0
        #
        # Calculate the mean equinox based on Table 26.A or 26.B:
        jde0 = (
       	jde_table[opt,0] +
       	jde_table[opt,1]*y +
       	jde_table[opt,2]*y*y +
       	jde_table[opt,3]*y*y*y +
       	jde_table[opt,4]*y*y*y*y
        )
        #
        # Calculate the other three terms:
        t = (jde0 - 2451545.0)/36525.0
        w = (35999.373*t) - 2.47
        dl = (
            1.0 + 
            (0.0334*self.dcos(w)) + 
            (0.0007*self.dcos(2*w))
        )
        #
        # Compute the sum of the 24 periodic terms:
        s = 0
        j = 0
        for i in xrange(24):
            s += periodic_terms[j]*self.dcos(
                (periodic_terms[j+1] + (periodic_terms[j+2]*t))
            )
            j += 3
        #
        # Calculate the JDE (Julian Ephemeris Day) of the equinox/solstice:
        jde = jde0 + ((s*0.00001)/dl)
        #
        # Calcuate the JDE for the first day of this year:
        jde_year = self.julian_day(year, 1, 1)
        #
        # Return the day of year:
        return (jde - jde_year + 1)
    #
    def correct_kepler(self, lon, t):
        """
        Name:     EVAP_G.correct_kepler
        Input:    - longitudes w.r.t. the perihelion; at least two points 
                    before and one point after the 180 degree cross-over (lon)
                  - days associates with the longitudes (t)
        Output:   days
        Features: Corrects the array of days output from simplified_kepler due 
                  to the asymptote in the arctan at 180 degrees
        """
        # Find indexes of positive and negative values:
        neg_points = numpy.where(t<0)[0]
        pos_points = numpy.where(t>=0)[0]
        #
        # Linearly extrapolate the positive position of the first negative 
        # value for:
        #   xa ... ya  <-- known value pair a
        #   xb ... yb  <-- known value pair b
        #   xi ... yi  <-- extrapolated value pair i
        # where all variables are known except for yi, such that:
        #   (yi-ya)/(yb-ya) = (xi-xa)/(xb-xa)
        # and solving for yi:
        #   yi = ya + (yb-ya)*(xi-xa)/(xb-xa)
        xa = pos_points[-2]
        xb = pos_points[-1]
        xi = neg_points[0]
        #
        ya = t[xa]
        yb = t[xb]
        yi = ya + (yb - ya)*(xi - xa)/(xb - xa)
        #
        # Calculate the difference between the pos. and neg. values at y:
        diff_y = yi - t[xi]
        #
        # Add the difference to the negative values to push them to the proper
        # positive position, i.e.:
        #
        #      +y                     +y
        #       |     o               |           o
        #       |   o                 |         o  
        #       | o                   |       o    
        #       |------------- +x  => |     o      
        #       |           o         |   o        
        #       |         o           | o          
        #       |       o             |------------ +x
        #
        #       (a) Original Data     (b) Corrected
        #
        y = numpy.copy(t)
        y[neg_points] = y[neg_points] + diff_y
        #
        return (y)
    #
    def simplified_kepler(self,lon):
        """
        Name:     EVAP_G.simplified_kepler
        Input:    float, longitude w.r.t. the perihelion (lon)
        Output:   float, days traveled from the perihelion (t)
        Features: Calculates the time (days) of earth's orbit around the sun
                  for a given longitude using a simplified Kepler approach
        Depends:  correct_kepler
                  Global constants:
                  - ke
        Ref:      Kepler's Second Law of Planetary Motion
        """
        # Make lon into numpy array, if it is not already:
        if isinstance(lon, numpy.float) or isinstance(lon, numpy.int):
            lon = numpy.array([lon,])
        elif not isinstance(lon, numpy.ndarray):
            lon = numpy.array(lon)
        #
        # Calculate the days
        # NOTE: arctan(inf) = pi/2
        t = (
            self.kN/numpy.pi*(
                numpy.arctan(
                    self.dsin(lon)/(self.dcos(lon) + 1)
                ) - ke*self.dsin(lon)
            )
        )
        #
        # Correct the days, if necessary:
        if len(lon) > 1:
            t = self.correct_kepler(lon, t)
        #
        return t
    #
    def map_days(self, n, y):
        """
        Name:     EVAP_G.map_days
        Input:    - int, day of the year (n)
                  - int, year (y)
        Output:   float, longitude relative to the vernal equinox (lamda_doy)
        Features: Computes earth's longitude relative to the vernal equinox 
                  for a given day of the year
        Depends:  Functions:
                  - earth_period
                  - earth_velocity
                  - equinox
                  - simplified_kepler
                  Global constants:
                  - komega
        """
        # Create longitude field and compute the days for orbit:
        lon_nu = numpy.array([1.0*i for i in xrange(361)])
        day_nu = self.simplified_kepler(lon_nu)
        #
        # Compute the angle of the vernal equinox w.r.t. the perihelion
        # i.e., the explementary angle of komega:
        wp = (360.0 - komega)
        #
        # Calculate the length of time it takes earth to travel from the
        # perihelion (t=0) to the vernal equinox:
        if wp < 180.0:
            days_to_ve = self.simplified_kepler(wp)[0]
        else:
            lon_temp = numpy.array([179.0, 180.0, 181.0, wp])
            days_to_ve = self.simplified_kepler(lon_temp)[3]
        #
        # Get day of year of vernal equinox
        if y == 0:
            day_ve = 80.0
        else:
            day_ve = self.equinox(y)
        #
        # Calculate the offset between days to and day of vernal equinox:
        offset = (day_ve - days_to_ve)
        #
        # Calculate the calendar days and set between 0 and kN:
        calendar_days = (day_nu + offset)
        calendar_days[numpy.where(calendar_days >= self.kN)] -= self.kN
        #
        # Check to see if n is listed in calendar:
        if n in calendar_days:
            icalendar = numpy.where(calendar_days==n)[0][0]
            nu_doy = lon_nu[icalendar]
        else:
            # Find the calendar day the precedes doy:
            calendar_temp = numpy.sort(numpy.append(calendar_days, [n,]))
            dbefore = calendar_temp[numpy.where(calendar_temp==n)[0][0]-1]
            #
            # Get the index of the preceding day:
            ibefore = numpy.where(calendar_days == dbefore)[0][0]
            #
            # Get the angular velocity for the longitude of the preceding day:
            vbefore = self.earth_velocity(lon_nu[ibefore])
            #
            # Calculate the delta time
            dt = (n - dbefore)
            #
            # Calculate the new longitude, degrees:
            nu_doy = lon_nu[ibefore] + (vbefore*dt)*180.0/numpy.pi
        #
        # Convert nu to lambda:
        lamda_doy = nu_doy + komega
        if lamda_doy >= 360:
            lamda_doy -= 360
        #
        return (nu_doy, lamda_doy)
    #
    def sat_slope(self, tc):
        """
        Name:     EVAP_G.sat_slope
        Input:    numpy nd.array, air temperatures (tc), degrees C
        Output:   numpy nd.array, slopes of sat vap press temp curve (s)
        Features: Calculates the slopes of the sat pressure temp curve, Pa/K
        Ref:      Eq. 13, Allen et al. (1998)
        """
        s = (17.269)*(237.3)*(610.78)*(
            numpy.exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2.0)
        )
        return s
    #
    def enthalpy_vap(self, tc):
        """
        Name:     EVAP_G.enthalpy_vap
        Input:    numpy nd.array, air temperatures (tc), degrees C
        Output:   numpy nd.array, latent heats of vaporization
        Features: Calculates the enthalpy of vaporization, J/kg
        Ref:      Eq. 8, Henderson-Sellers (1984)
        """
        return (1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2.0)
    #
    def elv2pres(self, z):
        """
        Name:     EVAP_G.elv2pres
        Input:    numpy nd.array, elevations above sea level (z), m
        Output:   numpy nd.array, atmospheric pressures, Pa
        Features: Calculates atm. pressure for a given elevation
        Depends:  Global constants
                  - kPo
                  - kTo
                  - kL
                  - kMa
                  - kG
                  - kR
        Ref:      Cavcar (2000)
        """
        p = kPo*(1.0 - (kL*z)/kTo)**(kG*kMa/(kR*kL))
        return p
    #
    def density_h2o(self, tc, p):
        """
        Name:     EVAP_G.density_h2o
        Input:    - numpy nd.array, air temperatures (tc), degrees C
                  - numpy nd.array, atmospheric pressures (p), Pa
        Output:   numpy nd.array, densities of water, kg/m^3
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
            (1.0e3)*po*(ko + ca*pbar + cb*(pbar**2.0))/(
                ko + ca*pbar + cb*(pbar**2.0) - pbar
            )
        )
        return pw
    #
    def psychro(self, tc, p):
        """
        Name:     EVAP_G.psychro
        Input:    - numpy nd.array, air temperatures (tc), degrees C
                  - numpy nd.array, atm. pressures (p), Pa
        Output:   numpy nd.array, psychrometric constant, Pa/K
        Features: Calculates the psychrometric constant for a given temperature
                  and pressure
        Depends:  Global constants:
                  - kMa
                  - kMv
        Refs:     Allen et al. (1998)
                  Tsilingiris (2008) 
        """
        # Calculate the specific heat capacity of water, J/kg/K
        # Eq. 47, Tsilingiris (2008)
        cp = (1.0e3)*(
            1.0045714270 +
            (2.050632750e-3)*tc -
            (1.631537093e-4)*tc*tc +
            (6.212300300e-6)*tc*tc*tc -
            (8.830478888e-8)*tc*tc*tc*tc +
            (5.071307038e-10)*tc*tc*tc*tc*tc
        )
        #
        # Calculate latent heat of vaporization, J/kg
        lv = self.enthalpy_vap(tc)
        #
        # Calculate psychrometric constant, Pa/K
        # Eq. 8, Allen et al. (1998)
        return (cp*kMa*p/(kMv*lv))
    #

###############################################################################
## FUNCTIONS 
###############################################################################
def add_one_day(dt0):
    """
    Name:     add_one_day
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt1)
    Features: Adds one day to datetime
    """
    dt1 = dt0 + datetime.timedelta(days=1) 
    return dt1
    
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

def add_one_year(dt0):
    """
    Name:     add_one_year
    Input:    datetime date
    Output:   datetime date
    Features: Adds one year to the datetime preserving calendar date, if it 
              exists, otherwise uses the following day (i.e., February 29 
              becomes March 1)
    Ref:      G. Rees (2013) Stack Overflow
              http://stackoverflow.com/questions/15741618/add-one-year-in-
              current-date-python
    """
    try:
        return dt0.replace(year = dt0.year + 1)
    except ValueError:
        return dt0 + (
            datetime.date(dt0.year + 1, 1, 1) - datetime.date(dt0.year, 1 ,1)
        )

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

def get_year_days(ts):
    """
    Name:     get_year_days
    Input:    datetime date
    Output:   int
    Features: Returns the total number of days in the year
    Depends:  add_one_year
    """
    ts1 = datetime.date(ts.year, 1 , 1)
    ts2 = add_one_year(ts1)
    return (ts2 - ts1).days

def get_elevation(d):
    """
    Name:     get_elevation
    Input:    char, directory to CRU netcdf file (d)
    Output:   numpy nd.array
    Features: Reads CRU TS 3.00 0.5 deg. elevation data from netcdf file
    """
    # Read directory for elevation file:
    my_file = glob.glob(d + "*dat")[0]
    #
    # Open file and read data:
    # NOTE: array is shape: 360 x 720
    #      'lat' goes from -89.75 -- 89.75 (south to north)
    #      'lon' goes from -179.75 -- 179.75 (east to west)
    #      missing data = -999.0
    f = numpy.loadtxt(my_file)
    return f

def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime date, base timestamp
              - datetime date, current timestamp
              - numpy nd.array, days since base timestamp
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
        #       where v is 'tmp', 'pre', 'cld'
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

def save_to_file(d, f):
    """
    Name:     save_to_file
    Input:    - numpy nd.array (d)
              - string, file name with path (f)
    Output:   None
    Features: Writes data to file in ASCII raster format (1000 x value) with 
              missing values set to kerror
    Depends:  - writeout
              - kerror
    """
    # Define header line for ASCII raster:
    header = (
        "NCOLS 720\n"
        "NROWS 360\n"
        "XLLCORNER -180.0\n"
        "YLLCORNER -90.0\n"
        "CELLSIZE 0.5\n"
        "NODATA_VALUE %s\n"
        ) % (kerror)
    #
    # Save header line to file:
    writeout(f, header)
    #
    # Iterate through data:
    for i in xrange(360):
        # Reverse column ordering:
        y = 359 - i
        #
        # Reset outline for next row:
        outline = ""
        #
        for x in xrange(720):
            val = d[y,x]
            #
            if val != kerror:
                val = (1e3)*val
                #
            # Add value to outline:
            outline = "%s%d " % (outline, int(val))
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(f, 'a')
        OUT.write(outline)
        OUT.close()

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

def spin_up(y, w_init, thresh, elv, lc, ef, tmp_d, cld_d, pre_d):
    """
    Name:     spin_up
    Input:    - int, year (y)
              - float, soil moisture initial value (w_init)
              - float, error value (thresh)
              - numpy nd.array, elevation array (elv)
              - numpy nd.array, land clip (lc)
              - numpy nd.array, error field (ef)
              - string, CRU TS tmp netcdf directory (tmp_d)
              - string, CRU TS cld netcdf directory (cld_d)
              - string, CRU TS pre netcdf directory (pre_d)
    Output:   numpy nd.array
    Features: Returns soil moisture array following a spin-up on CRU TS 
              meteorological data of a given year
    """
    # Initialize daily soil moisture for one year:
    w = w_init*(numpy.ones(shape=(366,360,720)))
    #
    # Spin up iteration counter and max iterations:
    spin_count = 0
    spin_max = 10
    #
    # Define start and end dates based on given year:
    start_date = datetime.date(y, 1, 1)
    end_date = datetime.date(y+1, 1, 1)
    #
    # Initialize error field mean:
    w_err = (1.0 + thresh)
    #
    # Run yearly spin if errors are still large and max iterations not reached:
    while w_err > thresh and spin_count < spin_max:
        # Increment spin counter:
        spin_count += 1
        #
        # Initialize iteration:
        cur_date = start_date
        ny = get_year_days(start_date)
        while cur_date < end_date:
            nm = get_month_days(cur_date)
            #
            tmp = get_monthly_cru(tmp_d, cur_date, 'tmp')
            pre = get_monthly_cru(pre_d, cur_date, 'pre')
            cld = get_monthly_cru(cld_d, cur_date, 'cld')
            #
            tair = (tmp*lc) + ef
            ppt = (pre*lc)
            sf = (1.0 - cld/100.0)*lc
            #
            cur_day = 0
            while cur_day < nm:
                n = (cur_date.timetuple().tm_yday + cur_day)
                idx = int(n-1)-1
                if idx < 0:
                    idx = int(ny-1)
                sw = (kCw/kWm)*w[idx,:, :]
                my_evap = EVAP_G(n, elv, sf, tair, sw, y=cur_date.year)
                ro = w[idx,:, :] + ppt/(1.0*nm) + my_evap.wc - my_evap.aet_d
                #
                ro_full = numpy.where(ro >= kWm)
                ro_empty = numpy.where(ro <= 0)
                ro_reg = numpy.where((ro < kWm) & (ro > 0))
                #
                for j in xrange(len(ro_full[0])):
                    (a,b) = (ro_full[0][j], ro_full[1][j])
                    w[(n-1),a,b] = kWm
                for j in xrange(len(ro_empty[0])):
                    (a,b) = (ro_empty[0][j], ro_empty[1][j])
                    w[(n-1),a,b] = 0.0
                for j in xrange(len(ro_reg[0])):
                    (a,b) = (ro_reg[0][j], ro_reg[1][j])
                    w[(n-1),a,b] = ro[a,b]
                #
                cur_day += 1
            #
            cur_date = add_one_month(cur_date)
        #
        # Calc mean error between 31 December and 1 January
        w_err = numpy.abs((w[-1,:,:] - w[0,:,:]).mean())
    #
    #print 'spun', spin_count, 'years'
    return w

def mean_monthly_w(w, m, y):
    """
    Name:     mean_monthly_w
    Inputs:   - numpy.ndarray, daily soil moisture: 366x360x720 (w)
              - int, month, 1..12 (m)
              - int, year (y)
    Outputs:  numpy.ndarray
    Features: Returns mean monthly soil moisture for a given month and year
    Depends:  get_month_days
    """
    # Initialize return value:
    mean_w = numpy.zeros(shape=(360,720))
    #
    # Create list of days of the year to process:
    try:
        my_date = datetime.date(y, m, 1)
    except ValueError:
        print "Month must be 1 to 12"
    else:
        nm = get_month_days(my_date)
        n = [my_date.timetuple().tm_yday + i for i in xrange(nm)]
        #
        # Begin summation over days:
        for i in n:
            mean_w += w[(i-1), :, :]
        #
        # Calculate mean:
        mean_w /= float(nm)
        #
        return mean_w
    #

###############################################################################
## DEFINITIONS 
###############################################################################
mac = 0   # 1: TRUE, 0: FALSE
if mac:
    # Mac:
    cld_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_21/"
    pre_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_21/"
    tmp_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_21/"
    elv_dir = "/Users/twdavis/Projects/data/cru/cru_ts_3_00/"
    output_dir = "/Users/twdavis/Projects/data/alpha/analysis/"
else:
    # Linux:
    cld_dir = "/usr/local/share/database/cru/"
    elv_dir = "/usr/local/share/database/cru/elv/"
    pre_dir = "/usr/local/share/database/cru/"
    swd_dir = "/usr/local/share/database/watch/netcdf/"
    tmp_dir = "/usr/local/share/database/cru/"
    output_dir = "/home/user/Projects/gepisat/data/stash/out/"

an_dic = {
    "AN1" : ('loutre', 'kepler', 'loutre'),  # Simplified Kepler
    "AN2" : ('loutre', 'berger', 'loutre'),  # Berger's Method
    "AN3" : ('loutre', 'woolf', 'loutre'),   # Woolf'd Method
    "AN4" : ('klein', 'kepler', 'loutre'),   # Klein's Method
    "AN5" : ('loutre', 'kepler', 'spencer'), # Spencer's Method
    "AN6" : ('loutre', 'kepler', 'cooper'),  # Cooper's Method
    "AN7" : ('loutre', 'kepler', 'circle')   # Circle Method
}

###############################################################################
## INITIALIZATIONS
###############################################################################
# Initialize elevation 360x720 data (meters):
# Note: missing data = -999
elv = get_elevation(elv_dir)

# Get an array to clip land values (all others set to zero):
land_clip = (elv - elv.min()).clip(max=1)

# Get error values (-9999 where clip is zero)
error_field = kerror*(1 - land_clip)

# Set 'ocean' elevations equal to zero:
elv *= land_clip

# Initialize monthly key outputs:
ppfd_mo = numpy.zeros(shape=(360,720))  # PPFD, mol/m^2
aet_mo = numpy.zeros(shape=(360,720))   # AET, mm
eet_mo = numpy.zeros(shape=(360,720))   # EET, mm
pet_mo = numpy.zeros(shape=(360,720))   # PET, mm

# Annual total top-of-atmosphere PPFD, Qo_ann
qo_ann = numpy.zeros(shape=(360,720))

# Initialize soil moisture:
# * w.shape: days, lat, lon
w = spin_up(2000, kWm, 0.05, elv, land_clip, 
            error_field, tmp_dir, cld_dir, pre_dir)
#w = kWm*(numpy.ones(shape=(366,360,720)))

# Save initialization (optional for mult. processing)
w_bak = numpy.copy(w)

# Starting/ending dates:
start_date = datetime.date(2001, 1, 1)
end_date = datetime.date(2002,1,1)

# Output date:
print_date = datetime.date(2000,1,1)

###############################################################################
## MAIN PROGRAM 
###############################################################################
# Initialize current monthly timestamp
cur_date = start_date
w = numpy.copy(w_bak)
out_no = "AN1"

# ~~~~~~~~~~~~~~~~~~~
# MONTHLY timestamps:
# ~~~~~~~~~~~~~~~~~~~
while cur_date < end_date:
    # Get number of days in the current month and year:
    nm = get_month_days(cur_date)
    ny = get_year_days(cur_date)
    #
    # Open and read tmp, pre, & cld data for this month:
    # NOTE: all the gridded data (CRU TS) have the same
    #       shape (360x720) and are indexed identically, 
    #       row-major from (-179.75, -89.57) to (179.75, 89.75)
    tmp = get_monthly_cru(tmp_dir, cur_date, 'tmp')
    pre = get_monthly_cru(pre_dir, cur_date, 'pre')
    cld = get_monthly_cru(cld_dir, cur_date, 'cld')
    #
    # Convert cloudiness to fractional sunshine, sf
    sf = (1.0 - cld/100.0)*land_clip
    #
    # Reduce error field on the temperature grid (i.e., 9.9e36 -> -9999)
    tair = (tmp*land_clip) + error_field
    #
    # Clip precipitation field:
    ppt = (pre*land_clip)
    #
    # Reset monthly fields for next monthly total:
    ppfd_mo *= 0.0
    aet_mo *= 0.0
    eet_mo *= 0.0
    pet_mo *= 0.0
    #
    # ~~~~~~~~~~~~~~~~
    # DAILY timesteps:
    # ~~~~~~~~~~~~~~~~
    cur_day = 0
    while cur_day < nm:
        # Get the Julian day of the year:
        n = (cur_date.timetuple().tm_yday + cur_day)
        #
        # Get the index for yesterday:
        idx = int(n-1)-1
        if idx < 0:
            idx = int(ny-1)
        #
        # Calculate the (360x720) evaportive supply rate (mm/hr)
        # Eq. 4, STASH 2.0 Documentation
        sw = (kCw/kWm)*w[idx,:, :]
        #
        # Calculate the (360x720) daily evaporation and radiation values
        my_evap = EVAP_G(n, elv, sf, tair, sw, 
                         y=0,          # cur_date.year
                         drm=an_dic[out_no][0],    # loutre, klein
                         lamm=an_dic[out_no][1],   # kepler, woolf, berger
                         delm=an_dic[out_no][2])   # loutre, circle, cooper, sp.
        #
        # Update daily soil moisture:
        # Eq. 83, STASH 2.0 Documentation
        # NOTE: assume precip is monthly total (i.e., mm/mo);
        #       divide by the number of days to get daily precip
        ro = w[idx,:, :] + ppt/(1.0*nm) + my_evap.wc - my_evap.aet_d
        #
        # Indexes where ro exceeds bucket capacity (full)
        ro_full = numpy.where(ro >= kWm)
        for j in xrange(len(ro_full[0])):
            (a,b) = (ro_full[0][j], ro_full[1][j])
            w[(n-1),a,b] = kWm
        #
        # Indexes where ro exceeds bucket capacity (empty)
        ro_empty = numpy.where(ro <= 0)
        for j in xrange(len(ro_empty[0])):
            (a,b) = (ro_empty[0][j], ro_empty[1][j])
            w[(n-1),a,b] = 0.0
        #
        # Indexes for regular ro
        ro_reg = numpy.where((ro < kWm) & (ro > 0))
        for j in xrange(len(ro_reg[0])):
            (a,b) = (ro_reg[0][j], ro_reg[1][j])
            w[(n-1),a,b] = ro[a,b]
        #
        # Increment daily aet & eet to monthly values:
        ppfd_mo += my_evap.ppfd_d
        aet_mo += my_evap.aet_d
        eet_mo += my_evap.eet_d
        pet_mo += my_evap.pet_d
        #
        # Output daily solar radiation (top of atmos. PPFD)
        if 0:
            qo_ann += (1e-6)*(kfFEC*my_evap.ra_d) # mol/m^2
        #
        # Update the day of the month:
        cur_day += 1
        #
    # Calculate Cramer-Prentice alpha:
    cpa_mo = numpy.zeros(shape=(360,720))
    #
    # Find where denominator is zero (can't divide):
    cp_iz = numpy.where(eet_mo == 0)
    for j in xrange(len(cp_iz[0])):
        (a,b) = (cp_iz[0][j], cp_iz[1][j])
        cpa_mo[a,b] = kerror
    #
    # Find where denominator is not zero (can divide)
    cp_nz = numpy.where(eet_mo != 0)
    for j in xrange(len(cp_nz[0])):
        (a,b) = (cp_nz[0][j], cp_nz[1][j])
        cpa_mo[a,b] = (aet_mo[a,b]/eet_mo[a,b])
    cpa_mo = (cpa_mo*land_clip) + error_field
    #
    # Prepare monthly outputs:
    ppfd_mo = (ppfd_mo*land_clip) + error_field
    pet_mo = (pet_mo*land_clip) + error_field
    eet_mo = (eet_mo*land_clip) + error_field
    cwd_mo = (pet_mo - aet_mo)*land_clip + error_field
    #
    # Save monthly values to file:
    if cur_date >= print_date:
        # Output file names:
        ppfd_out_file = "%s%s_%s_%s.txt" % (output_dir, out_no, 
                                            "PPFD_mo", cur_date)
        eet_out_file = "%s%s_%s_%s.txt" % (output_dir, out_no, 
                                           "EET_mo", cur_date)
        pet_out_file = "%s%s_%s_%s.txt" % (output_dir, out_no, 
                                           "PET_mo", cur_date)
        cpa_out_file = "%s%s_%s_%s.txt" % (output_dir, out_no, 
                                           "CPA_mo", cur_date)
        cwd_out_file = "%s%s_%s_%s.txt" % (output_dir, out_no, 
                                           "CWD_mo", cur_date)
        #
        save_to_file(ppfd_mo, ppfd_out_file)
        save_to_file(eet_mo, eet_out_file)
        save_to_file(pet_mo, pet_out_file)
        save_to_file(cpa_mo, cpa_out_file)
        save_to_file(cwd_mo, cwd_out_file)
    #
    # Update date field to next month:
    cur_date = add_one_month(cur_date)
#
if 0:
    qo_present = numpy.copy(qo_ann)
    qo_holocene = numpy.copy(qo_ann)
    #
    qo_out_file = "%s%s.txt" % (output_dir, "Qo_Annual_11000YA")
    qo_ann = (1e-3)*(qo_holocene*land_clip) + error_field
    save_to_file(qo_ann, qo_out_file)
