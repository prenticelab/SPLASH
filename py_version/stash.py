#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# stash.py
# * based on cramer_prentice.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-30 -- created
# 2014-11-25 -- last updated
#
# ------------
# description:
# ------------
# This script runs the STASH 2.0 model.
#
# ----------
# changelog:
# ----------
# 00. based on cramer_prentice.py created [14.01.30]
# 01. global constants [14.08.26]
# 02. EVAP class [14.08.26]
# 03. moved class constants to global constants
# 04. updated komega to float --- may influence cooper delta [14.09.29]
# 05. added 'berger' lamm method [14.09.30]
# 06. added check to woolf's method for lambda > 360 [14.10.01]
# 07. added Spencer method for declination [14.10.10]
# 08. replaced tau with Allen (1996); removed kZ [14.10.10]
# 09. distinguished shortwave from visible light albedo [14.10.16]
# 10. updated value and reference for semi-major axis, a [14.10.31]
# 11. fixed Cooper's and Spencer's declination equations [14.11.25]
# 12. replaced simplified kepler with full kepler [14.11.25]
#
# -----
# todo:
# -----
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import numpy

###############################################################################
## GLOBAL CONSTANTS:
###############################################################################
kA = 107       # constant for Rnl (Monteith & Unsworth, 1990)
ka = 1.49598e8 # semi-major axis, km (Allen, 1973)
kalb_sw = 0.17 # shortwave albedo (Federer, 1968)
kalb_vis = 0.03 # visible light albedo (Sellers, 1985)
kb = 0.20      # constant for Rnl (Linacre, 1968)
kc = 0.25      # cloudy transmittivity (Linacre, 1968)
kCw = 1.05     # supply constant, mm/hr (Federer, 1982)
kd = 0.50      # angular coefficient of transmittivity (Linacre, 1968)
ke = 0.0167    # eccentricity for 2000 CE (Berger, 1978)
keps = 23.44   # obliquity for 2000 CE, degrees (Berger, 1978)
kfFEC = 2.04   # from flux to energy conversion, umol/J (Meek et al., 1984)
kG = 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
kGM = 1.32712e11 # standard gravity of the sun, km^3/s^2
kGsc = 1360.8  # solar constant, W/m^2 (Kopp & Lean, 2011)
kL = 0.0065    # temperature lapse rate, K/m (Cavcar, 2000)
kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kMv = 0.01802  # molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
ksb = 5.670373e-8 # Stefan-Boltzman constant, W/m^2/K^4
kTo = 298.15   # base temperature, K (Prentice, unpublished)
kWm = 150      # soil moisture capacity, mm (Cramer & Prentice, 1988)
kw = 0.26      # entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
komega = 283.0 # longitude of perihelion for 2000 CE, degrees (Berger, 1978)

###############################################################################
## CLASSES
###############################################################################
class EVAP:
    """
    Name:     EVAP
    Features: This class calculates daily radiation and evapotranspiration
              quantities
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
              Cooper, P. I. (1969), The absorption of radiation in solar 
                stills, Solar Energy, vol. 12(3), pp. 333–346 
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
        self, lon, lat, n, elv=0.0, y=0, sf=1.0, tc=23.0, sw=1.0, 
        drm='loutre', lamm = 'kepler', delm = 'loutre'):
        """
        Name:     EVAP.__init__
        Input:    - float, longitude, degrees (lon)
                  - float, latitude, degrees (lat)
                  - int, day of the year (n)
                  - float, elevation, m (elv)
                  - int, year (y)
                  - float, fraction of sunshine hours (sf)
                  - float, mean daily air temperature, C (tc)
                  - float, evaporative supply rate, mm/hr (sw)
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
        # Error handle and assign required public variables:
        if lon > 180.0 or lon < -180.0:
            print "Longitude outside range of validity (-180 to 180)!"
            exit(1)
        else:
            self.user_lon = lon
            #
        if lat > 90.0 or lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.user_lat = lat
            #
        if n < 1 or n > 366:
            print "Day of year outside range of validity (1 to 366)!"
            exit(1)
        else:
            self.user_day = n
            #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.user_year == 0:
            self.kN = 365
        else:
            self.kN = (
                self.julian_day((self.user_year+1),1,1) - 
                self.julian_day(self.user_year, 1, 1)
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if lamm == 'kepler':
            # Kepler Method:
            self.my_nu, self.my_lambda = self.map_days(self.user_day, 
                                                       self.user_year)
        elif lamm == 'woolf':
            # Woolf Method:
            woolf_b = (self.user_day - 1.0)*(360.0/self.kN)
            self.my_lambda = (
                279.9348 + woolf_b + 1.914827*self.dsin(woolf_b) - 
                0.079525*self.dcos(woolf_b) + 0.019938*self.dsin(2*woolf_b) - 
                0.00162*self.dcos(2*woolf_b)
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
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if drm == 'loutre':
            # Eq. 1, Loutre (2002)
            my_rho = (1.0 - ke**2.0)/(1.0 + ke*self.dcos(self.my_nu))
            self.dr = (1.0/my_rho)**2.0
        elif drm == 'klein':
            # Eq. 11, STASH 2.0 Documentation (Klein, 1977)
            self.dr = (
                1.0 + 2.0*ke*numpy.cos(2.0*numpy.pi*self.user_day/self.kN)
            )
        else:
            print "Distance factor method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if delm == 'loutre':
            # Eq. 14, STASH 2.0 Documentation (Loutre, 2002)
            self.delta = numpy.arcsin(self.dsin(self.my_lambda)*self.dsin(keps))
            self.delta *= (180.0/numpy.pi)
        elif delm == 'cooper':
            # Eq. 16, STASH 2.0 Documentation (Cooper, 1969)
            self.delta = keps*numpy.sin(
                2.0*numpy.pi*(self.user_day + 284.0)/self.kN
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
                          0.002697*numpy.cos(3.0*spencer_b) + 
                          0.00148*numpy.sin(3.0*spencer_b))
            self.delta *= (180.0/numpy.pi)
        else:
            print "Declination angle method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(self.delta)*self.dsin(self.user_lat)
        rv = self.dcos(self.delta)*self.dcos(self.user_lat)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Note: u/v == tan(delta)*tan(lat)
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            # Polar day (no sunset)
            self.hs = 180.0 
        elif (ru/rv) <= -1.0:
            # Polar night (no sunrise)
            self.hs = 0.0
        else:
            self.hs = numpy.arccos(-1.0*ru/rv)
            self.hs *= (180.0/numpy.pi)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 1.10.3, Duffy & Beckman (1993)
        self.ra_d = (86400.0/numpy.pi)*kGsc*self.dr*(
            ru*(numpy.pi/180.0)*self.hs + rv*self.dsin(self.hs)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate transmittivity (tau), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Linacre (1968)
        tau_o = (kc + kd*self.user_sf)
        #
        # Eq. 2, Allen (1996)
        tau = tau_o*(1.0 + (2.67e-5)*self.user_elv)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 57, STASH 2.0 Documentation
        self.ppfd_d = (1.0e-6)*kfFEC*(1.0-kalb_vis)*tau*self.ra_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation (rnl), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Prentice et al. (1993);
        # Eq. 5 and 6, Linacre (1968)
        self.rnl = (
            (kb + (1.0 - kb)*self.user_sf)*(kA - self.user_tc)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate variable substitute (rw), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = (1.0-kalb_sw)*tau*kGsc*self.dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate net radiation cross-over hour angle (hn), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (self.rnl - rw*ru)/(rw*rv) >= 1.0:
            # Net radiation negative all day
            self.hn = 0
        elif (self.rnl - rw*ru)/(rw*rv) <= -1.0:
            # Net radiation positive all day
            self.hn = 180.0
        else:
            self.hn = (180.0/numpy.pi)*numpy.arccos((self.rnl - rw*ru)/(rw*rv))
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate daytime net radiation (rn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 53, STASH 2.0 Documentation
        self.rn_d = (86400.0/numpy.pi)*(self.hn*(numpy.pi/180.0)*(
            rw*ru - self.rnl) + rw*rv*self.dsin(self.hn)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate nighttime net radiation (rnn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 56, STASH 2.0 Documentation
        rnn_d = (86400.0/numpy.pi)*(
            rw*ru*(self.hs-self.hn)*(numpy.pi/180.0) + 
            rw*rv*(self.dsin(self.hs)-self.dsin(self.hn)) + self.rnl*(
                numpy.pi - 2.0*self.hs*(numpy.pi/180.0) + self.hn*(numpy.pi/180.0)
            )
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 15. Calculate water-to-energy conversion (econ), m^3/J
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
        self.econ = s/(lv*pw*(s + g))
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 16. Calculate daily condensation (wc), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 68, STASH 2.0 Documentation
        self.wc = 1000.0*self.econ*numpy.abs(rnn_d)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 17. Estimate daily EET (eet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 70, STASH 2.0 Documentation
        self.eet_d = 1000.0*self.econ*(self.rn_d)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 18. Estimate daily PET (pet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 72, STASH 2.0 Documentation
        self.pet_d = (1.0+kw)*self.eet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = 1000.0*3600.0*(1.0+kw)*self.econ
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 20. Calculate the intersection hour angle (hi), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cos_hi = (
            self.user_sw/(rw*rv*rx) + self.rnl/(rw*rv) - ru/rv
        )
        if cos_hi >= 1.0:
            # Supply exceeds demand:
            self.hi = 0.0
        elif cos_hi <= -1.0:
            # Supply limits demand everywhere:
            self.hi = 180.0
        else:
            self.hi = numpy.arccos(cos_hi)*(180.0/numpy.pi)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 21. Estimate daily AET (aet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 81, STASH 2.0 Documentation
        self.aet_d = (24.0/numpy.pi)*(
            self.user_sw*self.hi*(numpy.pi/180.0) + 
            rx*rw*rv*(self.dsin(self.hn) - self.dsin(self.hi)) + 
            (rx*rw*ru - rx*self.rnl)*(self.hn - self.hi)*(numpy.pi/180.0)
        )
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def dcos(self, x):
        """
        Name:     EVAP.dcos
        Input:    float, angle, degrees (x)
        Output:   float, cos(x*pi/180)
        Features: Calculates the cosine of an angle given in degrees
        """
        return numpy.cos(x*numpy.pi/180.0)
    #
    def dsin(self, x):
        """
        Name:     EVAP.dsin
        Input:    float, angle, degrees (x)
        Output:   float, sin(x*pi/180)
        Features: Calculates the sine of an angle given in degrees
        """
        return numpy.sin(x*numpy.pi/180.0)
    #
    def dtan(self, x):
        """
        Name:     EVAP.dtan
        Input:    float, angle, degrees (x)
        Output:   float, tan(x*pi/180)
        Features: Calculates the tangent of an angle given in degrees
        """
        return numpy.tan(x*numpy.pi/180.0)
    #
    def earth_velocity(self, lon):
        """
        Name:     EVAP.earth_velocity
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
                1.0+ke*self.dcos(lon)
                )**2.0*(1.0-ke**2)**(-1.5)
        )
        return w
    #
    def julian_day(self,y,m,i):
        """
        Name:     EVAP.julian_day
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
        jde = int(365.25*(y+4716)) + int(30.6001*(m+1)) + i + b - 1524.5
        return jde
    #
    def equinox(self, year, opt=0):
        """
        Name:     EVAP.equinox
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
            y = (year-2000.0)/1000.0
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
    def full_kepler(self, lon):
        """
        Name:     EVAP.full_kepler
        Input:    float, longitude w.r.t. the perihelion, deg (lon)
        Output:   float, days traveled from perihelion (t)
        Features: Returns the days traveled since the earth past the perihelion
        Depends:  -ke (ellipticity)
        Ref:      Kepler's Second Law
        """
        # Make lon into numpy array, if it is not already:
        if isinstance(lon, numpy.float) or isinstance(lon, numpy.int):
            lon = numpy.array([lon,])
        elif not isinstance(lon, numpy.ndarray):
            lon = numpy.array(lon)
        #
        xee = ke**2
        xte = ke**3
        xse = numpy.sqrt(1.0 - xee)
        xco = self.dcos(lon) + 1
        #
        A = (0.5*self.kN/numpy.pi)*(1 - xee)**1.5
        B = 2*numpy.arctan((ke - 1.)*self.dsin(lon)/(xse*xco))
        C = xse*(xee - 1)
        D = 2.*ke*self.dsin(lon)
        E = xco*((xte - xee - ke + 1)*self.dsin(lon)**2/xco**2 - 
                 xte - xee + ke + 1)
        kepler_t = A*(B/C - D/E)
        #
        # Correct negative points due to inverse tan function:
        neg_points = numpy.where(lon > 180)[0]
        kepler_t[neg_points] -= 2*numpy.pi*A/C
        #
        return (kepler_t)
    #
    def map_days(self, n, y):
        """
        Name:     EVAP.map_days
        Input:    - int, day of the year (n)
                  - int, year (y)
        Output:   float, longitude relative to the vernal equinox (lamda_doy)
        Features: Computes earth's longitude relative to the vernal equinox 
                  for a given day of the year
        Depends:  Functions:
                  - earth_period
                  - earth_velocity
                  - equinox
                  - full_kepler
                  Global constants:
                  - komega
        """
        # Create longitude field and compute the days for orbit:
        lon_nu = numpy.array([360.*i/365. for i in xrange(366)])
        day_nu = self.full_kepler(lon_nu)
        #
        # Compute the angle of the vernal equinox w.r.t. the perihelion
        # i.e., the explementary angle of komega:
        wp = (360.0 - komega)
        #
        # Calculate the length of time it takes earth to travel from the
        # perihelion (t=0) to the vernal equinox:
        if wp == 180:
            wp += 1e-6
        days_to_ve = self.full_kepler(wp)[0]
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
        if nu_doy >= 360:
            nu_doy -= 360
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
        Name:     EVAP.sat_slope
        Input:    float, air temperature (tc), degrees C
        Output:   float, slope of sat vap press temp curve (s)
        Features: Calculates the slope of the sat pressure temp curve, Pa/K
        Ref:      Eq. 13, Allen et al. (1998)
        """
        s = (17.269)*(237.3)*(610.78)*(
            numpy.exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2)
        )
        return s
    #
    def enthalpy_vap(self, tc):
        """
        Name:     EVAP.enthalpy_vap
        Input:    float, air temperature (tc), degrees C
        Output:   float, latent heat of vaporization
        Features: Calculates the enthalpy of vaporization, J/kg
        Ref:      Eq. 8, Henderson-Sellers (1984)
        """
        return (1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2.0)
    #
    def elv2pres(self, z):
        """
        Name:     EVAP.elv2pres
        Input:    float, elevation above sea level (z), m
        Output:   float, atmospheric pressure, Pa
        Features: Calculates atm. pressure for a givene elevation
        Depends:  Global constants
                  - kPo
                  - kTo
                  - kL
                  - kMa
                  - kG
                  - kR
        Ref:      Cavcar (2000)
        """
        p = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
        return p
    #
    def density_h2o(self, tc, p):
        """
        Name:     EVAP.density_h2o
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
            1000.0*po*(ko + ca*pbar + cb*pbar**2.0)/(
                ko + ca*pbar + cb*pbar**2.0 - pbar
            )
        )
        return pw
    #
    def psychro(self, tc, p):
        """
        Name:     EVAP.psychro
        Input:    - float, air temperature (tc), degrees C
                  - float, atm. pressure (p), Pa
        Output:   float, psychrometric constant, Pa/K
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
def julian_day(y, m, i):
    """
    Name:     julian_day
    Input:    - int, year (y)
              - int, month (m)
              - int, day of month (i)
    Output:   float, Julian Ephemeris Day
    Features: Converts Gregorian date (year, month, day) to Julian Ephemeris
              Day
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

###############################################################################
## MAIN PROGRAM 
###############################################################################
if 0:
    my_evap = EVAP(
        lon = -0.641,
        lat = 51.4,
        n = 172,
        elv = 74.0,
        y = 2001,
        sf = 0.43,
        tc = 17.3,
        sw = 0.5,
        drm = 'loutre',  # loutre or klein
        lamm = 'kepler',  # kepler, woolf or berger
        delm = 'loutre'  # loutre, spencer, cooper or circle
    )

# Example data (Met Office average climate)
data = numpy.array([
    (0.21, 4.80, 61.0),
    (0.27, 4.85, 41.2),
    (0.30, 7.10, 44.5), 
    (0.40, 9.10, 48.0), 
    (0.39, 12.4, 46.4), 
    (0.39, 15.3, 44.6), 
    (0.40, 17.6, 46.0), 
    (0.43, 17.3, 52.3), 
    (0.36, 14.6, 50.3), 
    (0.32, 11.2, 71.8), 
    (0.23, 7.55, 66.3),
    (0.19, 5.05, 62.9)
], dtype=[('sf',numpy.float),('tc',numpy.float),('ppt',numpy.float)])

# User definitions:
my_lon = -0.641 # longitude, degrees
my_lat = 51.4   # latitude, degrees
my_elv = 74.0   # elevation, m

# Initialize soil moisture
w = numpy.array([0.0*j for j in xrange(366)])

# Create empty daily results
daily_results = numpy.empty(
    shape=(366,), 
    dtype=[
        ('ra_d',numpy.float), 
        ('rn_d', numpy.float), 
        ('ppfd_d', numpy.float), 
        ('wc', numpy.float), 
        ('eet_d', numpy.float), 
        ('pet_d', numpy.float), 
        ('aet_d', numpy.float)
    ]
)

# Create empty monthly totals:
monthly_results = numpy.zeros(
    shape=(12,),
    dtype=[
        ('EET', numpy.float),
        ('PET', numpy.float),
        ('AET', numpy.float),
        ('CP', numpy.float),
        ('CWD', numpy.float),
        ('PPFD', numpy.float),
        ('RO', numpy.float)
    ]
)

# Define year and days in year:
y = 2000
ny = julian_day(y+1, 1, 1) - julian_day(y, 1, 1)

# Iterate through the months:
months = [j+1 for j in xrange(12)]
for m in months:
    nm = julian_day(y, m+1, 1) - julian_day(y, m, 1)
    days = [j+1 for j in xrange(int(nm))]
    for i in days:
        # Calculate the day of the year:
        n = julian_day(y, m, i) - julian_day(y, 1, 1) + 1
        #
        # Get index for yesterday (Note: zero indexing)
        idx = int(n-1)-1
        if idx < 0:
            idx = int(ny-1)
        #
        # Calculate evaporative supply rate
        sw = kCw*w[idx]/kWm
        #
        # Calculate radiation and evaporation quantities
        my_evap = EVAP(
            lon = my_lon,
            lat = my_lat,
            n = n,
            elv = my_elv,
            y = y,
            sf = data['sf'][m-1],
            tc = data['tc'][m-1],
            sw = sw,
            drm = 'loutre',  # loutre or klein
            lamm = 'kepler',  # kepler, berger or woolf
            delm = 'loutre'  # loutre, spencer, cooper or circle
        )
        #
        # Update soil moisture:
        ro = w[idx] + data['ppt'][m-1]/nm + my_evap.wc - my_evap.aet_d
        if ro > kWm:
            # Bucket is full 
            # * set soil moisture to capacity
            # * add remaining water to monthly runoff total
            w[n-1] = kWm
            monthly_results['RO'][m-1] += (ro - kWm)
        elif ro < 0:
            # Bucket is empty
            # * set soil moisture to zero
            w[n-1] = 0
        else:
            w[n-1] = ro
        #
        # Save the daily results:
        daily_results['ra_d'][n-1] = my_evap.ra_d
        daily_results['rn_d'][n-1] = my_evap.rn_d
        daily_results['ppfd_d'][n-1] = my_evap.ppfd_d
        daily_results['wc'][n-1] = my_evap.wc
        daily_results['eet_d'][n-1] = my_evap.eet_d
        daily_results['pet_d'][n-1] = my_evap.pet_d
        daily_results['aet_d'][n-1] = my_evap.aet_d
        #
        # Accumulate the monthly totals:
        monthly_results['EET'][m-1] += my_evap.eet_d
        monthly_results['PET'][m-1] += my_evap.pet_d
        monthly_results['AET'][m-1] += my_evap.aet_d
        monthly_results['PPFD'][m-1] += my_evap.ppfd_d
    #
    # Calculate monthly quantities
    monthly_results['CP'][m-1] = (
        monthly_results['AET'][m-1]/monthly_results['PET'][m-1]
    )
    monthly_results['CWD'][m-1] = (
    	monthly_results['PET'][m-1] - monthly_results['AET'][m-1]
    )
#
