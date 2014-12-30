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
# 2014-12-09 -- last updated
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
# 13. removed options for approximation methods not considering variable 
#     orbital velocity (e.g. Spencer, Woolf, Klein, Cooper, and Circle 
#     methods) [14.12.09]
# 14. reduced the list of constants and EVAP class functions [14.12.09]
# 15. added matplotlib to module list [14.12.09]
# 16. added plots for results [14.12.09]
#
# -----
# todo:
# -----
# 1. finish with STASH class
# 2. finish plots
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import matplotlib.pyplot as plt
import numpy
import csv

###############################################################################
## GLOBAL CONSTANTS:
###############################################################################
kA = 107       # constant for Rnl (Monteith & Unsworth, 1990)
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
kGsc = 1360.8  # solar constant, W/m^2 (Kopp & Lean, 2011)
kL = 0.0065    # temperature lapse rate, K/m (Cavcar, 2000)
kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kMv = 0.01802  # molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
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
    def __init__(self, lon, lat, n, elv=0.0, y=0, sf=1.0, tc=23.0, sw=1.0):
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
        if y == 0:
            self.kN = 365
            self.user_year = 2001
        else:
            self.kN = (
                self.julian_day((self.user_year+1),1,1) - 
                self.julian_day(self.user_year, 1, 1)
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        self.my_nu, self.my_lambda = self.berger_tls(n)
        # 
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        my_rho = (1.0 - ke**2)/(1.0 + ke*self.dcos(self.my_nu))
        self.dr = (1.0/my_rho)**2
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        self.delta = numpy.arcsin(self.dsin(self.my_lambda)*self.dsin(keps))
        self.delta *= (180.0/numpy.pi)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(self.delta)*self.dsin(lat)
        rv = self.dcos(self.delta)*self.dcos(lat)
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
        tau_o = (kc + kd*sf)
        #
        # Eq. 2, Allen (1996)
        tau = tau_o*(1.0 + (2.67e-5)*elv)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 57, STASH 2.0 Documentation
        self.ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*self.ra_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation (rnl), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        self.rnl = (kb + (1.0 - kb)*sf)*(kA - tc)
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
            rw*rv*(self.dsin(self.hs)-self.dsin(self.hn)) + 
            self.rnl*(numpy.pi - 2.0*self.hs*(numpy.pi/180.0) + 
                self.hn*(numpy.pi/180.0))
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 15. Calculate water-to-energy conversion (econ), m^3/J
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Slope of saturation vap press temp curve, Pa/K
        s = self.sat_slope(tc)
        # Enthalpy of vaporization, J/kg
        lv = self.enthalpy_vap(tc)
        # Density of water, kg/m^3
        pw = self.density_h2o(tc, self.elv2pres(elv))
        # Psychrometric constant, Pa/K
        g = self.psychro(tc, self.elv2pres(elv))
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
        cos_hi = sw/(rw*rv*rx) + self.rnl/(rw*rv) - ru/rv
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
    def berger_tls(self, n):
        """
        Name:     EVAP.berger_tls
        Input:    int, day of year
        Output:   tuple, 
                  - true anomaly, degrees
                  - true longitude, degrees
        Features: Returns true anomaly and true longitude for a given day
        Depends:  - ke
                  - komega
        Ref:      Berger, A. L. (1978), Long term variations of daily insolation
                  and quaternary climatic changes, J. Atmos. Sci., 35, 2362-
                  2367.
        """
        # Variable substitutes:
        xee = ke**2 
        xec = ke**3
        xse = numpy.sqrt(1.0 - xee)
        #
        # Mean longitude for vernal equinox:
        xlam = (
            (ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(komega) - 
            xee/4.0*(0.5 + xse)*self.dsin(2.0*komega) + 
            xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*komega)
            )
        xlam = numpy.degrees(2.0*xlam)
        #
        # Mean longitude for day of year:
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        #
        # Mean anomaly:
        anm = dlamm - komega
        ranm = numpy.radians(anm)
        #
        # True anomaly:
        ranv = (ranm + (2.0*ke - xec/4.0)*numpy.sin(ranm) + 
            5.0/4.0*xee*numpy.sin(2.0*ranm) + 
            13.0/12.0*xec*numpy.sin(3.0*ranm))
        anv = numpy.degrees(ranv)
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
        return (1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2)
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
        Ref:      Allen et al. (1998)
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
        Refs:     Allen et al. (1998); Tsilingiris (2008) 
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
class STASH:
    """
    Name:     STASH
    Features: This class maintains daily, monthly and annual quantities of 
              radiation, evapotranspiration, and soil moisture based on the 
              STASH methods
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization 
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lat, lon, elv):
        """
        Name:     STASH.__init__
        Input:    - float, latitude, degrees (lat)
                  - float, longitude, degrees (lon)
                  - float, elevation, meters (elv)
        """
        # Error handle and assign required public variables:
        self.elv = elv
        #
        if lon > 180.0 or lon < -180.0:
            print "Longitude outside range of validity (-180 to 180)!"
            exit(1)
        else:
            self.lon = lon
        #
        if lat > 90.0 or lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.lat = lat
        #
        # Initialize daily totals
        self.daily_totals = numpy.empty(
            shape=(366,), 
            dtype=[('ho',numpy.float),    # daily solar irradiation, J/m2
                   ('hn', numpy.float),   # daily net radiation, J/m2
                   ('qn', numpy.float),   # daily PPFD, mol/m2
                   ('cn', numpy.float),   # daily condensation water, mm
                   ('wn', numpy.float),   # daily soil moisture, mm
                   ('pn', numpy.float),   # daily precipitation, mm
                   ('ro', numpy.float),   # daily runoff, mm
                   ('eq_n', numpy.float), # daily equilibrium ET, mm
                   ('ep_n', numpy.float), # daily potential ET, mm
                   ('ea_n', numpy.float)] # daily actual ET, mm
        )
        #
        # Initialize monthly totals
        self.monthly_totals = numpy.zeros(shape=(12,),
                                          dtype=[('eq_m', numpy.float),
                                                 ('ep_m', numpy.float),
                                                 ('ea_m', numpy.float),
                                                 ('cpa', numpy.float),
                                                 ('cwd', numpy.float),
                                                 ('qm', numpy.float)])
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def spin_up(self, ppt, tc, sf):
        """
        Name:     STASH.spin
        Input:    - numpy.ndarray, monthly precip, mm (ppt) 
                  - numpy.ndarray, mean monthly temperature, deg C (tc) 
                  - numpy.ndarray, mean monthly sunshine fraction (sf)
        Output:   None.
        Features: Spins up the daily soil moisture
        """
        # Run once:
        self.run_one_year(2000, ppt, tc, sf)
        #
        # Check to see if 1 Jan soil moisture matches 31 Dec:
        start_sm = self.daily_totals['wn'][0]
        end_sm = self.daily_totals['wn'][-1]
        diff_sm = numpy.abs(end_sm - start_sm)
        #
        # If not, spin again!
        self.spin_count = 1
        while (diff_sm > 1e-4):
            self.run_one_year(2000, ppt, tc, sf)
            start_sm = self.daily_totals['wn'][0]
            end_sm = self.daily_totals['wn'][-1]
            diff_sm = numpy.abs(end_sm - start_sm)
            self.spin_count += 1

        # Once in equilibrium, write daily and monthly values 
        # of one year to file
        self.write_to_file()

    #
    def run_one_year(self, y, ppt, tc, sf):
        """
        Name:     STASH.run_one_year
        Input:    - int, year (y)
                  - numpy.ndarray, monthly total precip, mm (ppt)
                  - numpy.ndarray, mean monthly temp, deg C (tc)
                  - numpy.ndarray, mean fractional sunshine (sf)
        Output:   None.
        Features: Calculates daily and monthly quantities for one year
        """
        # Reset monthly totals:
        self.reset_monthly_totals()
        #
        # Calculates days in year:
        ny = self.julian_day(y+1, 1, 1) - self.julian_day(y, 1, 1)
        #
        # Iterate through the months:
        months = [j+1 for j in xrange(12)]
        for m in months:
            nm = self.julian_day(y, m+1, 1) - self.julian_day(y, m, 1)
            days = [j+1 for j in xrange(int(nm))]
            for i in days:
                # Calculate the day of the year:
                n = self.julian_day(y, m, i) - self.julian_day(y, 1, 1) + 1
                #
                # Get index for yesterday (Note: zero indexing)
                # xxx in effect, idx is TWO days before -- is this correct? xxx
                # Note that 'n' represents day of year and runs from 1 to 365 (366),
                # whereas idx is the index and runs from 0 to 364 (365) (zero-indexing 
                # in Python!)
                idx = int(n-1) - 1
                if idx < 0:
                    idx = int(ny - 1)
                #
                # Calculate evaporative supply rate, mm/h
                sw = kCw*self.daily_totals['wn'][idx]/kWm
                #
                # Calculate radiation and evaporation quantities
                my_evap = EVAP(self.lon, self.lat, n, self.elv, 2000, sf[m-1], 
                               tc[m-1], sw)
                #
                # Calculate daily precipitation:
                self.daily_totals['pn'][n-1] = ppt[m-1]/nm

                #
                # Update soil moisture:
                self.daily_totals['wn'][n-1] = (self.daily_totals['wn'][idx] +
                    self.daily_totals['pn'][n-1] + my_evap.wc - my_evap.aet_d)
                if self.daily_totals['wn'][n-1] > kWm:
                    # Bucket is full 
                    # * set soil moisture to capacity
                    # * add remaining water to monthly runoff total
                    self.daily_totals['ro'][n-1] = self.daily_totals['wn'][n-1]
                    self.daily_totals['ro'][n-1] -= kWm
                    self.daily_totals['wn'][n-1] = kWm
                elif self.daily_totals['wn'][n-1] < 0:
                    # Bucket is empty
                    # * set soil moisture to zero
                    self.daily_totals['wn'][n-1] = 0
                    self.daily_totals['ro'][n-1] = 0
                else:
                    self.daily_totals['ro'][n-1] = 0
                #
                # Save the daily totals:
                self.daily_totals['ho'][n-1] = my_evap.ra_d
                self.daily_totals['hn'][n-1] = my_evap.rn_d
                self.daily_totals['qn'][n-1] = my_evap.ppfd_d
                self.daily_totals['cn'][n-1] = my_evap.wc
                self.daily_totals['eq_n'][n-1] = my_evap.eet_d
                self.daily_totals['ep_n'][n-1] = my_evap.pet_d
                self.daily_totals['ea_n'][n-1] = my_evap.aet_d
                #
                # Update monthly totals:
                self.monthly_totals['eq_m'][m-1] += my_evap.eet_d
                self.monthly_totals['ep_m'][m-1] += my_evap.pet_d
                self.monthly_totals['ea_m'][m-1] += my_evap.aet_d
                self.monthly_totals['qm'][m-1] += my_evap.ppfd_d
                #
            # END LOOP ON DAYS
            # Calculate other monthly totals:
            self.monthly_totals['cpa'][m-1] = self.monthly_totals['ea_m'][m-1]
            self.monthly_totals['cpa'][m-1] /= self.monthly_totals['ep_m'][m-1]
            self.monthly_totals['cwd'][m-1] = self.monthly_totals['ep_m'][m-1]
            self.monthly_totals['cwd'][m-1] -= self.monthly_totals['ea_m'][m-1]
        # END LOOP ON MONTHS
    #
    def reset_monthly_totals(self):
        """
        Name:     STASH.reset_monthly_totals
        Input:    None.
        Output:   None.
        Features: Resets the monthly results to zero
        """
        self.monthly_totals = numpy.zeros(shape=(12,),
                                          dtype=[('eq_m', numpy.float),
                                                 ('ep_m', numpy.float),
                                                 ('ea_m', numpy.float),
                                                 ('cpa', numpy.float),
                                                 ('cwd', numpy.float),
                                                 ('qm', numpy.float)])
    #
    def write_to_file(self):
        """
        Name:     STASH.write_to_file
        Input:    None.
        Output:   None.
        Features: Writes daily and monthly values to files
        """
        #
        print "writing daily totals to files"
        numpy.savetxt('./output/ho.d.out', self.daily_totals['ho'])
        numpy.savetxt('./output/hn.d.out', self.daily_totals['hn'])
        numpy.savetxt('./output/qn.d.out', self.daily_totals['qn'])
        numpy.savetxt('./output/cn.d.out', self.daily_totals['cn'])
        numpy.savetxt('./output/eq_n.d.out', self.daily_totals['eq_n'])
        numpy.savetxt('./output/ep_n.d.out', self.daily_totals['ep_n'])
        numpy.savetxt('./output/ea_n.d.out', self.daily_totals['ea_n'])
        #
        print "writing monthly totals to files"
        numpy.savetxt('./output/eq_m.m.out',self.monthly_totals['eq_m'])
        numpy.savetxt('./output/ep_m.m.out',self.monthly_totals['ep_m'])
        numpy.savetxt('./output/ea_m.m.out',self.monthly_totals['ea_m'])
        numpy.savetxt('./output/qm.m.out',self.monthly_totals['qm'])
        numpy.savetxt('./output/cpa.m.out',self.monthly_totals['cpa'])
        numpy.savetxt('./output/cwd.m.out',self.monthly_totals['cwd'])
    #
    def julian_day(self, y, m, i):
        """
        Name:     STASH.julian_day
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
        sw = 0.5
        )

# Example data (Met Office average climate)
data = numpy.array([(0.21, 4.80, 61.0),
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
                    (0.19, 5.05, 62.9)], 
                   dtype=[('sf',numpy.float), 
                          ('tc',numpy.float), 
                          ('ppt',numpy.float)])

# User definitions:
my_lon = -0.641 # longitude, degrees
my_lat = 51.4   # latitude, degrees
my_elv = 74.0   # elevation, m
my_class = STASH(my_lat, my_lon, my_elv)
my_class.spin_up(data['ppt'], data['tc'], data['sf'])

# # Plot daily stuff:
# my_days = range(1,367)
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
# plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
# ax1.plot(my_days, my_class.daily_totals['wn'], 'g-', linewidth=2, 
#          label='Soil Moisture')
# ax1.set_xlabel('Day', fontsize=18)
# ax1.set_ylabel('Moisture, mm', fontsize=18)
# #ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
# #            ncol=4, mode="expand", borderaxespad=0., fontsize=18)
# plt.show()

# # Plot monthly ET
# my_months = range(1, 13)
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
# plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
# ax1.plot(my_months, my_class.monthly_totals['ep_m'], 'k-', linewidth=2, label='Potential')
# ax1.plot(my_months, my_class.monthly_totals['ea_m'], 'c--', linewidth=2, label='Actual')
# ax1.plot(my_months, my_class.monthly_totals['eq_m'], 'g-', linewidth=2, label='Equilibrium')
# ax1.plot(my_months, my_class.monthly_totals['cwd'], 'r:', linewidth=2, label='Deficit')
# ax1.set_xticks(my_months)
# ax1.set_xlabel('Months', fontsize=18)
# ax1.set_ylabel('Evapotranspiration, mm', fontsize=18)
# ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#             ncol=4, mode="expand", borderaxespad=0., fontsize=18)
# plt.xlim([1, 12])
# plt.show()

# # Plot monthly CPA
# my_months = range(1, 13)
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
# plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
# ax1.plot(my_months, my_class.monthly_totals['cpa'], 'k-', linewidth=2, 
#          label='Cramer-Prentice alpha')
# ax1.set_xticks(my_months)
# ax1.set_xlabel('Months', fontsize=18)
# ax1.set_ylabel('Evapotranspiration, mm', fontsize=18)
# ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#             ncol=1, mode="expand", borderaxespad=0., fontsize=18)
# plt.xlim([1, 12])
# plt.show()

# Plot monthly precip and runoff
#my_daily_precip = []
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
#plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
#ax1.plot(months, data['ppt'], 'b-', linewidth=2, label='Precip')
#ax1.plot(months, monthly_results['RO'], 'g--', linewidth=2, label='Runoff')
#ax1.plot(months, monthly_results['W'], 'c-', linewidth=2, label='Soil Moisture')
#ax1.set_xticks(months)
#ax1.set_xlabel('Months', fontsize=18)
#ax1.set_ylabel('Evapotranspiration, mm', fontsize=18)
#ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=3, mode="expand", borderaxespad=0., fontsize=18)
#plt.xlim([1, 12])
#plt.show()


# Plot daily radiation
#x_data = range(1,367)
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
#plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
#ax1.plot(x_data, daily_results['ra_d'], 'k-', linewidth=2, label='Ho')
#ax1.plot(x_data, daily_results['rn_d'], 'c--', linewidth=2, label='Hn')
#ax1.set_xlabel('Day', fontsize=18)
#ax1.set_ylabel('Radiation, J m$^{-2}$', fontsize=18)
#ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=4, mode="expand", borderaxespad=0., fontsize=18)
#plt.show()