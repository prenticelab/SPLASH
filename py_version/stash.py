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
# 2015-02-07 -- last updated
#
# ------------
# description:
# ------------
# This script runs the STASH 2.0 model.
#
# IMPORTANT NOTE: 
#   Global variables are defined inside a function at definition; therefore, 
#   you must re-run function definitions if you change global variable values

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
# 17. removed longitude from EVAP & STASH classes [15.01.13]
# 18. general housekeeping [15.01.13]
# 19. updated plots for results [15.01.16]
# 20. added example data CSV file & updated data for daily input [15.01.16]
# 21. fixed spin_up indexing in STASH class [15.01.16]
# 22. fixed Cramer-Prentice alpha definition [15.01.16]
# 23. updated plots [15.01.18]
# 24. updated reference to kL [15.01.29]
# 25. general housekeeping on EVAP class [15.02.07]
# 25. changed condensation variable name from 'wc' to 'cn' [15.02.07]
#
# -----
# todo:
# -----
# 1. Add a check to make certain the actual daily evapotranspiration does not 
#    exceed inputs (i.e., precipitation and condensation) plus reserves (i.e., 
#    yesterday's soil moisture content).
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import matplotlib.pyplot as plt
import numpy

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
kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
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
              - daily condensation (cn), mm
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
    def __init__(self, lat, n, elv=0.0, y=0, sf=1.0, tc=23.0, sw=1.0):
        """
        Name:     EVAP.__init__
        Input:    - float, latitude, degrees (lat)
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
            kN = 365
            self.user_year = 2001
        else:
            kN = self.julian_day((y+1),1,1) - self.julian_day(y, 1, 1)
        self.kN = kN
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        self.my_nu = my_nu
        self.my_lambda = my_lambda
        # 
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        kee = ke**2
        my_rho = (1.0 - kee)/(1.0 + ke*self.dcos(my_nu))
        dr = (1.0/my_rho)**2
        self.dr = dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        pir = (numpy.pi/180.0)
        delta = numpy.arcsin(self.dsin(my_lambda)*self.dsin(keps))
        delta /= pir
        self.delta = delta
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(delta)*self.dsin(lat)
        rv = self.dcos(delta)*self.dcos(lat)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            # Polar day (no sunset)
            hs = 180.0 
        elif (ru/rv) <= -1.0:
            # Polar night (no sunrise)
            hs = 0.0
        else:
            hs = -1.0*ru/rv
            hs = numpy.arccos(hs)
            hs /= pir
        self.hs = hs
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 1.10.3, Duffy & Beckman (1993)
        ra_d = (86400.0/numpy.pi)*kGsc*dr*(ru*pir*hs + rv*self.dsin(hs))
        self.ra_d = ra_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate transmittivity (tau), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
        tau_o = (kc + kd*sf)
        tau = tau_o*(1.0 + (2.67e-5)*elv)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra_d
        self.ppfd_d = ppfd_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation (rnl), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        rnl = (kb + (1.0 - kb)*sf)*(kA - tc)
        self.rnl = rnl
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate variable substitute (rw), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = (1.0 - kalb_sw)*tau*kGsc*dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate net radiation cross-over hour angle (hn), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rnl - rw*ru)/(rw*rv) >= 1.0:
            # Net radiation negative all day
            hn = 0
        elif (rnl - rw*ru)/(rw*rv) <= -1.0:
            # Net radiation positive all day
            hn = 180.0
        else:
            hn = (rnl - rw*ru)/(rw*rv)
            hn = numpy.arccos(hn)
            hn /= pir
        self.hn = hn
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate daytime net radiation (rn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rn_d = (86400.0/numpy.pi)*(hn*pir*(rw*ru - rnl) + rw*rv*self.dsin(hn))
        self.rn_d = rn_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate nighttime net radiation (rnn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rnn_d = rw*ru*(hs - hn)*pir
        rnn_d += rw*rv*(self.dsin(hs) - self.dsin(hn))
        rnn_d += rnl*(numpy.pi - 2.0*hs*pir + hn*pir)
        rnn_d *= (86400.0/numpy.pi)
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
        econ = s/(lv*pw*(s + g))
        self.econ = econ
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 16. Calculate daily condensation (cn), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cn = (1e3)*econ*numpy.abs(rnn_d)
        self.cn = cn
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 17. Estimate daily EET (eet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eet_d = (1e3)*econ*rn_d
        self.eet_d = eet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 18. Estimate daily PET (pet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pet_d = (1.0 + kw)*eet_d
        self.pet_d = pet_d
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = (3.6e6)*(1.0 + kw)*econ
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 20. Calculate the intersection hour angle (hi), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cos_hi = sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
        if cos_hi >= 1.0:
            # Supply exceeds demand:
            hi = 0.0
        elif cos_hi <= -1.0:
            # Supply limits demand everywhere:
            hi = 180.0
        else:
            hi = numpy.arccos(cos_hi)
            hi /= pir
        self.hi = hi
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 21. Estimate daily AET (aet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        aet_d = sw*hi*pir
        aet_d += rx*rw*rv*(self.dsin(hn) - self.dsin(hi))
        aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*pir
        aet_d *= (24.0/numpy.pi)
        self.aet_d = aet_d
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
        pir = numpy.pi/180.0
        #
        # Mean longitude for vernal equinox:
        xlam =(ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(komega)
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
        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
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
        po = 0.99983952
        po += (6.788260e-5)*tc
        po += -(9.08659e-6)*tc*tc
        po += (1.022130e-7)*tc*tc*tc
        po += -(1.35439e-9)*tc*tc*tc*tc
        po += (1.471150e-11)*tc*tc*tc*tc*tc
        po += -(1.11663e-13)*tc*tc*tc*tc*tc*tc
        po += (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc
        po += -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
        #
        # Calculate bulk modulus at 1 atm:
        ko = 19652.17
        ko += 148.1830*tc
        ko += -2.29995*tc*tc
        ko += 0.01281*tc*tc*tc
        ko += -(4.91564e-5)*tc*tc*tc*tc
        ko += (1.035530e-7)*tc*tc*tc*tc*tc
        #
        # Calculate temperature dependent coefficients:
        ca = 3.26138
        ca += (5.223e-4)*tc
        ca += (1.324e-4)*tc*tc
        ca += -(7.655e-7)*tc*tc*tc
        ca += (8.584e-10)*tc*tc*tc*tc
        #
        cb = (7.2061e-5)
        cb += -(5.8948e-6)*tc
        cb += (8.69900e-8)*tc*tc
        cb += -(1.0100e-9)*tc*tc*tc
        cb += (4.3220e-12)*tc*tc*tc*tc
        #
        # Convert atmospheric pressure to bar (1 bar = 100000 Pa)
        pbar = (1.0e-5)*p
        #
        pw = (ko + ca*pbar + cb*pbar**2.0)
        pw /= (ko + ca*pbar + cb*pbar**2.0 - pbar)
        pw *= (1e3)*po
        #
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
        cp = 1.0045714270
        cp += (2.050632750e-3)*tc
        cp += -(1.631537093e-4)*tc*tc
        cp += (6.212300300e-6)*tc*tc*tc
        cp += -(8.830478888e-8)*tc*tc*tc*tc
        cp += (5.071307038e-10)*tc*tc*tc*tc*tc
        cp *= (1e3)
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
    def __init__(self, lat, elv):
        """
        Name:     STASH.__init__
        Input:    - float, latitude, degrees (lat)
                  - float, elevation, meters (elv)
        """
        # Error handle and assign required public variables:
        self.elv = elv
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
    def spin_up(self, y, ppt, tc, sf):
        """
        Name:     STASH.spin
        Input:    - int, year (y)
                  - numpy.ndarray, daily precip, mm (ppt) 
                  - numpy.ndarray, daily air temperature, deg C (tc) 
                  - numpy.ndarray, daily sunshine fraction (sf)
        Output:   None.
        Features: Spins up the daily soil moisture
        """
        # Determine end-of-year index:
        ny = self.julian_day(y+1, 1, 1) - self.julian_day(y, 1, 1)
        idx = int(ny-1)
        #
        # Run once:
        self.run_one_year(y, ppt, tc, sf)
        #
        # Check to see if 1 Jan soil moisture matches 31 Dec:
        start_sm = self.daily_totals['wn'][0]
        end_sm = self.daily_totals['wn'][idx]
        diff_sm = numpy.abs(end_sm - start_sm)
        #
        # If not, spin again!
        self.spin_count = 1
        while (diff_sm > 1):
            self.run_one_year(y, ppt, tc, sf)
            start_sm = self.daily_totals['wn'][0]
            end_sm = self.daily_totals['wn'][idx]
            diff_sm = numpy.abs(end_sm - start_sm)
            self.spin_count += 1

        # Once in equilibrium, write daily and monthly values 
        # of one year to file
        #self.write_to_file()

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
                # Note that 'n' represents day of year and runs from 1 to 365 
                # (366), whereas idx is the index and runs from 0 to 364 (365) 
                # (zero-indexing in Python!)
                idx = int(n - 1) - 1
                if idx < 0:
                    idx = int(ny - 1)
                #
                # Calculate evaporative supply rate, mm/h
                sw = kCw*self.daily_totals['wn'][idx]/kWm
                #
                # Calculate radiation and evaporation quantities
                my_evap = EVAP(self.lat, n, self.elv, y, sf[n-1], tc[n-1], sw)
                #
                # Calculate daily precipitation:
                self.daily_totals['pn'][n-1] = ppt[n-1]
                #
                # Update soil moisture:
                self.daily_totals['wn'][n-1] = (self.daily_totals['wn'][idx] +
                    self.daily_totals['pn'][n-1] + my_evap.cn - my_evap.aet_d)
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
                self.daily_totals['cn'][n-1] = my_evap.cn
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
            self.monthly_totals['cpa'][m-1] /= self.monthly_totals['eq_m'][m-1]
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
    my_evap = EVAP(lat = 51.4,
                   n = 172,
                   elv = 74.0,
                   y = 2001,
                   sf = 0.43,
                   tc = 17.3,
                   sw = 0.5)

# Example data (San Francisco, 2000 CE)
my_file = 'example_data.csv'
data = numpy.loadtxt(my_file, 
                     dtype={'names': ('sf', 'tair', 'pn'),
                            'formats' : ('f4', 'f4', 'f4')},
                     delimiter=',',
                     skiprows=1)
my_lat = 37.7   # latitude, degrees
my_elv = 142.   # elevation, m

# Create STASH class and spin up soil moisture:
my_class = STASH(my_lat, my_elv)
my_class.spin_up(2000, data['pn'], data['tair'], data['sf'])
#my_class.run_one_year(2000, data['ppt'], data['tc'], data['sf'])

#
# Plot monthly ET results
#
my_months = range(1, 13)
fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
ax1.plot(my_months, my_class.monthly_totals['ep_m'], 
         'k-', linewidth=2, label='PET')
ax1.plot(my_months, my_class.monthly_totals['ea_m'], 
         'c--', linewidth=2, label='AET')
ax1.plot(my_months, my_class.monthly_totals['eq_m'], 
         'g-', linewidth=2, label='EET')
ax1.plot(my_months, my_class.monthly_totals['cwd'], 
         'r:', linewidth=2, label='CWD')
ax1.set_xticks(my_months)
ax1.set_xlabel('Months', fontsize=18)
ax1.set_ylabel('ET (mm)', fontsize=18)
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0., fontsize=18)
plt.xlim([1, 12])
plt.show()


#
# Plot monthly results
#
my_months = range(1, 13)
my_y_max = my_class.monthly_totals['ep_m'].max()
fig = plt.figure()
# [1]
ax1 = fig.add_subplot(411)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.plot(my_months, my_class.monthly_totals['ep_m'], 'k-', linewidth=2, 
         label='$E^p$')
ax1.plot(my_months, my_class.monthly_totals['ea_m'], 'k--', linewidth=2,
         label='$E^a$')
ax1.set_ylabel('$E_m$ (mm)', fontsize=14)
ax1.set_xticks(range(1, 13, 1))
ax1.set_yticks(range(0, 176, 50))
plt.xlim([1, 12])
plt.ylim([0, my_y_max])
ax1.text(1.25, 150, '(a)', fontsize=12)
g1 = ax1.legend(loc=1, fontsize=14)
f1 = g1.get_frame()
f1.set_linewidth(0)
# [2]
ax2 = fig.add_subplot(412)
plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.plot(my_months, my_class.monthly_totals['cwd'], 'k-', linewidth=2)
ax2.set_ylabel('$CWD$ (mm)', fontsize=14)
ax2.set_xticks(range(1, 13, 1))
ax2.set_yticks(range(0, 176, 50))
plt.xlim([1, 12])
plt.ylim([0, my_y_max])
ax2.text(1.25, 150, '(b)', fontsize=12)
# [3]
ax3 = fig.add_subplot(413)
plt.setp(ax3.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.plot(my_months, my_class.monthly_totals['eq_m'], 'k-', linewidth=2,
         label='$E^q$')
ax3.plot(my_months, my_class.monthly_totals['ea_m'], 'k--', linewidth=2,
         label='$E^a$')
ax3.set_ylabel('$E_m$ (mm)', fontsize=14)
ax3.set_xticks(range(1, 13, 1))
ax3.set_yticks(range(0, 176, 50))
plt.xlim([1, 12])
plt.ylim([0, my_y_max])
ax3.text(1.25, 150, '(c)', fontsize=12)
g3 = ax3.legend(loc=1, fontsize=14)
f3 = g3.get_frame()
f3.set_linewidth(0)
# [4]
ax4 = fig.add_subplot(414)
plt.setp(ax4.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax4.get_xticklabels(), rotation=0, fontsize=12)
ax4.plot(my_months, my_class.monthly_totals['cpa'], 'k-', linewidth=2)
ax4.set_ylabel('$\\alpha$', fontsize=14)
ax4.set_xticks(range(1, 13, 1))
ax4.set_yticks([0.3*i for i in range(0, 5, 1)])
plt.xlim([1, 12])
plt.ylim([0, 1.3])
ax4.text(1.25, 1.1, '(d)', fontsize=12)
plt.show()


#
# Plot daily ET results
#
my_days = range(1, 367)
fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=16)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=16)
ax1.plot(my_days, my_class.daily_totals['ep_n'], 
         'b-', linewidth=2, label='Potential')
ax1.plot(my_days, my_class.daily_totals['eq_n'], 
         'g-', linewidth=2, label='Equilibrium')
ax1.plot(my_days, my_class.daily_totals['ea_n'], 
         'r--', linewidth=2, label='Actual')
ax1.set_ylabel('Evapotranspiration, mm d$^{-1}$', fontsize=18)
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0., fontsize=18)
plt.show()


#
# Plot daily results
#
my_days = range(1, 367)
my_xtks = range(0, 370, 60)
fig = plt.figure()
# [1]
ax1 = fig.add_subplot(811)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.plot(my_days, data['sf'], 'k-', linewidth=2)
ax1.set_ylabel('$S_f$', fontsize=14)
ax1.set_xticks(my_xtks)
ax1.set_yticks([0.1*i for i in range(4, 8, 1)])
plt.xlim([-15, 380])
ax1.text(-5, 0.7, '(a)', fontsize=12)
# [2]
ax2 = fig.add_subplot(812)
plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.plot(my_days, (1e-6)*my_class.daily_totals['hn'], 'k-', linewidth=2)
ax2.set_ylabel('$H_N$ (MJ m$^{-2}$)', fontsize=14)
ax2.set_xticks(my_xtks)
ax2.set_yticks(range(6, 19, 3))
plt.xlim([-15, 380])
ax2.text(-5, 16, '(b)', fontsize=12)
# [3]
ax3 = fig.add_subplot(813)
plt.setp(ax3.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.plot(my_days, my_class.daily_totals['cn'], 'k-', linewidth=2)
ax3.set_ylabel('$C_n$ (mm)', fontsize=14)
ax3.set_xticks(my_xtks)
ax3.set_yticks([0.1*i for i in range(5, 9, 1)])
plt.xlim([-15, 380])
ax3.text(-5, 0.75, '(c)', fontsize=12)
# [4]
ax4 = fig.add_subplot(814)
plt.setp(ax4.get_yticklabels(), rotation=0, fontsize=12)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.plot(my_days, data['ppt'], 'k-', linewidth=2)
ax4.set_ylabel('$P_n$ (mm)', fontsize=14)
ax4.set_xticks(my_xtks)
ax4.set_yticks(range(0, 26, 5))
plt.xlim([-15, 380])
ax4.text(-5, 20, '(d)', fontsize=12)
# [5]
ax5 = fig.add_subplot(815)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), rotation=0, fontsize=12)
ax5.plot(my_days, my_class.daily_totals['wn'], 'k-', linewidth=2)
ax5.set_ylabel('$W_n$ (mm)', fontsize=14)
ax5.set_xticks(my_xtks)
ax5.set_yticks(range(30, 151, 30))
plt.xlim([-15, 380])
ax5.text(-5, 120, '(e)', fontsize=12)
# [6]
ax6 = fig.add_subplot(816)
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), rotation=0, fontsize=12)
ax6.plot(my_days, my_class.daily_totals['ro'], 'k-', linewidth=2)
ax6.set_ylabel('$RO$ (mm)', fontsize=14)
ax6.set_xticks(my_xtks)
ax6.set_yticks(range(0, 21, 5))
plt.xlim([-15, 380])
ax6.text(-5, 18, '(f)', fontsize=12)
# [7]
ax7 = fig.add_subplot(817)
plt.setp(ax7.get_xticklabels(), visible=False)
plt.setp(ax7.get_yticklabels(), rotation=0, fontsize=12)
ax7.plot(my_days, data['tc'], 'k-', linewidth=2)
ax7.set_ylabel('$T_{air}$ ($^{\circ}C$)', fontsize=14)
ax7.set_xticks(my_xtks)
ax7.set_yticks(range(10, 26, 5))
plt.xlim([-15, 380])
ax7.text(-5, 23, '(g)', fontsize=12)
# [8]
ax8 = fig.add_subplot(818)
plt.setp(ax8.get_xticklabels(), rotation=0, fontsize=12)
plt.setp(ax8.get_yticklabels(), rotation=0, fontsize=12)
ax8.plot(my_days, my_class.daily_totals['ep_n'], 'k-', linewidth=2)
ax8.plot(my_days, my_class.daily_totals['ea_n'], 'k--', linewidth=2)
ax8.set_ylabel('$E_n$ (mm)', fontsize=14)
ax8.set_xticks(my_xtks)
ax8.set_yticks([1.5*i for i in range(5)])
plt.xlim([-15, 380])
ax8.text(-5, 5, '(h)', fontsize=12)

plt.show()
