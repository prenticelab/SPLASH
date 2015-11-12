#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# evap.py
#
# 2014-01-30 -- created
# 2015-11-11 -- last updated
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
# led algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
# transpiration and plant-available moisture, Geoscientific Model Development,
# 2015 (in progress)

###############################################################################
## IMPORT MODULES:
###############################################################################
import numpy
from const import (ke, keps, kGsc, kA, kb, kc, kd, kfFEC, kalb_vis, kalb_sw,
                   kw, komega, kG, kL, kR, kPo, kTo, kMa, kMv)


###############################################################################
## CLASSES
###############################################################################
class EVAP:
    """
    Name:     EVAP
    Features: This class calculates daily radiation and evapotranspiration
              quantities:
              - PPFD (ppfd_d), mol/m^2/day
              - EET (eet_d), mm/day
              - PET (pet_d), mm/day
              - AET (aet_d), mm/day
              - condensation (cn), mm/day
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
            kN = self.julian_day((y+1), 1, 1) - self.julian_day(y, 1, 1)
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
        self.cond = cn
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

    def dsin(self, x):
        """
        Name:     EVAP.dsin
        Input:    float, angle, degrees (x)
        Output:   float, sin(x*pi/180)
        Features: Calculates the sine of an angle given in degrees
        """
        return numpy.sin(x*numpy.pi/180.0)

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
        Ref:      Berger, A. L. (1978), Long term variations of daily
                  insolation and quaternary climatic changes, J. Atmos. Sci.,
                  35, 2362-2367.
        """
        # Variable substitutes:
        xee = ke**2
        xec = ke**3
        xse = numpy.sqrt(1.0 - xee)
        pir = numpy.pi/180.0
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

    def julian_day(self, y, m, i):
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

    def enthalpy_vap(self, tc):
        """
        Name:     EVAP.enthalpy_vap
        Input:    float, air temperature (tc), degrees C
        Output:   float, latent heat of vaporization
        Features: Calculates the enthalpy of vaporization, J/kg
        Ref:      Eq. 8, Henderson-Sellers (1984)
        """
        return (1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2)

    def elv2pres(self, z):
        """
        Name:     EVAP.elv2pres
        Input:    float, elevation above sea level (z), m
        Output:   float, atmospheric pressure, Pa
        Features: Calculates atm. pressure for a given elevation
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
