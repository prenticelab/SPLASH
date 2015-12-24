#!/usr/bin/python
#
# solar.py
#
# 2015-12-22 -- last updated
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
# IMPORT MODULES:
###############################################################################
import logging
import logging.handlers
import sys

import numpy
from const import (ke, keps, kGsc, kA, kb, kc, kd, kfFEC, kalb_vis, kalb_sw,
                   komega, pir)


###############################################################################
# CLASSES
###############################################################################
class SOLAR:
    """
    Name:     SOLAR
    Features: This class calculates the daily radiation fluxes.
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lat, n, elv=0.0, y=0, sf=1.0, tc=23.0):
        """
        Name:     SOLAR.__init__
        Input:    - float, latitude, degrees (lat)
                  - int, day of the year (n)
                  - float, elevation, m (elv)
                  - int, year (y)
                  - float, fraction of sunshine hours (sf)
                  - float, mean daily air temperature, C (tc)
        Depends:  - julian_day
                  - berger_tls
                  - dcos
                  - dsin
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("SOLAR class called")

        # Assign default public variables:
        self.elv = elv
        self.logger.info("elevation set to %0.3f m", elv)

        # Error handle and assign required public variables:
        if lat > 90.0 or lat < -90.0:
            self.logger.error(
                "Latitude outside range of validity, (-90 to 90)!")
            sys.exit(1)
        else:
            self.logger.info("latitude set to %0.3f degrees", lat)
            self.lat = lat

        if n < 1 or n > 366:
            self.logger.error(
                "Day of year outside range of validity, (1 to 366)!")
            sys.exit(1)
        else:
            self.logger.info("day of year set to %d", n)
            self.day = n

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            kN = 365
            self.year = 2001
            self.logger.info("year set to %d", 2001)
        else:
            kN = self.julian_day((y+1), 1, 1) - self.julian_day(y, 1, 1)
            self.year = y
            self.logger.info("year set to %d", y)
        self.kN = kN
        self.logger.info("number of days in year set to %d", kN)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        self.my_nu = my_nu
        self.my_lambda = my_lambda
        self.logger.info("nu set to %f degrees", my_nu)
        self.logger.info("lambda set to %f degrees", my_lambda)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        kee = ke**2
        my_rho = (1.0 - kee)/(1.0 + ke*self.dcos(my_nu))
        dr = (1.0/my_rho)**2
        self.dr = dr
        self.logger.info("rho set to %f", my_rho)
        self.logger.info("distance factor, dr, set to %f", dr)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        delta = numpy.arcsin(self.dsin(my_lambda)*self.dsin(keps))
        delta /= pir
        self.delta = delta
        self.logger.info("declination, delta, set to %f", delta)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(delta)*self.dsin(lat)
        rv = self.dcos(delta)*self.dcos(lat)
        self.ru = ru
        self.rv = rv
        self.logger.info("variable substitute, ru, set to %f", ru)
        self.logger.info("variable substitute, rv, set to %f", rv)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            # Polar day (no sunset)
            self.logger.debug("polar day---no sunset")
            hs = 180.0
        elif (ru/rv) <= -1.0:
            # Polar night (no sunrise)
            self.logger.debug("polar night---no sunrise")
            hs = 0.0
        else:
            hs = -1.0*ru/rv
            hs = numpy.arccos(hs)
            hs /= pir
        self.hs = hs
        self.logger.info("sunset angle, hs, set to %f", hs)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 1.10.3, Duffy & Beckman (1993)
        ra_d = (86400.0/numpy.pi)*kGsc*dr*(ru*pir*hs + rv*self.dsin(hs))
        self.ra_d = ra_d
        self.logger.info("daily ET radiation set to %f MJ/m^2", (1.0e-6)*ra_d)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate transmittivity (tau), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
        tau_o = (kc + kd*sf)
        tau = tau_o*(1.0 + (2.67e-5)*elv)
        self.tau = tau
        self.logger.info("base transmittivity set to %f", tau_o)
        self.logger.info("transmittivity set to %f", tau)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra_d
        self.ppfd_d = ppfd_d
        self.logger.info("daily PPFD set to %f mol/m^2", ppfd_d)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation (rnl), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        rnl = (kb + (1.0 - kb)*sf)*(kA - tc)
        self.rnl = rnl
        self.logger.info("net longwave radiation set to %f W/m^2", rnl)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate variable substitute (rw), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = (1.0 - kalb_sw)*tau*kGsc*dr
        self.rw = rw
        self.logger.info("variable substitute, rw, set to %f", rw)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate net radiation cross-over hour angle (hn), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rnl - rw*ru)/(rw*rv) >= 1.0:
            # Net radiation negative all day
            self.logger.debug("net radiation negative all day")
            hn = 0
        elif (rnl - rw*ru)/(rw*rv) <= -1.0:
            # Net radiation positive all day
            self.logger.debug("net radiation positive all day")
            hn = 180.0
        else:
            hn = (rnl - rw*ru)/(rw*rv)
            hn = numpy.arccos(hn)
            hn /= pir
        self.hn = hn
        self.logger.info("net radiation cross-over angle set to %f", hn)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate daytime net radiation (rn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rn_d = (86400.0/numpy.pi)*(hn*pir*(rw*ru - rnl) + rw*rv*self.dsin(hn))
        self.rn_d = rn_d
        self.logger.info("daytime net radiation set to %f MJ/m^2", (1.0e-6)*rn_d)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate nighttime net radiation (rnn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rnn_d = rw*ru*(hs - hn)*pir
        rnn_d += rw*rv*(self.dsin(hs) - self.dsin(hn))
        rnn_d += rnl*(numpy.pi - 2.0*hs*pir + hn*pir)
        rnn_d *= (86400.0/numpy.pi)
        self.rnn_d = rnn_d
        self.logger.info("nighttime net radiation set to %f MJ/m^2", (1.0e-6)*rnn_d)

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
        self.logger.debug("calculating cosine of %f degrees", x)
        return numpy.cos(x*pir)

    def dsin(self, x):
        """
        Name:     EVAP.dsin
        Input:    float, angle, degrees (x)
        Output:   float, sin(x*pi/180)
        Features: Calculates the sine of an angle given in degrees
        """
        self.logger.debug("calculating sine of %f degrees", x)
        return numpy.sin(x*pir)

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
        self.logger.debug("calculating heliocentric longitudes for day %d", n)

        # Variable substitutes:
        xee = ke**2
        xec = ke**3
        xse = numpy.sqrt(1.0 - xee)

        # Mean longitude for vernal equinox:
        xlam = (ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(komega)
        xlam -= xee/4.0*(0.5 + xse)*self.dsin(2.0*komega)
        xlam += xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*komega)
        xlam *= 2.0
        xlam /= pir
        self.logger.debug("mean longitude for vernal equinox set to %f", xlam)

        # Mean longitude for day of year:
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        self.logger.debug("mean longitude for day of year set to %f", dlamm)

        # Mean anomaly:
        anm = (dlamm - komega)
        ranm = (anm*pir)
        self.logger.debug("mean anomaly set to %f", ranm)

        # True anomaly:
        ranv = ranm
        ranv += (2.0*ke - xec/4.0)*numpy.sin(ranm)
        ranv += 5.0/4.0*xee*numpy.sin(2.0*ranm)
        ranv += 13.0/12.0*xec*numpy.sin(3.0*ranm)
        anv = ranv/pir

        # True longitude:
        my_tls = anv + komega
        if my_tls < 0:
            my_tls += 360.0
        elif my_tls > 360:
            my_tls -= 360.0
        self.logger.debug("true longitude set to %f", my_tls)

        # True anomaly:
        my_nu = (my_tls - komega)
        if my_nu < 0:
            my_nu += 360.0
        self.logger.debug("true anomaly set to %f", my_nu)

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
        self.logger.debug("calculating Julian day")
        if m <= 2.0:
            y -= 1.0
            m += 12.0

        a = int(y/100)
        b = 2 - a + int(a/4)

        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
        return jde

###############################################################################
## MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    fh = logging.handlers.RotatingFileHandler("solar.log", backupCount=9)
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    fh.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(fh)

    # Test one-year of SPLASH:
    my_lat = 37.7
    my_elv = 142.
    my_day = 172
    my_year = 2000
    my_sf = 1.0
    my_temp = 23.0
    # def __init__(self, lat, n, elv=0.0, y=0, sf=1.0, tc=23.0):
    my_class = SOLAR(my_lat, my_day, my_elv, my_year, my_sf, my_temp)
    fh.doRollover()