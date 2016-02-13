#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# solar_grid.py
#
# LAST UPDATED: 2016-02-06
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
import datetime
import logging

import numpy

from const import (ke, keps, kGsc, kA, kb, kc, kd, kfFEC, kalb_vis, kalb_sw,
                   komega, pir)
from data_grid import DATA_G
from utilities import calculate_latitude
from utilities import dcos
from utilities import dsin
from utilities import get_x_y


###############################################################################
# CLASSES
###############################################################################
class SOLAR_G:
    """
    Name:     SOLAR_G
    Features: This class calculates the daily radiation fluxes.
    History:  Version 1.0.0-dev
              - added print vals functions [16.01.29]
              - moved dcos and dsin to utilites [16.02.06]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, elv, error_val=numpy.inf):
        """
        Name:     SOLAR_G.__init__
        Inputs:   - ndarray, elevation, m (elv)
                  - error value (error_val)
        Features: Initializes grid-based solar radiation class
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("SOLAR_G class called")

        # Assign default public variables:
        self.elv = elv
        self.logger.info("saving %d elevation points", elv.size)

        self.good_idx = numpy.where(elv != error_val)
        self.noval_idx = numpy.where(elv == error_val)

        # Calculate latitude array:
        self.logger.info("calculating latitude array")
        my_y = numpy.array([j for j in range(360)])
        lat_array = calculate_latitude(my_y, 0.5)
        self.lat = numpy.reshape(numpy.repeat(lat_array, 720), (360, 720), 'C')

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def calculate_daily_fluxes(self, n, y, sf, tc):
        """
        Name:     SOLAR_G.calculate_daily_fluxes
        Inputs:   - int, day of the year (n)
                  - int, year (y)
                  - nd.array, fraction of sunshine hours (sf)
                  - nd.array, mean daily air temperature, C (tc)
        Features: Calculates the daily radiation fluxes
        Depends:  - julian_day
                  - berger_tls
                  - dcos
                  - dsin
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Validate day of year
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if n < 1 or n > 366:
            self.logger.error(
                "Day of year outside range of validity, (1 to 366)!")
            raise ValueError(
                "Day of year outside range of validity (1 to 366)!")
        else:
            self.day = n

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            self.logger.debug("set default year to 2000")
            kN = 365
            self.year = 2001
        elif y < 0:
            self.logger.error("year set out of range")
            raise ValueError(
                "Please use a valid Julian or Gregorian calendar year")
        else:
            self.logger.debug("calculating days in year %d", y)
            kN = self.julian_day((y+1), 1, 1) - self.julian_day(y, 1, 1)
            self.year = y
        self.kN = kN
        self.logger.info(
            ("calculating daily radiation fluxes for day %d of %d "
             "for year %d") % (n, kN, self.year))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        self.my_nu = my_nu
        self.my_lambda = my_lambda
        self.logger.info("true anomaly, nu, set to %f degrees", my_nu)
        self.logger.info("true lon, lambda, set to %f degrees", my_lambda)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        kee = ke**2
        my_rho = (1.0 - kee)/(1.0 + ke*dcos(my_nu))
        dr = (1.0/my_rho)**2
        self.dr = dr
        self.logger.info("relative Earth-Sun distance, rho, set to %f", my_rho)
        self.logger.info("distance factor, dr, set to %f", dr)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        delta = numpy.arcsin(dsin(my_lambda)*dsin(keps))
        delta /= pir
        self.delta = delta
        self.logger.info("declination, delta, set to %f", delta)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.logger.debug("calculating variable substitutes ru and rv")
        ru = dsin(self.lat)
        ru *= dsin(delta)
        rv = dcos(self.lat)
        rv *= dcos(delta)
        self.ru = ru
        self.rv = rv

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the sunset hour angles (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 3.22, Stine & Geyer (2001)
        self.logger.debug("calculating sunset hour angles")
        hs = numpy.zeros(shape=(360, 720))

        # Indexes of pixels under polar day conditions:
        # Note: hs == 180 degrees
        hs_pos = numpy.where(ru/rv >= 1.0)
        hs[hs_pos] += 180.0

        # Indexes of pixels for regular conditions:
        # Note: hs = acos(-u/v)
        hs_reg = numpy.where((ru/rv < 1.0) & (ru/rv > -1.0))
        hs[hs_reg] -= 1.0
        hs[hs_reg] *= ru[hs_reg]
        hs[hs_reg] /= rv[hs_reg]
        hs[hs_reg] = numpy.arccos(hs[hs_reg])
        hs[hs_reg] /= pir
        self.hs = hs

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 1.10.3, Duffy & Beckman (1993)
        self.logger.debug("calculating daily ET radiation")
        ra_d = (ru*pir)*hs
        ra_d += rv*dsin(hs)
        ra_d *= (86400.0/numpy.pi)*kGsc*dr
        self.ra_d = ra_d

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate transmittivity (tau), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
        # NOTE: noval indexes have transmittivity of 1.0
        tau_o = numpy.ones(shape=(360, 720))
        tau = numpy.ones(shape=(360, 720))

        self.logger.debug("calculating transmittivity")
        tau_o[self.good_idx] = (kc + kd*sf[self.good_idx])
        tau[self.good_idx] = tau_o[self.good_idx]
        tau[self.good_idx] *= (1.0 + (2.67e-5)*self.elv[self.good_idx])
        self.tau = tau

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate daily PPFD (ppfd_d), mol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.logger.debug("calculating daily PPFD")
        ppfd_d = tau*ra_d
        ppfd_d *= (1.0e-6)*kfFEC
        ppfd_d *= (1.0 - kalb_vis)
        self.ppfd_d = ppfd_d

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Estimate net longwave radiation (rnl), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        # NOTE: noval indexes have 0 longwave radiation
        self.logger.debug("calculating longwave radiation")
        rnl = numpy.zeros(shape=(360, 720))
        rnl[self.good_idx] += kA
        rnl[self.good_idx] -= tc[self.good_idx]
        rnl[self.good_idx] *= (kb + (1.0 - kb)*sf[self.good_idx])
        self.rnl = rnl

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate variable substitute (rw), W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.logger.debug("calculating variable substitute, rw")
        rw = tau
        rw *= kGsc*dr
        rw *= (1.0 - kalb_sw)
        self.rw = rw

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate net radiation cross-over hour angle (hn), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create zeroed cross-over angle:
        self.logger.debug("calculating cross-over hour angle")
        hn = numpy.zeros(shape=(360, 720))

        # Define cosine of hn:
        cos_hn = numpy.copy(rnl)
        cos_hn -= (rw*ru)
        cos_hn /= rw
        cos_hn /= rv

        # Indexes of pixels where Rnl is all-day positive:
        hn_pos = numpy.where(cos_hn <= -1.0)
        hn[hn_pos] += 180.0

        # Indexes of pixels for regular Rnl:
        hn_reg = numpy.where((cos_hn < 1.0) & (cos_hn > -1.0))
        hn[hn_reg] = numpy.arccos(cos_hn[hn_reg])
        hn[hn_reg] /= pir
        self.hn = hn

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate daytime net radiation (rn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # NOTE: noval indexes of daytime net radiation set to 0
        self.logger.debug("calculating daytime net radiation")
        rn_d = numpy.zeros(shape=(360, 720))
        rn_d[self.good_idx] += rw[self.good_idx]
        rn_d[self.good_idx] *= ru[self.good_idx]
        rn_d[self.good_idx] -= rnl[self.good_idx]
        rn_d[self.good_idx] *= hn[self.good_idx]*pir
        rn_d[self.good_idx] += (
            rw[self.good_idx]*rv[self.good_idx]*dsin(hn[self.good_idx]))
        rn_d[self.good_idx] *= (86400.0/numpy.pi)
        self.rn_d = rn_d

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate nighttime net radiation (rnn_d), J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.logger.debug("calculating nighttime net radiation")
        rnn_d = numpy.zeros(shape=(360, 720))
        rnn_d[self.good_idx] += hs[self.good_idx]
        rnn_d[self.good_idx] -= hn[self.good_idx]
        rnn_d[self.good_idx] *= pir
        rnn_d[self.good_idx] *= rw[self.good_idx]
        rnn_d[self.good_idx] *= ru[self.good_idx]
        rnn_d[self.good_idx] += (
            rw[self.good_idx]*rv[self.good_idx]*dsin(hs[self.good_idx]))
        rnn_d[self.good_idx] -= (
            rw[self.good_idx]*rv[self.good_idx]*dsin(hn[self.good_idx]))
        rnn_d[self.good_idx] += rnl[self.good_idx]*numpy.pi
        rnn_d[self.good_idx] -= rnl[self.good_idx]*2.0*pir*hs[self.good_idx]
        rnn_d[self.good_idx] += rnl[self.good_idx]*hn[self.good_idx]*pir
        rnn_d[self.good_idx] *= (86400.0/numpy.pi)
        self.rnn_d = rnn_d

    def berger_tls(self, n):
        """
        Name:     SOLAR_G.berger_tls
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
        xlam = (ke/2.0 + xec/8.0)*(1.0 + xse)*dsin(komega)
        xlam -= xee/4.0*(0.5 + xse)*dsin(2.0*komega)
        xlam += xec/8.0*(1.0/3.0 + xse)*dsin(3.0*komega)
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
        Name:     SOLAR_G.julian_day
        Input:    - int, year (y)
                  - int, month (m)
                  - int, day of month (i)
        Output:   float, Julian Ephemeris Day
        Features: Converts Gregorian date (year, month, day) to Julian
                  Ephemeris Day
        Ref:      Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," Astronomical
                  Algorithms
        """
        self.logger.debug("calculating Julian day for %d-%d-%d" % (y, m, i))
        if m <= 2.0:
            y -= 1.0
            m += 12.0

        a = int(y/100)
        b = 2 - a + int(a/4)

        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
        return jde

    def print_vals(self, lon, lat):
        """
        Name:     SOLAR_G.print_vals
        Inputs:   - float, longitude, degrees (lon)
                  - float, latitude, degrees (lat)
        Outputs:  None.
        Features: Prints daily radiation fluxes
        Depends:  get_x_y
        """
        x, y = get_x_y(lon, lat)

        print("year: %d" % (self.year))
        print("day of year, n: %d" % (self.day))
        print("days in year, kN: %d" % (self.kN))
        print("true anomaly, nu: %0.6f degrees" % (self.my_nu))
        print("true lon, lambda: %0.6f degrees" % (self.my_lambda))
        print("distance factor, dr: %0.6f" % (self.dr))
        print("declination, delta: %0.6f" % (self.delta))
        print("variable substitute, ru: %0.6f" % (self.ru[y, x]))
        print("variable substitute, rv: %0.6f" % (self.rv[y, x]))
        print("sunset angle, hs: %0.6f degrees" % (self.hs[y, x]))
        print("daily ET radiation: %0.6f MJ/m^2" % ((1.0e-6)*self.ra_d[y, x]))
        print("transmittivity, tau: %0.6f" % (self.tau[y, x]))
        print("daily PPFD: %0.6f mol/m^2" % (self.ppfd_d[y, x]))
        print("net longwave radiation: %0.6f W/m^2" % (self.rnl[y, x]))
        print("variable substitute, rw: %0.6f" % (self.rw[y, x]))
        print("cross-over hour angle: %0.6f degrees" % (self.hn[y, x]))
        print("daytime net radiation: %0.6f MJ/m^2" % (
            (1.0e-6)*self.rn_d[y, x]))
        print("nighttime net radiation: %0.6f MJ/m^2" % (
            (1.0e-6)*self.rnn_d[y, x]))


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("solar_grid.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    cru_dir = "/usr/local/share/data/cru"
    root_logger.info("creating data class and loading data files")
    data = DATA_G()
    data.find_cru_files(cru_dir)
    data.read_elv()
    root_logger.info("finished loading data files")

    n = 172
    y = 2000
    my_date = datetime.date(y, 1, 1)
    my_date += datetime.timedelta(days=(n-1))
    root_logger.info("current date set to %s", my_date)

    root_logger.info("reading monthly data")
    data.read_monthly_clim(my_date)
    root_logger.info("finished reading monthly data")

    root_logger.info("creating solar class")
    my_class = SOLAR_G(data.elv, data.error_val)

    root_logger.info("calculating daily fluxes")
    my_lon = -122.25
    my_lat = 37.75
    my_class.calculate_daily_fluxes(172, 2000, data.sf, data.tair)

    print("Data:")
    data.print_vals(my_lon, my_lat)
    print("")
    print("Results:")
    my_class.print_vals(my_lon, my_lat)
