#!/usr/bin/python
#
# evap.py
#
# LAST UPDATED: 2016-02-13
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
# led algorithms for simulating habitats (SPLASH): Robust indices of radiation,
# evapotranspiration and plant-available moisture, Geoscientific Model
# Development, 2015 (in progress)

###############################################################################
# IMPORT MODULES:
###############################################################################
import logging

import numpy

from const import (kw, kG, kL, kR, kPo, kTo, kMa, kMv, pir)
from solar import SOLAR
from utilities import dsin


###############################################################################
# CLASSES
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
    Version:  1.0.0-dev
              - replaced radiation methods with SOLAR class [15.12.29]
              - implemented logging [15.12.29]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lat, elv=0.0):
        """
        Name:     EVAP.__init__
        Input:    - float, latitude, degrees (lat)
                  - float, elevation, m (elv)
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("EVAP class called")

        # Assign default public variables:
        self.elv = elv
        self.logger.info("elevation set to %f m", elv)

        # Create SOLAR class:
        try:
            self.solar = SOLAR(lat, elv)
        except:
            self.logger.exception("failed to initialize SOLAR class")
        else:
            self.logger.debug("initialized solar class")

        # Initialize class variables:
        self.sat = None    # slope of saturation vap press temp curve, Pa/K
        self.lv = None     # enthalpy of vaporization, J/kg
        self.pw = None     # density of water, kg/m^3
        self.psy = None    # psychrometric constant, Pa/K
        self.econ = None   # water-to-energy conversion factor
        self.cond = None   # daily condensation, mm
        self.eet_d = None  # daily EET, mm
        self.pet_d = None  # daily PET, mm
        self.rx = None     # variable substitute, (mm/hr)/(W/m^2)
        self.hi = None     # intersection hour angle (hi), degrees
        self.aet_d = None  # daily AET (aet_d), mm

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def calculate_daily_fluxes(self, sw, n, y=0, sf=1.0, tc=23.0):
        """
        Name:     EVAP.calculate_daily_fluxes
        Input:    - float, evaporative supply rate, mm/hr (sw)
                  - int, day of the year (n)
                  - [optional] int, year (y)
                  - [optional] float, fraction of sunshine hours (sf)
                  - [optional] float, mean daily air temperature, C (tc)
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Validate supply rate
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if sw < 0:
            self.logger.error("supply rate is outside range of validity")
            raise ValueError("Please provide a valid evporative supply rate")

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate radiation
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        try:
            self.solar.calculate_daily_fluxes(n, y, sf, tc)
        except:
            self.logger.exception("failed to calculate daily radiation fluxes")
            raise
        else:
            ru = self.solar.ru
            rv = self.solar.rv
            rw = self.solar.rw
            rnl = self.solar.rnl
            hn = self.solar.hn
            rn_d = self.solar.rn_d
            rnn_d = self.solar.rnn_d
            self.logger.info(
                ("calculating daily evaporative fluxes for day %d of %d for "
                 "year %d with sunshine fraction %f, air temperature %f "
                 "Celcuis, and supply rate %f") % (
                    n, self.solar.kN, self.solar.year, sf, tc, sw))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate water-to-energy conversion (econ), m^3/J
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Slope of saturation vap press temp curve, Pa/K
        s = self.sat_slope(tc)
        self.sat = s
        self.logger.info("slope of saturation, s, set to %f Pa/K", s)

        # Enthalpy of vaporization, J/kg
        lv = self.enthalpy_vap(tc)
        self.lv = lv
        self.logger.info("enthalpy of vaporization set to %f MJ/kg", (1e-6)*lv)

        # Density of water, kg/m^3
        pw = self.density_h2o(tc, self.elv2pres(self.elv))
        self.pw = pw
        self.logger.info("density of water set to %f kg/m^3", pw)

        # Psychrometric constant, Pa/K
        g = self.psychro(tc, self.elv2pres(self.elv))
        self.psy = g
        self.logger.info("psychrometric constant set to %f Pa/K", g)

        econ = s/(lv*pw*(s + g))
        self.econ = econ
        self.logger.info("Econ set to %f mm^3/J", (1e9)*econ)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate daily condensation (cn), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cn = (1e3)*econ*numpy.abs(rnn_d)
        self.cond = cn
        self.logger.info("daily condensation set to %f mm", cn)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Estimate daily EET (eet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eet_d = (1e3)*econ*rn_d
        self.eet_d = eet_d
        self.logger.info("daily EET set to %f mm", eet_d)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Estimate daily PET (pet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pet_d = (1.0 + kw)*eet_d
        self.pet_d = pet_d
        self.logger.info("daily PET set to %f mm", pet_d)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = (3.6e6)*(1.0 + kw)*econ
        self.rx = rx
        self.logger.info("variable substitute, rx, set to %f", rx)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the intersection hour angle (hi), degrees
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
        self.logger.info("intersection hour angle, hi, set to %f", hi)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Estimate daily AET (aet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        aet_d = sw*hi*pir
        aet_d += rx*rw*rv*(dsin(hn) - dsin(hi))
        aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*pir
        aet_d *= (24.0/numpy.pi)
        self.aet_d = aet_d
        self.logger.info("daily AET set to %f mm", aet_d)

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
        self.logger.debug(
            "calculating temperature dependency at %f degrees", tc)
        return s

    def enthalpy_vap(self, tc):
        """
        Name:     EVAP.enthalpy_vap
        Input:    float, air temperature (tc), degrees C
        Output:   float, latent heat of vaporization
        Features: Calculates the enthalpy of vaporization, J/kg
        Ref:      Eq. 8, Henderson-Sellers (1984)
        """
        self.logger.debug(
            "calculating temperature dependency at %f degrees", tc)
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
        self.logger.debug("estimating atmospheric pressure at %f m", z)
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
        self.logger.debug(
            ("calculating density of water at temperature %f Celcius and "
             "pressure %f Pa") % (tc, p))

        # Calculate density at 1 atm (kg/m^3):
        po = 0.99983952
        po += (6.788260e-5)*tc
        po += -(9.08659e-6)*tc*tc
        po += (1.022130e-7)*tc*tc*tc
        po += -(1.35439e-9)*tc*tc*tc*tc
        po += (1.471150e-11)*tc*tc*tc*tc*tc
        po += -(1.11663e-13)*tc*tc*tc*tc*tc*tc
        po += (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc
        po += -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
        self.logger.info("water density at 1 atm calculated as %f kg/m^3", po)

        # Calculate bulk modulus at 1 atm (bar):
        ko = 19652.17
        ko += 148.1830*tc
        ko += -2.29995*tc*tc
        ko += 0.01281*tc*tc*tc
        ko += -(4.91564e-5)*tc*tc*tc*tc
        ko += (1.035530e-7)*tc*tc*tc*tc*tc
        self.logger.info("bulk modulus at 1 atm calculated as %f bar", ko)

        # Calculate temperature dependent coefficients:
        ca = 3.26138
        ca += (5.223e-4)*tc
        ca += (1.324e-4)*tc*tc
        ca += -(7.655e-7)*tc*tc*tc
        ca += (8.584e-10)*tc*tc*tc*tc
        self.logger.info("temperature coef, Ca, calculated as %f", ca)

        cb = (7.2061e-5)
        cb += -(5.8948e-6)*tc
        cb += (8.69900e-8)*tc*tc
        cb += -(1.0100e-9)*tc*tc*tc
        cb += (4.3220e-12)*tc*tc*tc*tc
        self.logger.info("temperature coef, Cb, calculated as %f bar^-1", cb)

        # Convert atmospheric pressure to bar (1 bar = 100000 Pa)
        pbar = (1.0e-5)*p
        self.logger.info("atmospheric pressure calculated as %f bar", pbar)

        pw = (ko + ca*pbar + cb*pbar**2.0)
        pw /= (ko + ca*pbar + cb*pbar**2.0 - pbar)
        pw *= (1e3)*po
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
        self.logger.debug(
            ("calculating psychrometric constant at temperature %f Celcius "
             "and pressure %f Pa") % (tc, p))

        # Calculate the specific heat capacity of water, J/kg/K
        # Eq. 47, Tsilingiris (2008)
        cp = 1.0045714270
        cp += (2.050632750e-3)*tc
        cp += -(1.631537093e-4)*tc*tc
        cp += (6.212300300e-6)*tc*tc*tc
        cp += -(8.830478888e-8)*tc*tc*tc*tc
        cp += (5.071307038e-10)*tc*tc*tc*tc*tc
        cp *= (1e3)
        self.logger.info("specific heat capacity calculated as %f J/kg/K", cp)

        # Calculate latent heat of vaporization, J/kg
        lv = self.enthalpy_vap(tc)
        self.logger.info(
            "enthalpy of vaporization calculated as %f MJ/kg", (1e-6)*lv)

        # Calculate psychrometric constant, Pa/K
        # Eq. 8, Allen et al. (1998)
        return (cp*kMa*p/(kMv*lv))

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    root_handler = logging.StreamHandler()
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # Test one-year of SPLASH:
    my_lat = 37.7
    my_elv = 142.
    my_day = 172
    my_year = 2000
    my_sf = 1.0
    my_temp = 23.0
    my_sw = 0.9

    my_class = EVAP(my_lat, my_elv)
    my_class.calculate_daily_fluxes(my_sw, my_day, my_year, my_sf, my_temp)
