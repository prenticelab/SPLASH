#!/usr/bin/python
#
# evap.py
#
# VERSION: 1.1-dev
# LAST UPDATED: 2016-02-19
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
#
# This file is part of the SPLASH model.
#
# SPLASH is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# SPLASH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
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
#from scipy.stats import nanmean

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
        self.logger.debug("EVAP class called")

        # Assign default public variables:
        self.elv = elv
        self.logger.debug("elevation set to %f m", elv)
        # Calculate the size of the output variables:
        if numpy.ndim(elv) == 2:
            lond = 1
        elif numpy.ndim(elv) == 1:
            lond = 2
        else:
            lond = 3
        
        no_lons = elv.shape[lond]
        print no_lons
        print lond
        
        # Create SOLAR class:
        try:
            self.solar = SOLAR(lat, elv)
        except:
            self.logger.exception("failed to initialize SOLAR class")
        else:
            self.logger.debug("initialized solar class")

         # Initialize class variables:
        self.sat = numpy.zeros([no_lons,])     # slope of saturation vap press temp curve, Pa/K
        self.lv = numpy.zeros([no_lons,])    # enthalpy of vaporization, J/kg
        self.pw = numpy.zeros([no_lons,])    # density of water, kg/m^3
        self.psy = numpy.zeros([no_lons,])   # psychrometric constant, Pa/K
        self.econ = numpy.zeros([no_lons,])  # water-to-energy conversion factor
        self.cond = numpy.zeros([no_lons,])  # daily condensation, mm
        self.eet_d = numpy.zeros([no_lons,]) # daily EET, mm
        self.pet_d = numpy.zeros([no_lons,]) # daily PET, mm
        self.rx = numpy.zeros([no_lons,])    # variable substitute, (mm/hr)/(W/m^2)
        self.hi = numpy.zeros([no_lons,])    # intersection hour angle (hi), degrees
        self.aet_d = numpy.zeros([no_lons,]) # daily AET (aet_d), mm

        # Initialize class variables:
        #self.sat = None     # slope of saturation vap press temp curve, Pa/K
        #self.lv = None    # enthalpy of vaporization, J/kg
        #self.pw = None    # density of water, kg/m^3
        #self.psy = None   # psychrometric constant, Pa/K
        #self.econ = None  # water-to-energy conversion factor
        #self.cond = None  # daily condensation, mm
        #self.eet_d = None # daily EET, mm
        #self.pet_d = None # daily PET, mm
        #self.rx = None    # variable substitute, (mm/hr)/(W/m^2)
        #self.hi = None    # intersection hour angle (hi), degrees
        #self.aet_d = None # daily AET (aet_d), mm

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
        sw_neg = numpy.where(sw < 0)
        
        if numpy.sum(sw_neg[0]) > 0:
            self.logger.error("supply rate is outside range of validity")
            raise ValueError("Please provide a valid evporative supply rate")
        else:
            self.logger.debug("Supply rate valid")

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate radiation
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        try:
            self.solar.calculate_daily_fluxes(n, y, sf, tc)
        except:
            self.logger.exception("failed to calculate daily radiation fluxes")
            raise
        else:
            ru = self.solar.ru.copy()
            rv = self.solar.rv.copy()
            rw = self.solar.rw.copy()
            rnl = self.solar.rnl.copy()
            hn = self.solar.hn.copy()
            self.netrad = self.solar.hn.copy()
            rn_d = self.solar.rn_d.copy()
            rnn_d = self.solar.rnn_d.copy()
            self.ppfd_d = self.solar.ppfd_d.copy()
            self.logger.debug(
                ("calculating daily evaporative fluxes for day %d of %d for "
                 "year %d ") % (
                    n, self.solar.kN, self.solar.year))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate water-to-energy conversion (econ), m^3/J
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Slope of saturation vap press temp curve, Pa/K
        s = self.sat_slope(tc)
        self.sat = s
        self.logger.debug("slope of saturation, s, set to %f Pa/K", numpy.nanmean(s))

        # Enthalpy of vaporization, J/kg
        lv = self.enthalpy_vap(tc)
        self.lv = lv
        self.logger.debug("enthalpy of vaporization set to %f MJ/kg", (1e-6)*numpy.nanmean(lv))

        # Density of water, kg/m^3
        pw = self.density_h2o(tc, self.elv2pres(self.elv))
        self.pw = pw
        self.logger.debug("density of water set to %f kg/m^3", numpy.nanmean(pw))

        # Psychrometric constant, Pa/K
        g = self.psychro(tc, self.elv2pres(self.elv))
        self.psy = g
        self.logger.debug("psychrometric constant set to %f Pa/K", numpy.nanmean(g))

        econ = s/(lv*pw*(s + g))
        self.econ = econ
        self.logger.debug("Econ set to %f mm^3/J", (1e9)*numpy.nanmean(econ))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate daily condensation (cn), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cn = (1e3)*econ*numpy.abs(rnn_d)
        self.cond = cn
        self.logger.debug("daily condensation set to %f mm", numpy.nanmean(cn))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Estimate daily EET (eet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eet_d = (1e3)*econ*rn_d
        self.eet_d = eet_d
        self.logger.debug("daily EET set to %f mm", numpy.nanmean(eet_d))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Estimate daily PET (pet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pet_d = (1.0 + kw)*eet_d
        self.pet_d = pet_d
        self.logger.debug("daily PET set to %f mm", numpy.nanmean(pet_d))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = (3.6e6)*(1.0 + kw)*econ
        self.rx = rx
        self.logger.debug("variable substitute, rx, set to %f", numpy.nanmean(rx))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the intersection hour angle (hi), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cos_hi = sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
        
        self.supply_ex_dem = numpy.where(cos_hi >= 1.0)
        self.supply_lim_dem = numpy.where(cos_hi <= -1.0)
        #self.other_sup_dem = numpy.where(((cos_hi  < numpy.float64(1.0)) & (cos_hi > numpy.float64(-1.0)))
        #                                 | numpy.isnan(cos_hi))
        
        
        hi = (numpy.arccos(cos_hi))/pir
        
        if type(hi) == numpy.float64:
             if cos_hi.all >= 1.0:
                 # Supply exceeds demand:
                 hi = 0.0
             elif cos_hi.all <= -1.0:
                 # Supply limits demand everywhere:
                 hi = 180.0
             else:
                 hi = numpy.arccos(cos_hi)
                 hi /= pir
             
        else:
            hi[self.supply_ex_dem] = 0.0
            hi[self.supply_lim_dem] = 180.0
        

        self.hi = hi

        #
        # self.hi = hi
        self.logger.debug("intersection hour angle, hi, set to %f", numpy.nanmean(hi))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Estimate daily AET (aet_d), mm
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        aet_d = sw*hi*pir
        aet_d += rx*rw*rv*(dsin(hn) - dsin(hi))
        aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*pir
        aet_d *= (24.0/numpy.pi)
        self.aet_d = aet_d

        self.logger.debug("daily AET set to %f mm", numpy.nanmean(aet_d))

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
            "calculating temperature dependency at %f degrees", numpy.nanmean(tc))
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
            "calculating temperature dependency at %f degrees", numpy.nanmean(tc))
        enth_out = (1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2)
        return enth_out

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
        self.logger.debug("estimating atmospheric pressure at %f m", numpy.nanmean(z))
        p = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
        #p = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(8.3143*kL))
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
        # self.logger.debug(
        #     ("calculating density of water at temperature %f Celcius and "
        #      "pressure %f Pa") % (tc, p))

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
        # self.logger.debug("water density at 1 atm calculated as %f kg/m^3", po)

        # Calculate bulk modulus at 1 atm (bar):
        ko = 19652.17
        ko += 148.1830*tc
        ko += -2.29995*tc*tc
        ko += 0.01281*tc*tc*tc
        ko += -(4.91564e-5)*tc*tc*tc*tc
        ko += (1.035530e-7)*tc*tc*tc*tc*tc
        # self.logger.debug("bulk modulus at 1 atm calculated as %f bar", ko)

        # Calculate temperature dependent coefficients:
        ca = 3.26138
        ca += (5.223e-4)*tc
        ca += (1.324e-4)*tc*tc
        ca += -(7.655e-7)*tc*tc*tc
        ca += (8.584e-10)*tc*tc*tc*tc
        # self.logger.debug("temperature coef, Ca, calculated as %f", ca)

        cb = (7.2061e-5)
        cb += -(5.8948e-6)*tc
        cb += (8.69900e-8)*tc*tc
        cb += -(1.0100e-9)*tc*tc*tc
        cb += (4.3220e-12)*tc*tc*tc*tc
        # self.logger.debug("temperature coef, Cb, calculated as %f bar^-1", cb)

        # Convert atmospheric pressure to bar (1 bar = 100000 Pa)
        pbar = (1.0e-5)*p
        # self.logger.debug("atmospheric pressure calculated as %f bar", pbar)

        pw = (ko + ca*pbar + cb*pbar**2.0)
        pw /= (ko + ca*pbar + cb*pbar**2.0 - pbar)
        pw *= (1e3)*po
        return pw

    def viscosity_h2o(self, tc, p):
        """
        Name:     LUE.viscosity_h2o
        Input:    - float, ambient temperature (tc), degrees C
                  - float, ambient pressure (p), Pa
        Return:   float, viscosity of water (mu), Pa s
        Features: Calculates viscosity of water at a given temperature and
                  pressure.
        Depends:  density_h2o
        Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
                  Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
                  international formulation for the viscosity of H2O, J. Phys.
                  Chem. Ref. Data, Vol. 38(2), pp. 101-125.
        """
        # Define reference temperature, density, and pressure values:
        tk_ast = 647.096      # Kelvin
        rho_ast = 322.0       # kg/m^3
        mu_ast = (1e-6)       # Pa s

        # Get the density of water, kg/m^3
        rho = self.density_h2o(tc, p)

        # Calculate dimensionless parameters:
        tbar = (tc + 273.15)/tk_ast
        tbarx = tbar**(0.5)
        tbar2 = tbar**2
        tbar3 = tbar**3
        rbar = rho/rho_ast

        # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
        mu0 = 1.67752
        mu0 += 2.20462/tbar
        mu0 += 0.6366564/tbar2
        mu0 += -0.241605/tbar3
        mu0 = 1e2*tbarx/mu0

        # Create Table 3, Huber et al. (2009):
        hj0 = (0.520094, 0.0850895, -1.08374, -0.289555, 0., 0.)
        hj1 = (0.222531, 0.999115, 1.88797, 1.26613, 0., 0.120573)
        hj2 = (-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.)
        hj3 = (0.161913,  0.257399, 0., 0., 0., 0.)
        hj4 = (-0.0325372, 0., 0., 0.0698452, 0., 0.)
        hj5 = (0., 0., 0., 0., 0.00872102, 0.)
        hj6 = (0., 0., 0., -0.00435673, 0., -0.000593264)
        h = hj0 + hj1 + hj2 + hj3 + hj4 + hj5 + hj6
        h_array = numpy.reshape(numpy.array(h), (7, 6))

        # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
        mu1 = 0.
        ctbar = (1./tbar) - 1.
        for i in range(6):
            coef1 = numpy.power(ctbar, i)
            coef2 = 0.
            for j in range(7):
                coef2 += h_array[j][i]*numpy.power((rbar - 1.), j)
            mu1 += coef1*coef2
        mu1 = numpy.exp(rbar*mu1)

        # Calculate mu_bar (Eq. 2, Huber et al., 2009)
        #   assumes mu2 = 1
        mu_bar = mu0*mu1

        # Calculate mu (Eq. 1, Huber et al., 2009)
        mu = mu_bar*mu_ast    # Pa s

        return mu


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
        # self.logger.debug(
        #     ("calculating psychrometric constant at temperature %f Celcius "
        #      "and pressure %f Pa") % (tc, p))

        # Calculate the specific heat capacity of water, J/kg/K
        # Eq. 47, Tsilingiris (2008)

       
        #t_pos = numpy.where(tc >= 0.0)

        if type(tc) == numpy.float64:
            if tc < 0.0:
                cp = 1.013 * 1e3 #J/kg/K
            else:
                cp = 1.0045714270
                cp += (2.050632750e-3)*tc
                cp += -(1.631537093e-4)*tc*tc
                cp += (6.212300300e-6)*tc*tc*tc
                cp += -(8.830478888e-8)*tc*tc*tc*tc
                cp += (5.071307038e-10)*tc*tc*tc*tc*tc
                cp *= (1e3)
        
        else:
             t_neg = numpy.where(tc < 0.0)
             cp = 1.0045714270
             cp += (2.050632750e-3)*tc
             cp += -(1.631537093e-4)*tc*tc
             cp += (6.212300300e-6)*tc*tc*tc
             cp += -(8.830478888e-8)*tc*tc*tc*tc
             cp += (5.071307038e-10)*tc*tc*tc*tc*tc
             cp *= (1e3)
        
        #Assuming the temeprature dependance of Cp is negligable when t < 0 C (MS revision: Appendix B)
             cp[t_neg] = 1.013 * 1e3 #J/kg/K
        # self.logger.debug("specific heat capacity calculated as %f J/kg/K", cp)

        # Calculate latent heat of vaporization, J/kg
        lv_4_g = self.enthalpy_vap(tc)
        # self.logger.debug(
        #     "enthalpy of vaporization calculated as %f MJ/kg", (1e-6)*lv)

        # Calculate psychrometric constant, Pa/K
        # Eq. 8, Allen et al. (1998)
        g = (cp*kMa*p/(kMv*lv_4_g))
        phyc = g.copy()

        neg_g = numpy.where(g < 0.0)
        
        #phyc[neg_g] *= 0.0
        
        return phyc

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.debug)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("evap.log")
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
