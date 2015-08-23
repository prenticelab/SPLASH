#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# const.py
#
# 2014-01-30 -- created
# 2015-08-22 -- last updated
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J. 
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led 
# algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
# transpiration and plant-available moisture, Geoscientific Model Development, 
# 2015 (in progress)

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
kR = 8.31447   # universal gas constant, J/mol/K (Moldover et al., 1988)
kTo = 288.15   # base temperature, K (Berberan-Santos et al., 1997)
kWm = 150      # soil moisture capacity, mm (Cramer & Prentice, 1988)
kw = 0.26      # entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
komega = 283.0 # longitude of perihelion for 2000 CE, degrees (Berger, 1978)