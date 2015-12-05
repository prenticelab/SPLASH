# R version 3.2.2 (2015-08-14) -- "Fire Safety"
#
# const.R
#
# last updated: 2015-12-04
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Robust indices of radiation
# evapo-transpiration and plant-available moisture, Geoscientific Model
# Development, 2015 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains the global constants defined in SPLASH.
#
# NOTE: orbital parameters: eccentricity, obliquity, and longitude of the
# perihelion, are assumed constant while they infact vary slightly over time.
# There are methods for their calculation (e.g., Meeus, 1991). Eccentricity
# varies 0.005--0.072 and is decreasing at rate of 0.00004 per century.
# Obliquity varies 22.1--24.5 degrees with a period of ~41000 years.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - updated values and references for ka and kR [14.10.31]
# - reduced list of constants [15.01.13]
# - updated kR and kTo values and references [15.03.24]
#
kA <- 107           # constant for Rl (Monteith & Unsworth, 1990)
kalb_sw <- 0.17     # shortwave albedo (Federer, 1968)
kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
kb <- 0.20          # constant for Rl (Linacre, 1968; Kramer, 1957)
kc <- 0.25          # constant for Rs (Linacre, 1968)
kCw <- 1.05         # supply constant, mm/hr (Federer, 1982)
kd <- 0.50          # constant for Rs (Linacre, 1968)
ke <- 0.01670       # eccentricity of earth's orbit, 2000CE (Berger 1978)
keps <- 23.44       # obliquity of earth's elliptic, 2000CE (Berger 1978)
kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
kGsc <- 1360.8      # solar constant, W/m^2 (Kopp & Lean, 2011)
kL <- 0.0065        # adiabatic lapse rate, K/m (Cavcar, 2000)
kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kMv <- 0.01802      # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
komega <- 283       # lon. of perihelion, degrees, 2000CE (Berger, 1978)
kSecInDay <- 86400  # number of seconds in a day
kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)
kWm <- 150          # soil moisture capacity, mm (Cramer-Prentice, 1988)
kw <- 0.26          # PET entrainment, (1+kw)*EET (Priestley-Taylor, 1972)
pir <- pi/180       # pi in radians
