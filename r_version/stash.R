# R-Studio 0.98
#
# stash.R
#
# written by Tyler W. Davis
# Imperial College London
#
# last updated: 2015-01-30
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script runs the STASH 2.0 code for point-based data.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 01. updated values and references for ka and kR [14.10.31]
# 02. fixed Cooper's and Spencer's declination angle equations [14.11.25]
# 03. replaced simplified_kepler with full_kepler method [14.11.25]
# 04. added berger_tls function [15.01.13]
# 05. reduced list of constants [15.01.13]
# 06. updated evap function (similar to stash.py EVAP class) [15.01.13]
# 07. updated monthly and daily results & process [15.01.13]
# 08. updated plots of results [15.01.16]
# 09. added example data CSV file [15.01.16]
# 10. fixed Cramer-Prentice alpha definition [15.01.16]
# 11. updated monthly results plot [15.01.22]
# 12. added fix to daily soil moisture when n and ny are 365 [15.01.27]
# 13. added write out for daily/monthly results [15.01.27]
# 14. added example of yearly looping [15.01.27]
# 15. updated year extraction from filename in yearly loop example [15.01.30]
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: evap
# *
# * Input: double, latitude, degrees (lat)
# *        double, day of year (n)
# *        double, elevation (elv)  *optional
# *        double, year (y)         *optional
# *        double, fraction of sunshine hours (sf)        *optional
# *        double, mean daily air temperature, deg C (tc) *optional
# *        double, evaporative supply rate, mm/hr (sw)    *optional
# *
# * Return: list object (et.srad)
# *         $nu_deg ............ true anomaly, degrees
# *         $lambda_deg ........ true longitude, degrees
# *         $dr ................ distance factor, unitless
# *         $delta_deg ......... declination angle, degrees
# *         $hs_deg ............ sunset angle, degrees
# *         $ra_j.m2 ........... daily extraterrestrial radiation, J/m^2
# *         $tau ............... atmospheric transmittivity, unitless
# *         $ppfd_mol.m2 ....... daily PPFD, mol/m^2
# *         $hn_deg ............ net radiation hour angle, degrees
# *         $rn_j.m2 ........... daily net radiation, J/m^2
# *         $rnn_j.m2 .......... daily nighttime net radiation, J/m^2
# *         $econ_m3.j ......... water to energy conversion, m^3/J
# *         $cond_mm ........... daily condensation, mm
# *         $eet_mm ............ daily equilibrium ET, mm
# *         $pet_mm ............ daily potential ET, mm
# *         $hi_deg ............ intersection hour angle, degrees
# *         $aet_mm ............ daily actual ET, mm
# *
# * Features: This function calculates daily evaporation rates.
# *
# *           Depends on:
# *           - kalb_sw ........ shortwave albedo
# *           - kalb_vis ....... visible light albedo
# *           - kb ............. empirical constant for longwave rad
# *           - kc ............. empirical constant for shortwave rad
# *           - kd ............. empirical constant for shortwave rad
# *           - ke ............. eccentricity
# *           - keps ........... obliquity
# *           - kfFEC .......... from-flux-to-energy conversion, umol/J
# *           - kGsc ........... solar constant
# *           - kw ............. entrainment factor for PET
# *           - berger_tls() ... calc true anomaly and longitude
# *           - dcos() ......... cos(x*pi/180), where x is in degrees
# *           - dsin() ......... sin(x*pi/180), where x is in degrees
# *           - density_h2o() .. density of water
# *           - elv2pres() ..... elevation dependent atm. pressure
# *           - enthalpy_vap() . latent heat of vaporization
# *           - julian_day() ... date to julian day
# *           - psychro() ...... psychrometric constant
# *           - sat_slope() .... slope of sat. pressure temp curve
# *
# * Ref's:
# * Allen, R.G. (1996), Assessing integrity of weather data for reference
# *   evapotranspiration estimation, Journal of Irrigation and Drainage
# *   Engineering, vol. 122, pp. 97--106.
# * Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998), 'Meteorological 
# *   data,' Crop evapotranspiration - Guidelines for computing crop water 
# *   requirements - FAO Irrigation and drainage paper 56, Food and 
# *   Agriculture Organization of the United Nations, Available:
# *   http://www.fao.org/docrep/x0490e/x0490e07.htm
# * Berger, A.L. (1978), Long-term variations of daily insolation and 
# *   quarternary climatic changes, Journal of Atmospheric Sciences, 
# *   vol. 35, pp. 2362--2367.
# * Berger, A.L., M.F. Loutre, and C. Tricot (1993), Insolation and 
# *   Earth's orbital periods, J. Geophys. Res., 98, 10341--10362.
# * Duffie, J. A. and W. A. Beckman (1991). Solar engineering of thermal 
# *   processes. 4th ed. New Jersey: John Wiley and Sons
# * Federer (1982), Transpirational supply and demand: plant, soil, and 
# *   atmospheric effects evaluated by simulation, Water Resources 
# *   Research, vol. 18, no. 2, pp. 355--362.
# * Ge, S., R.G. Smith, C.P. Jacovides, M.G. Kramer, R.I. Carruthers 
# *   (2011), Dynamics of photosynthetic photon flux density (PPFD) and 
# *   estimates in coastal northern California, Theoretical and Applied 
# *   Climatology, vol. 105, pp. 107--118.
# * Henderson-Sellers, B. (1984), A new formula for latent heat of
# *   vaporization of water as a function of temperature, Quarterly
# *   Journal of the Royal Meteorological Society 110, pp. 1186–1190
# * Linacre (1968), Estimating the net-radiation flux, Agricultural 
# *   Meteorology, vol. 5, pp. 49--63.
# * Prentice, I.C., M.T. Sykes, W. Cramer (1993), A simulation model for 
# *   the transient effects of climate change on forest landscapes, 
# *   Ecological Modelling, vol. 65, pp. 51--70.
# * Priestley, C.H.B. and R.J. Taylor (1972), On the assessment of surface 
# *   heat flux and evaporation using large-scale parameters, Monthly  
# *   Weather Review, vol. 100 (2), pp. 81--92.
# * Spencer, J. W. (1971). “Fourier series representation of the position 
# *   of the sun”. In: Search 2, p. 172.
# * Stine, W. B. and M. Geyer (2001). “Power from the Sun”. In: Available 
# *   online: http://www.powerfromthesun.net/Book/chapter03/chapter03.
# * Wetherald, R.T., S. Manabe (1972), Response to joint ocean-atmosphere 
# *   model to the seasonal variation of the solar radiation, Monthly 
# *   Weather Review, vol. 100 (1), pp. 42--59.
# * Woolf, H. M. (1968). On the computation of solar evaluation angles 
# *   and the determination of sunrise and sunset times. Tech. rep. 
# *   NASA-TM-X-164. National Aeronautics and Space Administration (NASA).
# *
# ************************************************************************
evap <- function(lat, n, elv=0, y=0, sf=1, tc=23.0, sw=1.0){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (lat > 90 || lat < -90){
    stop("Warning: Latitude outside range of validity (-90 to 90)!")
  }
  if (n < 1 || n > 366){
    stop("Warning: Day outside range of validity (1 to 366)!")
  }
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  et.srad <- list()
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the number of days in yeark (kN), days
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (y == 0){
    kN = 365
  } else {
    kN = (julian_day(y+1, 1, 1) - julian_day(y, 1, 1))
  }
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate heliocentric longitudes (nu and lambda), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my_helio <- berger_tls(n, kN)
  my_nu <- my_helio[1]
  my_lambda <- my_helio[2]
  #
  et.srad$nu_deg <- my_nu
  et.srad$lambda_deg <- my_lambda
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate distance factor (dr), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my_rho <- (1 - ke^2)/(1 + ke*dcos(my_nu))
  dr <- (1/my_rho)^2
  #
  et.srad$dr <- dr
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate the declination angle (delta), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  delta <- asin(dsin(my_lambda)*dsin(keps))
  delta <- delta*(180/pi)
  #
  et.srad$delta_deg <- delta
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate variable substitutes (u and v), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ru <- dsin(delta)*dsin(lat)
  rv <- dcos(delta)*dcos(lat)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Calculate the sunset hour angle (hs), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note: u/v == tan(delta)*tan(lat) 
  if (ru/rv >= 1.0){
    hs <- 180  # Polar day (no sunset)
  } else if (ru/rv <= -1.0){ 
    hs <- 0 # Polar night (no sunrise)
  } else {
    hs <- acos(-1.0*ru/rv)
    hs <- hs*(180.0/pi)
  }
  et.srad$hs_deg <- hs
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 07. Calculate daily extraterrestrial radiation (ra_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 1.10.3, Duffy & Beckman (1993)
  ra_d <- (86400.0/pi)*kGsc*dr*(ru*(pi/180)*hs + rv*dsin(hs))
  et.srad$ra_j.m2 <- ra_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 08. Calculate transmittivity (tau), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref:  Eq. 11, Linacre (1968): tau_o
  #       Eq. 2, Allen (1996): tau
  tau_o <- (kc + kd*sf)
  tau <- tau_o*(1 + (2.67e-5)*elv)
  #
  et.srad$tau <- tau
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 09. Calculate daily PPFD (ppfd_d), mol/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ppfd_d <- (1.0e-6)*kfFEC*(1 - kalb_vis)*tau*ra_d
  et.srad$ppfd_mol.m2 <- ppfd_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 10. Estimate net longwave radiation (rnl), W/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rnl <- (kb + (1.0 - kb)*sf)*(kA - tc)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 11. Calculate variable substitue (rw), W/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rw <- (1 - kalb_sw)*tau*kGsc*dr
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 12. Calculate net radiation cross-over angle (hn), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ((rnl - rw*ru)/(rw*rv) >= 1.0){
    hn <- 0  # Net radiation is negative all day
  } else if ((rnl - rw*ru)/(rw*rv) <= -1.0){ 
    hn <- 180 # Net radiation is positive all day
  } else {
    hn <- acos((rnl - rw*ru)/(rw*rv))
    hn <- hn*(180/pi)
  }
  et.srad$hn_deg <- hn
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 13. Calculate daytime net radiation (rn_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rn_d <- (86400.0/pi)*(hn*(pi/180)*(rw*ru - rnl) + rw*rv*dsin(hn))
  et.srad$rn_j.m2 <- rn_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 14. Calculate nighttime net radiation (rnn_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rnn_d <- (86400.0/pi)*(
    rw*ru*(hs - hn)*(pi/180) + rw*rv*(dsin(hs) - dsin(hn)) + 
      rnl*(pi - 2*hs*(pi/180) + hn*(pi/180))
  )
  et.srad$rnn_j.m2 <- rnn_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 15. Calculate water-to-energy conversion (econ), m^3/J
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Slope of saturation vap press temp curve, Pa/K
  s <- sat_slope(tc)
  # Enthalpy of vaporization, J/kg
  lv <- enthalpy_vap(tc)
  # Density of water, kg/m^3
  pw <- density_h2o(tc, elv2pres(elv))
  # Psychrometric constant, Pa/K
  g <- psychro(tc, elv2pres(elv))
  #
  econ <- s/(lv*pw*(s + g))
  et.srad$econ_m3.j <- econ
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 16. Calculate daily condensation (wc), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  wc <- (1e3)*econ*abs(rnn_d)
  et.srad$cond_mm <- wc
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 17. Estimate daily equilibrium evapotranspiration (eet_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eet_d <- (1e3)*econ*rn_d
  et.srad$eet_mm <- eet_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 18. Estimate daily potential evapotranspiration (pet_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pet_d <- (1.0 + kw)*eet_d
  et.srad$pet_mm <- pet_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rx <- (3.6e6)*(1.0 + kw)*econ
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 20. Calculate the intersection hour angle (hi), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cos_hi <- sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
  if (cos_hi >= 1.0){
    hi <- 0.0       # supply exceeds demand
  } else if (cos_hi <= -1.0){
    hi <- 180.0     # supply limits demand everywhere
  } else {
    hi <- acos(cos_hi)
    hi <- hi*(180/pi)
  }
  et.srad$hi_deg <- hi
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 21. Estimate daily actual evapotranspiration (aet_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  aet_d <- (24/pi)*(
    sw*hi*(pi/180) + 
      rx*rw*rv*(dsin(hn) - dsin(hi)) + 
      (rx*rw*ru - rx*rnl)*(hn - hi)*(pi/180)
  )
  #
  et.srad$aet_mm <- aet_d
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  et.srad
}

# ************************************************************************
# * Name: berger_tls
# *
# * Input: double, day of year (n)
# *        double, days in year (N)
# *
# * Return: numeric list, true anomaly and true longitude
# *
# * Features: Returns true anomaly and true longitude for a given day.
# *
# *           Depends on:
# *           - ke ............. eccentricity of earth's orbit, unitless
# *           - komega ......... longitude of perihelion
# *
# * Ref: Berger, A. L. (1978), Long term variations of daily insolation
# *      and quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
# ************************************************************************
berger_tls <- function(n, N){
  # Variable substitutes:
  xee <- ke^2
  xec <- ke^3
  xse <- sqrt(1 - ke^2)
  pir <- pi/180.0
  #
  # Mean longitude for vernal equinox:
  xlam <- (ke/2.0 + xec/8.0)*(1 + xse)*sin(komega*pir) - 
    xee/4.0*(0.5 + xse)*sin(2.0*komega*pir) + 
    xec/8.0*(1.0/3.0 + xse)*sin(3.0*komega*pir)
  xlam <- 2.0*xlam/pir
  #
  # Mean longitude for day of year:
  dlamm <- xlam + (n - 80.0)*(360.0/N)
  #
  # Mean anomaly:
  anm <- dlamm - komega
  ranm <- anm*pir
  #
  # True anomaly (uncorrected):
  ranv <- ranm + (2.0*ke - xec/4.0)*sin(ranm) +
    5.0/4.0*xee*sin(2.0*ranm) + 
    13.0/12.0*xec*sin(3.0*ranm)
  anv <- ranv/pir
  #
  # True longitude:
  my_tls <- anv + komega
  if (my_tls < 0){
    my_tls <- my_tls + 360
  } else if (my_tls > 360) {
    my_tls <- my_tls - 360
  }
  #
  # True anomaly:
  my_nu <- my_tls - komega
  if (my_nu < 0){
    my_nu <- my_nu + 360
  }
  return (c(my_nu, my_tls))
}

# ************************************************************************
# * Name: julian_day
# *
# * Input: double, year (y) 
# *        double, month (m) 
# *        double, day of month (i)
# *
# * Return: double, Julian day
# *
# * Features: This function converts a date in the Gregorian calendar
# *           to a Julian day number (i.e., a method of consecutative 
# *           numbering of days---does not have anything to do with 
# *           the Julian calendar!)
# *
# * Ref: Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", 
# *      Astronomical Algorithms
# ************************************************************************
julian_day <- function(y, m, i){
  if(m <= 2){
    y <- y - 1
    m <- m + 12
  }
  a <- floor(y/100)
  b <- 2 - a + floor(a/4)
  # 
  jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
  return(jde)
}

# ************************************************************************
# * Name: dcos
# *
# * Input: double (d), angle in degrees
# *
# * Return: double, cosine of angle
# *
# * Features: This function calculates the cosine of an angle (d) given
# *           in degrees.
# *
# * Note: This script is based on the Javascript function written by
# *       C Johnson, Theoretical Physicist, Univ of Chicago
# *
# *       'Equation of Time' URL: http://mb-soft.com/public3/equatime.html
# *       Javascript URL: http://mb-soft.com/believe/txx/astro22.js
# ************************************************************************
dcos <- function(d) {
  cos(d*pi/180)
}

# ************************************************************************
# * Name: dsin
# *
# * Input: double (d), angle in degrees
# *
# * Return: double, sine of angle
# *
# * Features: This function calculates the sine of an angle (d) given
# *           in degrees.
# ************************************************************************
dsin <- function(d) {
  sin(d*pi/180)
}

# ************************************************************************
# * Name: elv2pres
# *
# * Input: double (z), meters
# *
# * Return: double, Pa
# *
# * Features: This function calculates the elevation dependent 
# *           atmospheric pressure, Pa
# *
# *           Depends on:
# *           - kPo ............ base pressure, Pa
# *           - kTo ............ base temperature, C
# *           - kL ............. temperature lapse rate, K/m
# *           - kR ............. universal gas constant, J/mol/K
# *           - kMa ............ molecular weight of dry air, kg/mol
# *           - kG ............. gravity, m/s^2
# *
# * Ref: Cavcar (2000), The International Standard Atmosphere (ISA), 
# *        Anadolu University, Turkey.
# ************************************************************************
elv2pres <- function(z){
  kPo*(1 - kL*z/kTo)^(kG*kMa/(kR*kL))
}

# ************************************************************************
# * Name: sat_slope
# *
# * Input: double (tc), degrees C
# *
# * Return: double, Pa/K
# *
# * Features: This function calculates the temperature-dependent slope of
# *           the saturation pressure temperature curve using the 
# *           methodology presented in the eMast energy.cpp script
# *
# * Ref: Eq. 6, Prentice et al. (1993); 
# *      Eq. 13, Allen et al. (1998)
# ************************************************************************
sat_slope <- function(tc){
  (17.269)*(237.3)*(610.78)*exp(17.269*tc/(237.3 + tc))/(237.3 + tc)^2
}

# ************************************************************************
# * Name: enthalpy_vap
# *
# * Input: double (tc), degrees C
# *
# * Return: double, J/kg
# *
# * Features: This function calculates the temperature-dependent enthalpy
# *           of vaporization (latent heat of vaporization)
# *
# * Ref: Eq. 8, Henderson-Sellers (1984), A new formula for latent heat
# *      of vaporization of water as a function of temperature, Quarterly 
# *      Journal of the Royal Meteorological Society, vol. 110, pp. 1186--
# *      1190.
# ************************************************************************
enthalpy_vap <- function(tc){
  1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))^2
}

# ************************************************************************
# * Name: density_h2o
# *
# * Input: double (tc), air temperature, degrees C
# *        double (pa), atm pressure, Pa
# *
# * Return: double, kg/m^3
# *
# * Features: This function calculates the temperature and pressure 
#             dependent density of pure water
# *
# * Ref: Chen, C.T., R.A. Fine, and F.J. Millero (1977), The equation 
#          of state of pure water determined from sound speeds, The 
#          Journal of Chemical Physics 66, 2142; 
#          doi: 10.1063/1.434179
# ************************************************************************
density_h2o <- function(tc, pa){
  # Calculate density of water at 1 atm, g/cm^3
  po <- 0.99983952 + 
    (6.788260e-5)*tc + 
    -(9.08659e-6)*tc*tc +
    (1.022130e-7)*tc*tc*tc + 
    -(1.35439e-9)*tc*tc*tc*tc +
    (1.471150e-11)*tc*tc*tc*tc*tc +
    -(1.11663e-13)*tc*tc*tc*tc*tc*tc + 
    (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc + 
    -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
  #
  # Calculate the bulk modulus of water at 1 atm, atm
  ko <- 19652.17 +
    148.1830*tc + 
    -2.29995*tc*tc + 
    0.01281*tc*tc*tc + 
    -(4.91564e-5)*tc*tc*tc*tc + 
    (1.035530e-7)*tc*tc*tc*tc*tc
  #
  # Calculate temperature-dependend coefficients
  ca <- 3.26138 + 
    (5.223e-4)*tc + 
    (1.324e-4)*tc*tc + 
    -(7.655e-7)*tc*tc*tc + 
    (8.584e-10)*tc*tc*tc*tc
  cb <- (7.2061e-5) +
    -(5.8948e-6)*tc + 
    (8.69900e-8)*tc*tc + 
    -(1.0100e-9)*tc*tc*tc + 
    (4.3220e-12)*tc*tc*tc*tc
  #
  #
  # Convert pressure to bar (1 bar = 100000 Pa)
  pbar <- (1e-5)*pa
  #
  pw <- (1e3)*po*(ko + ca*pbar + cb*pbar^2)/(ko + ca*pbar + cb*pbar^2 - pbar)
  return(pw)
}

# ************************************************************************
# * Name: psychro
# *
# * Input: double (tc), air temperature, degrees C
# *        double (pa), atm pressure, Pa
# *
# * Return: double, Pa/K
# *
# * Features: This function calculates the temperature and pressure 
# *           dependent psychrometric constant
# *
# *           Depends on:
# *           - enthalpy_vap() . calculates latent heat of vaporization
# *
# * Ref: Temperature dependency of psychrometric constant (Eq. 8):
# *      Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998), 
# *        'Meteorological data,' Crop evapotranspiration - Guidelines 
# *        for computing crop water requirements - FAO Irrigation and 
# *        drainage paper 56, Food and Agriculture Organization of the 
# *        United Nations, Available: 
# *        http://www.fao.org/docrep/x0490e/x0490e07.htm
# *
# *     Pressure dependency of psychrometric constant (i.e., Cp):
# *     Tsilingris (2008), Thermophysical and transport properties of 
# *       humid air at temperature range between 0 and 100 °C, Energy 
# *       Conversion and Management, vol. 49, pp. 1098--1110.
# ************************************************************************
psychro <- function(tc, pa){
  # Calculate the specific heat capacity of water, J/kg/K
  cp <- 1.0045714270 +
    (2.050632750e-3)*tc -
    (1.631537093e-4)*tc*tc +
    (6.212300300e-6)*tc*tc*tc -
    (8.830478888e-8)*tc*tc*tc*tc +
    (5.071307038e-10)*tc*tc*tc*tc*tc
  cp <- (1e3)*cp
  #
  # # Calculate latent heat of vaporization, J/kg
  lv <- enthalpy_vap(tc)
  # 
  # # Calculate psychrometric constant, Pa/K
  return(cp*kMa*pa/(kMv*lv))
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define constants #########################################################
# /////////////////////////////////////////////////////////////////////////////
# NOTE: orbital parameters: eccentricity, obliquity, and longitude of the 
# perihelion, are assumed constant while they infact vary slightly over time.
# There are methods for their calculation (e.g., Meeus, 1991). Eccentricity 
# varies 0.005--0.072 and is decreasing at rate of 0.00004 per century. 
# Obliquity varies 22.1--24.5 degrees with a period of ~41000 years.
# 
kA <- 107             # constant for Rl (Monteith & Unsworth, 1990)
kalb_sw <- 0.17       # shortwave albedo (Federer, 1968)
kalb_vis <- 0.03      # visible light albedo (Sellers, 1985)
kb <- 0.20            # constant for Rl (Linacre, 1968; Kramer, 1957)
kc <- 0.25            # constant for Rs (Linacre, 1968)
kCw <- 1.05           # supply constant, mm/hr (Federer, 1982)
kd <- 0.50            # constant for Rs (Linacre, 1968)
ke <- 0.01670         # eccentricity of earth's orbit, 2000CE (Berger 1978)
keps <- 23.44         # obliquity of earth's elliptic, 2000CE (Berger 1978)
kfFEC <- 2.04         # from-flux-to-energy, umol/J (Meek et al., 1984)
kG <- 9.80665         # gravitational acceleration, m/s^2 (Allen, 1973)
kGsc <- 1360.8        # solar constant, W/m^2 (Kopp & Lean, 2011)
kL <- 0.0065          # adiabatic lapse rate, K/m (Cavcar, 2000)
kMa <- 0.028963       # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kMv <- 0.01802        # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
komega <- 283         # lon. of perihelion, degrees, 2000CE (Berger, 1978)
kSecInDay <- 86400    # number of seconds in a day
kPo <- 101325         # standard atmosphere, Pa (Allen, 1973)
kR <- 8.3143          # universal gas constant, J/mol/K (Allen, 1973)
kTo <- 298.15         # base temperature (25 deg C), K (I.C. Prentice)
kWm <- 150            # soil moisture capacity, mm (Cramer-Prentice, 1988)
kw <- 0.26            # PET entrainment, (1+kw)*EET (Priestley-Taylor, 1972)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### STASH 2.0 ################################################################
# /////////////////////////////////////////////////////////////////////////////
if (0){
  my_evap <- evap(51.4, 172, 74, 2001, 0.43, 17.3, 0.5)
}

# Initialize daily results:
daily_totals <- matrix(data=rep(0, 3660), nrow=366, ncol=10)
daily_totals <- as.data.frame(daily_totals)
names(daily_totals) <- c('ho',   # daily solar irradiation, J/m2 
                         'hn',   # daily net radiation, J/m2
                         'qn',   # daily PPFD, mol/m2
                         'cn',   # daily condensation, mm
                         'wn',   # daily soil moisture, mm
                         'ro',   # daily runoff, mm
                         'eq_n', # daily equilibrium ET, mm
                         'ep_n', # daily potential ET, mm
                         'ea_n') # daily actual ET, mm

# Initialize monthly results:
monthly_totals <- matrix(data=rep(0, 72), nrow=12, ncol=6)
monthly_totals <- as.data.frame(monthly_totals)
names(monthly_totals) <- c("eq_m",  # monthly equilibrium ET, mm 
                           "ep_m",  # monthly potential ET, mm
                           "ea_m",  # monthly actual ET, mm 
                           "cpa",   # Cramer-Prentice alpha, unitless 
                           "cwd",   # climatic water deficit, mm 
                           "q_m")   # monthly PPFD, mol/m2

# Location constants:
my_lat <- 37.7
my_elv <- 142

# Calculate days in the year
y <- 2000
ny <- julian_day(y+1, 1, 1) - julian_day(y, 1, 1)

# Example data (San Francisco, 2000 CE)
#   $sf, fractional sunshine hours, unitless
#   $tair, air temperature, deg C
#   $pn, daily precipitation, mm
my_file <- 'example_data.csv'
DATA <- read.csv(my_file)

# monthly:
all_months <- seq(from=1, to=12, by=1)
monthly_totals <- monthly_totals*0
for (m in all_months){
  # Calculate days of current month:
  nm <- julian_day(y, m+1, 1) - julian_day(y, m, 1)
  #
  # daily:
  for (i in seq(from=1, to=nm, by=1)){
    # Calculate day of year:
    n <- julian_day(y, m, i) - julian_day(y, 1 , 1) + 1
    #
    # Calculate evaporative supply (mm/hr)
    # * based on yesterday's available soil moisture
    #   ref: Federer (1982); Eq. 4, Prentice et al. (1993)
    idx <- (n - 1)
    if (idx < 1){
      idx <- ny
    }
    sw <- kCw*daily_totals$wn[idx]/kWm
    #
    # Compute daily radiation and evaporations values:
    ET <- evap(my_lat, n, my_elv, y, DATA$sf[n], DATA$tair[n], sw)
    #
    # Update daily soil moisture:
    daily_totals$wn[n] <- daily_totals$wn[idx] + DATA$pn[n] + ET$cond_mm - ET$aet_mm
    #
    if (daily_totals$wn[n] > kWm){
      # Bucket is full:
      # - set soil moisture to capacity
      # - add remaining water to runoff
      daily_totals$ro[n] <- daily_totals$wn[n] - kWm
      daily_totals$wn[n] <- kWm
    } else if (daily_totals$wn[n] < 0){
      # Bucket is empty:
      # - set runoff and soil moisture equal to zero
      daily_totals$ro[n] <- 0
      daily_totals$wn[n] <- 0
    } else{
      daily_totals$ro[n] <- 0
    }
    #
    if ((ny == 365) & (n == 365)){
      daily_totals$wn[n+1] <- daily_totals$wn[n]
    }
    #
    # Save daily results:
    daily_totals$ho[n] <- ET$ra_j.m2
    daily_totals$hn[n] <- ET$rn_j.m2
    daily_totals$qn[n] <- ET$ppfd_mol.m2
    daily_totals$cn[n] <- ET$cond_mm
    daily_totals$eq_n[n] <- ET$eet_mm
    daily_totals$ep_n[n] <- ET$pet_mm
    daily_totals$ea_n[n] <- ET$aet_mm
    #
    # Update monthly totals:
    monthly_totals$eq_m[m] <- monthly_totals$eq_m[m] + ET$eet_mm
    monthly_totals$ep_m[m] <- monthly_totals$ep_m[m] + ET$pet_mm
    monthly_totals$ea_m[m] <- monthly_totals$ea_m[m] + ET$aet_mm
    monthly_totals$q_m[m] <- monthly_totals$q_m[m] + ET$ppfd_mol.m2
  } # end daily
  monthly_totals$cpa[m] <- monthly_totals$ea_m[m]/monthly_totals$eq_m[m]
  monthly_totals$cwd[m] <- monthly_totals$ep_m[m] - monthly_totals$ea_m[m]
} # end monthly

#### Save Results ####
daily_outfile <- '/home/user/Desktop/out/stash_results_daily.csv'
write.csv(daily_totals, file=daily_outfile)

monthly_outfile <- '/home/user/Desktop/out/stash_results_monthly.csv'
write.csv(monthly_totals, file=monthly_outfile)

#### View Results ####
##
## Plot monthly ET results
##
par(mfrow=c(1,1))
par(mar=c(4,4.5,1,1))
plot(monthly_totals$ep_m,type='l',lwd=2, col='black',
     ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
lines(monthly_totals$ea_m, lty=2, lwd=2, col='cyan')
lines(monthly_totals$eq_m, lty=1, lwd=2, col='green')
lines(monthly_totals$cwd, lty=3, lwd=2, col='red')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=1, to=12, by=1))
axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=1, to=12, by=1))
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, 'Month', line=2)
mtext(side=2, as.expression(ET~(mm)), line=3)
legend('topleft',legend=c("PET", "AET", "EET", "CWD"),
       col=c('black', 'cyan', 'green', 'red'),
       lty=c(1,2,1,3), lwd=c(2,2,2,2), inset=0.01,
       y.intersp=1.2, horiz=FALSE, bty='n')


##
## Plot monthly results
##
tiff(out_file, width=900, height=500, units='px', 
     compression='none', pointsize=18, res=72)

par(mfrow=c(4,1))
# [1]
par(mar=c(1,5,1,1))
plot(monthly_totals$ep_m,type='l',lwd=2, col='black',
     ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
lines(monthly_totals$ea_m, lty=2, lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-50, to=175, by=25))
axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-50, to=200, by=50), 
     cex.axis=1.6)
mtext(side=2, expression(italic(E[m])~(mm)), line=3, cex=1.1)
legend('topright', legend=c(expression(italic(E[m]^{p})),
                            expression(italic(E[m]^{a}))),
       col=c('black', 'black'), lty=c(1, 2), cex=1.6, inset=0.02,
       adj=c(0.5, 0.5), lwd=c(2, 2), horiz=TRUE, bty='n', seg.len=1)
text(0.6, 150, "(a)", pos=4, cex=1.6)
# [2]
par(mar=c(1,5,1,1))
plot(monthly_totals$cwd,type='l',lwd=2, col='black', 
     ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-50, to=175, by=25))
axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-50, to=150, by=50), 
     cex.axis=1.6)
mtext(side=2, expression(Delta*italic(E[m])~(mm)), line=3, cex=1.1)
text(0.6, 150, "(b)", pos=4, cex=1.6)
# [3]
par(mar=c(1,5,1,1))
plot(monthly_totals$eq_m,type='l',lwd=2, col='black',
     ylim=c(0, max(monthly_totals$ep_m)), xlab=NA, ylab=NA, axes=F)
lines(monthly_totals$ea_m, lty=2, lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-50, to=175, by=25))
axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-50, to=150, by=50), 
     cex.axis=1.6)
mtext(side=2, expression(italic(E[m])~(mm)), line=3, cex=1.1)
legend('topright', legend=c(expression(italic(E[m]^{q})),
                            expression(italic(E[m]^{a}))),
       col=c('black', 'black'), lty=c(1, 2), cex=1.6, inset=0.02,
       adj=c(0.5, 0.5), lwd=c(2, 2), horiz=TRUE, bty='n', seg.len=1)
text(0.6, 150, "(c)", pos=4, cex=1.6)
# [4]
par(mar=c(2,5,1,1))
plot(monthly_totals$cpa,type='l',lwd=2, col='black', 
     ylim=c(0, 1.3), xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=12, by=1))
axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=1, to=12, by=1),
     cex.axis=1.6)
axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=-0.3, to=1.2, by=0.3))
axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=-0.3, to=1.2, by=0.3), 
     cex.axis=1.6)
mtext(side=2, expression(alpha), line=3, cex=1.1)
text(0.6, 1.1, "(d)", pos=4, cex=1.6)

dev.off()


##
## Plot daily ET results
##
par(mfrow=c(1,1))
par(mar=c(2.5,4.5,1,1))
plot(daily_totals$ep_n, type='l', lwd=2, col='blue', 
     ylim=c(0,1.25*max(daily_totals$ep_n)),
     xlab=NA, ylab=NA, axes=F)
lines(daily_totals$eq_n, lty=1, lwd=2, col='green')
lines(daily_totals$ea_n, lty=2, lwd=2, col='red')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=360, by=60))
axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=0, to=360, by=60), cex.axis=1.2)
axis(side=2, las=1, tck=-0.02, labels=NA, at=seq(from=0, to=8, by=1))
axis(side=2, las=1, lwd=0, line=-0.4, at=seq(from=0, to=8, by=1), cex.axis=1.2)
mtext(side=2, expression(list(Evapotranspiration, mm~d^{-1})), line=2, cex=1.2)
legend('top',legend=c("Potential", 
                      "Equilibrium   ", 
                      "Actual"),
       col=c('blue', 'green', 'red'),
       lty=c(1,1,2), lwd=c(2,2,2), inset=0.01, x.intersp=1.1,
       y.intersp=2.0, horiz=TRUE, bty='n', cex=1.2)


##
## Plot daily results
##
tiff(out_file, width=900, height=1000, units='px', 
     compression='none', pointsize=16, res=72)

par(mfrow=c(8,1))
# [1]
par(mar=c(1,5,1,1))
plot(DATA$sf, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0.3, to=0.7, by=0.1))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6, 
     at=seq(from=0.3, to=0.7, by=0.1))
mtext(side=2, expression(italic(S[f])), line=3, cex=1.1)
text(-12, 0.65, "(a)", pos=4, cex=1.7)
# [2]
par(mar=c(1,5,1,1))
plot(1e-6*daily_totals$hn, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=3, to=18, by=3))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
     at=seq(from=3, to=18, by=3))
mtext(side=2, expression(italic(H[N])~(MJ~m^{-2})), line=3, cex=1.1)
text(-12, 17, "(b)", pos=4, cex=1.7)
# [3]
par(mar=c(1,5,1,1))
plot(daily_totals$cn, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0.4, to=0.8, by=0.1))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6, 
     at=seq(from=0.4, to=0.8, by=0.1))
mtext(side=2, expression(italic(C[n])~(mm)), line=3, cex=1.1)
text(-12, 0.75, "(c)", pos=4, cex=1.7)
# [4]
par(mar=c(1,5,1,1))
plot(DATA$pn, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=-5, to=25, by=5))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6, 
     at=seq(from=-5, to=25, by=5))
mtext(side=2, expression(italic(P[n])~(mm)), line=3, cex=1.1)
text(-12, 22, "(d)", pos=4, cex=1.7)
# [5]
par(mar=c(1,5,1,1))
plot(daily_totals$wn, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0, to=150, by=30))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
     , at=seq(from=0, to=150, by=30))
mtext(side=2, expression(italic(W[n])~(mm)), line=3, cex=1.1)
text(-12, 130, "(e)", pos=4, cex=1.7)
# [6]
par(mar=c(1,5,1,1))
plot(daily_totals$ro, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=-5, to=20, by=5))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
     at=seq(from=-5, to=20, by=5))
mtext(side=2, expression(italic(RO)~(mm)), line=3, cex=1.1)
text(-12, 17, "(f)", pos=4, cex=1.7)
# [7]
par(mar=c(1,5,1,1))
plot(DATA$tair, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=0, to=25, by=5))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6, 
     at=seq(from=0, to=25, by=5))
mtext(side=2, expression(italic(T[air])~(degree*C)), line=3, cex=1.1)
text(-12, 23, "(g)", pos=4, cex=1.7)
# [8]
par(mar=c(2,5,1,1))
plot(daily_totals$ep_n, type='l', lwd=2, xlab=NA, ylab=NA, axes=F,
     ylim=c(0, max(daily_totals$ep_n)))
lines(daily_totals$ea_n, lty=2, lwd=2)
axis(side=1, las=1, tck=-0.03, labels=NA, at=seq(from=-60, to=720, by=60))
axis(side=1, las=1, lwd=0, line=-0.4, at=seq(from=-60, to=720, by=60), 
     cex.axis=1.6)
axis(side=2, las=1, tck=-0.03, labels=NA, at=seq(from=-1.5, to=6, by=1.5))
axis(side=2, las=1, lwd=0, line=-0.4, cex.axis=1.6,
     at=seq(from=-1.5, to=6, by=1.5))
mtext(side=2, expression(italic(E[n])~(mm)), line=3, cex=1.1)
text(-12, 5, "(h)", pos=4, cex=1.7)

dev.off()

#### Example Yearly Loop ####
data_dir <- '/home/user/Projects/STASH/'
my_files <- list.files(data_dir, pattern='^stash_data_.*csv$')
for (f in sort(my_files)){
  #
  my_file <- paste(data_dir, f, sep='')
  DATA <- read.csv(my_file)
  # 
  # Extract year from file name (the only numbers in the file name!)
  y <- as.numeric(gsub("[^0-9]", "", f))
  ny <- julian_day(y+1, 1, 1) - julian_day(y, 1, 1)
  #
  monthly_totals <- monthly_totals*0
  for (m in all_months){
    # Calculate days of current month:
    nm <- julian_day(y, m+1, 1) - julian_day(y, m, 1)
    #
    # daily:
    for (i in seq(from=1, to=nm, by=1)){
      # Calculate day of year:
      n <- julian_day(y, m, i) - julian_day(y, 1 , 1) + 1
      #
      # Calculate evaporative supply (mm/hr)
      # * based on yesterday's available soil moisture
      #   ref: Federer (1982); Eq. 4, Prentice et al. (1993)
      idx <- (n - 1)
      if (idx < 1){
        idx <- ny
      }
      sw <- kCw*daily_totals$wn[idx]/kWm
      #
      # Compute daily radiation and evaporations values:
      ET <- evap(my_lat, n, my_elv, y, DATA$sf[n], DATA$tair[n], sw)
      #
      # Update daily soil moisture:
      daily_totals$wn[n] <- daily_totals$wn[idx] + DATA$pn[n] + ET$cond_mm - ET$aet_mm
      #
      if (daily_totals$wn[n] > kWm){
        # Bucket is full:
        # - set soil moisture to capacity
        # - add remaining water to runoff
        daily_totals$ro[n] <- daily_totals$wn[n] - kWm
        daily_totals$wn[n] <- kWm
      } else if (daily_totals$wn[n] < 0){
        # Bucket is empty:
        # - set runoff and soil moisture equal to zero
        daily_totals$ro[n] <- 0
        daily_totals$wn[n] <- 0
      } else{
        daily_totals$ro[n] <- 0
      }
      #
      if ((ny == 365) & (n == 365)){
        daily_totals$wn[n+1] <- daily_totals$wn[n]
      }
      #
      # Save daily results:
      daily_totals$ho[n] <- ET$ra_j.m2
      daily_totals$hn[n] <- ET$rn_j.m2
      daily_totals$qn[n] <- ET$ppfd_mol.m2
      daily_totals$cn[n] <- ET$cond_mm
      daily_totals$eq_n[n] <- ET$eet_mm
      daily_totals$ep_n[n] <- ET$pet_mm
      daily_totals$ea_n[n] <- ET$aet_mm
      #
      # Update monthly totals:
      monthly_totals$eq_m[m] <- monthly_totals$eq_m[m] + ET$eet_mm
      monthly_totals$ep_m[m] <- monthly_totals$ep_m[m] + ET$pet_mm
      monthly_totals$ea_m[m] <- monthly_totals$ea_m[m] + ET$aet_mm
      monthly_totals$q_m[m] <- monthly_totals$q_m[m] + ET$ppfd_mol.m2
    } # end daily
    monthly_totals$cpa[m] <- monthly_totals$ea_m[m]/monthly_totals$eq_m[m]
    monthly_totals$cwd[m] <- monthly_totals$ep_m[m] - monthly_totals$ea_m[m]
  } # end monthly
  #
  # Save results to file:
  daily_outfile <- paste('/home/user/Desktop/out/stash_results_', 
                         as.character(y), '-daily.csv', sep='')
  write.csv(daily_totals, file=daily_outfile)
  monthly_outfile <- paste('/home/user/Desktop/out/stash_results_', 
                           as.character(y), '-monthly.csv', sep='')
  write.csv(monthly_totals, file=monthly_outfile)
}
