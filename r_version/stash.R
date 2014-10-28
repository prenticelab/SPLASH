# R-Studio
#
# stash.R
#
# written by Tyler W. Davis
# Imperial College London
#
# last updated: 2014-10-28
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script runs the STASH 2.0 code for point-based data.
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: evap
# *
# * Input: double, longitude, degrees (user.lon)
# *        double, latitude, degrees (user.lat)
# *        double, day of year (user.day)
# *        double, year (user.year) *optional
# *        double, fraction of sunshine hours (user.fsun)
# *        double, air temperature, deg C (user.tair)
# *        double, elevation (user.elv)
# *        double, soil moisture supply, mm/hr (user.supply)
# *        char,   distance method (user.dr)
# *                -> loutre -- (1/pd)^2
# *                -> klein -- 1+2e cos(2pin/N)
# *        char,   lambda method (user.lambda)
# *                -> kepler -- Simplified Kepler Method
# *                -> woolf -- Woolf 1968 method
# *                -> berger -- Berger'78 method
# *        char,   delta method (user.delta)
# *                -> loutre -- asin(sin(lambda)*sin(e))
# *                -> circle -- -eps cos(2pi/N*(n+10))
# *                -> cooper -- eps sin(2pi/N*(omega+n))
# *                -> spencer -- Fourier series
# *
# * Return: list object (et.srad)
# *         $time_hr ........... half-hourly timestamps, hr
# *         $dr ................ Ra distance factor, unitless
# *         $delta_deg ......... Declination angle, degrees
# *         $eot_min ........... Equation of time, min
# *         $lc_hr ............. Longitudinal correction factor, hr
# *         $hs_deg ............ Sunset angle, degrees
# *         $hs_lct ............ Sunset in local clock time, hr
# *         $ds_hr ............. Daylight hours, hr
# *         $dn_hr ............. Hours of positive Rn, hr
# *         $ra_w.m2 ........... Half-hourly extraterr. solar rad, W/m^2
# *         $rs_w.m2 ........... Half-hourly shortwave solar rad, W/m^2
# *         $rl_w.m2 ........... Constant daily longwave rad, W/m^2
# *         $rn_w.m2 ........... Half-hourly net radiation, W/m^2
# *         $hn_deg ............ Rn hour angle, degrees
# *         $rn_j.m2 ........... Daily net radiation, J/m^2
# *         $ppfd_umol.s.m2 .... Half-hourly PPFD, umol/s/m^2
# *         $ppfd_mol.m2 ....... Daily PPFD, mol/m^2
# *         $econ_m3.j ......... Water to energy conversion, m^3/J
# *         $cond_mm ........... Daily condensation, mm
# *         $eet_mm.hr ......... Half-hourly equilibrium ET, mm/hr
# *         $eet_mm ............ Daily equilibrium ET, mm
# *         $pet_mm.hr ......... Half-hourly potential ET, mm/hr
# *         $di_mm ............. Daily water demand, mm
# *         $si_mm ............. Daily water supply, mm
# *         $aet_mm.hr ......... Half-hourly actual ET, mm/hr
# *         $aet_mm ............ Daily actual ET, mm
# *
# * Features: This function calculates daily evaporation rates.
# *
# *           Depends on:
# *           - ka ............. length of earth's semi-major axis
# *           - kalb_sw ........ shortwave albedo
# *           - kalb_vis ....... visible light albedo
# *           - kb ............. empirical constant for longwave rad
# *           - kc ............. empirical constant for shortwave rad
# *           - kd ............. empirical constant for shortwave rad
# *           - ke ............. eccentricity
# *           - keps ........... obliquity
# *           - kfFEC .......... from-flux-to-energy conversion, umol/J
# *           - kGsc ........... solar constant
# *           - kGM ............ standard gravity of the sun
# *           - kZ ............. height of atmosphere, m
# *           - kw ............. entrainment factor for PET
# *           - dcos() ......... cos(x*pi/180), where x is in degrees
# *           - dsin() ......... sin(x*pi/180), where x is in degrees
# *           - dtan() ......... tan(x*pi/180), where x is in degrees
# *           - density_h2o() .. density of water
# *           - elv2pres() ..... elevation dependent atm. pressure
# *           - enthalpy_vap() . latent heat of vaporization
# *           - julian_day() ... date to julian day
# *           - map_days() ..... day to orbital longitude
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
# * Cooper, P. I. (1969). “The absorption of radiation in solar stills”. 
# *   In: Solar Energy 12.3, pp. 333–346 
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
# * Loutre, M.F. (2003) Ice ages (Milankovitch theory), Encyclopedia of 
# *   Atmospheric Sciences, edited by J.R. Holton, J.A. Curry, and J.A. 
# *   Pyle, pp. 995–1003, Elsevier Ltd.
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
evap <- function(user.lon, user.lat, user.day, user.year=0, user.fsun=1, 
                 user.tair=23.0, user.elv=0, user.supply=1.0,
                 user.dr='loutre', user.lambda='kepler', user.delta='loutre'){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (user.lon > 180 || user.lon < -180){
    stop("Warning: Longitude outside range of validity (-180 to 180)!")
  }
  if (user.lat > 90 || user.lat < -90){
    stop("Warning: Latitude outside range of validity (-90 to 90)!")
  }
  if (user.day < 1 || user.day > 366){
    stop("Warning: Day outside range of validity (1 to 366)!")
  }
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  et.srad <- list()
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # The number of days in the year:
  if (user.year == 0){
    kN = 365
  } else {
    kN = (julian_day(user.year+1, 1, 1) - julian_day(user.year, 1, 1))
  }
  #
  # Half-hourly time series:
  local.hh <- (seq(48) - 1.0)*0.5
  et.srad$time_hr <- local.hh
  #
  # One-second time series (for calculating total daylight hours):
  local.sec <- (seq(86400) - 1.0)/3600.0
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Heliocentric longitudes: nu and lambda (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (user.lambda == 'kepler'){
    # Kepler method:
    mdoy <- map_days(user.day, user.year, kN)
    my_nu<- mdoy[1]
    my_lambda<- mdoy[2]
  } else if (user.lambda == 'woolf'){
    # Woolf method:
    B <- (user.day-1)*360/kN
    my_lambda <- 279.9348 + B + 1.914827*dsin(B) - 0.079525*dcos(B) + 
      0.019938*dsin(2*B) - 0.00162*dcos(2*B)
    if (my_lambda < 0){
      my_lambda <- my_lambda + 360
    } else if (my_lambda > 360) {
      my_lambda <- my_lambda - 360
    }
    my_nu <- my_lambda - komega
    if (my_nu < 0){
      my_nu <- my_nu + 360
    }
  } else if (user.lambda == 'berger') {
    # Berger'78 method:
    pir <- pi/180.0
    xse <- sqrt(1 - ke^2)
    xee <- ke^2
    xec <- ke^3
    xlam <- (ke/2.0 + xec/8.0)*(1 + xse)*sin(komega*pir) - 
      xee/4.0*(0.5 + xse)*sin(2.0*komega*pir) + 
      xec/8.0*(1.0/3.0 + xse)*sin(3.0*komega*pir)
    xlam <- 2.0*xlam/pir
    dlamm <- xlam + (user.day - 80.0)*(360.0/kN)
    anm <- dlamm - komega
    ranm <- anm*pir
    ranv <- ranm + (2.0*ke - xec/4.0)*sin(ranm) +
      5.0/4.0*xee*sin(2.0*ranm) + 
      13.0/12.0*xec*sin(3.0*ranm)
    anv <- ranv/pir
    my_lambda <- anv + komega
    if (my_lambda < 0){
      my_lambda <- my_lambda + 360
    } else if (my_lambda > 360) {
      my_lambda <- my_lambda - 360
    }
    my_nu <- my_lambda - komega
    if (my_nu < 0){
      my_nu <- my_nu + 360
    }
  } else {
    stop("Warning: lambda method not recognized!")
  }
  et.srad$lambda_deg <- my_lambda
  et.srad$nu_deg <- my_nu
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate time zone hour (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Correction for local time zone; based on one hour per 15 degrees longitude
  # * positive or negative deviation from UTC/GMT
  if (user.lon < 0){
    temp.lon <- -1.0*user.lon
    temp.tzh <- floor(temp.lon/15)
    tz.hour <- -1.0*temp.tzh
  } else {
    tz.hour <- floor(user.lon/15)
  }
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate distance factor (unitless)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (user.dr == 'loutre'){
    # ref: Eq. 1, Loutre (2003)
    rho <- (1 - ke^2)/(1 + ke*dcos(my_nu))
    dr <- (1/rho)^2
  } else if (user.dr == 'klein'){
    # Klein (1977) method:
    dr <- 1 + 2*ke*cos(2*pi*user.day/kN)
  } else {
    stop("Warning: distance factor method not recognized!")
  }
  et.srad$dr <- dr
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate the declination angle (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Delta ranges from -23.45 in the winter to 23.45 in the summer
  if (user.delta == 'loutre'){
    # Loutre (2002)
    delta <- asin(dsin(my_lambda)*dsin(keps))*180/pi
  } else if (user.delta == 'cooper'){
    # Cooper method:
    delta <- keps*sin(2*pi*(komega + user.day)/kN)
  } else if (user.delta == 'circle'){
    # Circle method:
    delta <- -keps*cos(2*pi*(user.day + 10)/kN) 
  } else if (user.delta == 'spencer') {
    # Spencer method:
    B <- 2*pi*(user.day - 1)/kN
    delta <- 180.0/pi*(0.006918 - 0.399912*cos(B) + 0.070257*sin(B) - 
                         0.006758*cos(2*B) + 0.000907*sin(2*B) - 
                         0.0022697*cos(3*B) + 0.00148*sin(3*B))
  } else {
    stop("Warning: declination angle method not recognized!")
  }
  et.srad$delta_deg <- delta
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate the equation of time (minutes)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The equation of time calculates the difference between apparent solar time 
  # and mean solar time
  # ref: equation by Spencer (1971); corrected by Oglesbay (1998); unit 
  #      conversion by Iqbal (1983)
  B <- 2*pi*(user.day - 1)/kN
  EOT <- 4*180/pi*(
    7.5e-6 + 
      1.868e-3*cos(B) - 
      3.2077e-2*sin(B) - 
      1.4615e-2*cos(2*B) - 
      4.0849e-2*sin(2*B)
  )
  et.srad$eot_min <- EOT
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Longitudinal correction factor (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Correction factor for longitude
  # ref: Eq. 3.6, Stine and Geyer (2001)
  LC <- (15*tz.hour - user.lon)/15
  et.srad$lc_hr <- LC
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 07. Calculate solar time (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the solar time based on the local time
  # ref: Eq. 3.5, Stine and Geyer (2001)
  ts.hh <- local.hh + EOT/60.0 - LC
  ts.sec <- local.sec + EOT/60.0 - LC
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 08. Calculate the hour angle (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Solar zenith angle
  # ref: Eq. 3.1, Stine and Geyer (2001)
  # * solar noon is at 0 degrees (i.e., vertical line)
  # * based once again on 15 deg movement per hour
  w.hh <- 15*(ts.hh - 12)
  w.sec <- 15*(ts.sec - 12)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 09. Calculate the sunset angle (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 3.22, Stein & Geyer (2001)
  # - represents the number of hours from solar noon when the sun sets (i.e., 
  #   the Ra curve dips below zero)
  # - when using Hs to integrate over daylight hours, multiply by 2
  if (dtan(delta)*dtan(user.lat) >= 1.0){
    Hs <- 180  # Polar day (no sunset)
  } else if (dtan(delta)*dtan(user.lat) <= -1.0){ 
    Hs <- 0 # Polar night (no sunrise)
  } else {
    Hs <- (180/pi)*acos(-dtan(delta)*dtan(user.lat))
  }
  et.srad$hs_deg <- Hs
  et.srad$hs_lct <- (Hs/15.0) - (EOT/60.0) + LC + 12   # local clock time (hr)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 10. Calculate daylight hours (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 34, Allen et al. (1998)
  ds <- (24.0/180*Hs)
  et.srad$ds_hr <- (ds)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 11. Calculate the inclination factor (unitless)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Inclination for sub-daily solar radiation
  # ref: Eq. 1, Wetherald & Manabe (1972); 
  #      Eq. 1.6.5, Duffy & Beckman (1991); 
  #      Eq. 3, Loutre (2003)
  cosz.hh <- dcos(user.lat)*dcos(delta)*dcos(w.hh) + dsin(user.lat)*dsin(delta)
  cosz.sec <- (
    dcos(user.lat)*dcos(delta)*dcos(w.sec) + dsin(user.lat)*dsin(delta)
  )
  #
  # Inclination for daily solar radiation
  # ref: Eq. 1.10.3, Duffy & Beckman (1991); Eq. 5, Loutre (2003)
  # * assumes declination angle is constant throughout the day 
  #   (max rate of change is ~0.5 degree/day)
  cosz.day <- (
    dsin(delta)*dsin(user.lat)*(pi/180)*Hs + 
      dcos(delta)*dcos(user.lat)*dsin(Hs)
  )
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 12. Calculate Ra (W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extraterrestrial solar radiation w.r.t. a flat horizon
  # ref: Eq. 1.10.2 in Duffie and Beckman (1991)
  Ra.hh <- kGsc*dr*cosz.hh
  Ra.sec <- kGsc*dr*cosz.sec
  #
  # Set negative values equal to zero (i.e., no sun):
  Ra.hh[Ra.hh < 0 ] <- 0
  Ra.sec[Ra.sec < 0] <- 0
  #
  # Save half-hourly time series to output list:
  et.srad$ra_w.m2 <- (Ra.hh)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 13. Calculate atmospheric transmittivity (unitless)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Atmospheric transmittivity is an elevation-adjusted mean sea-level value,
  # such that, tau = tau_o * exp(z/H), where tau_o is the transmittivity at
  # mean sea level and exp(z/H) is a pseudo-Beer's law adjustment factor for
  # elevation effects.
  # ref:  Eq. 11, Linacre (1968): tau_o
  #       Eq. 2, Allen (1996): tau
  tau_o <- (kc + kd*user.fsun)
  tau <- tau_o*(1+2.67e-5*user.elv)
  et.srad$tau <- tau
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 14. Calculate net shortwave radiation, Rns (W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Net shortwave radiation, Rns = (1-albedo)*Rs where Rs = tau*Ra
  # ref: Eq. 4, Linacre (1968); Eq. 7, Prentice et al. (1993)
  Rs.hh <- (1 - kalb_sw)*(tau*Ra.hh)
  Rs.sec <- (1 - kalb_sw)*(tau*Ra.sec)
  et.srad$rs_w.m2 <- (Rs.hh)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 15. Calculate net longwave radiation, Rnl (W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Net longwave radiation (assumed constant throughout the day)
  # ref: Eq. 11, Prentice et al. (1993); Eq. 5 & 6, Linacre (1968)
  # - provides an upper limit value of Rnl due to the assumption that soil 
  #   and air temperature are equal (i.e., well-water soil conditions)
  # - dry soils would result in appreciably smaller Rnl
  Rl <- (kb + (1.0 - kb)*user.fsun)*(kA - user.tair)
  et.srad$rl_w.m2 <- (Rl)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 16. Calculate net radiation, Rn (W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Net radiation, Rn = (Rns - Rnl)
  # ref: Eq. 4, Linacre (1968)
  Rn.hh <- (Rs.hh - Rl)
  Rn.sec <- (Rs.sec - Rl)
  et.srad$rn_w.m2 <- Rn.hh
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 17. Calculate daily Ra (J/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 1.10.3, Duffy & Beckman (1991)
  et.srad$ra_j.m2 <- (86400.0/pi)*kGsc*dr*cosz.day
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 18. Calculate Rn cross-over angle, Hn (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The cross-over angle is based on the hour when Rn goes from negative to
  # positive. The following is the derivation for the cross-over angle, Hn. 
  # ref: original STASH FORTRAN code in-line comments (LBX-bern)
  #
  # Net radiation: 
  #   Rn = (Rns - Rnl). 
  # Substitution for Rns:
  #   Rn = (1-albedo)*tau*Ra - Rnl
  # Substitution for Ra:
  #   Rn = (1-albedo)*tau*Gsc*dr*cosz - Rl
  # Substitution for cosz:
  #   Rn = (1-albedo)*tau*Gsc*dr*(
  #        sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(w)) - Rl
  # Let 
  #   u = sin(phi)*sin(delta) 
  #   v = cos(phi)*cos(delta)
  #   w = (1-albedo)*tau*Gsc*dr
  # Such that:
  #   Rn = w*u + w*v*cos(w) - Rl
  # Set Rn=0 and solve for the angle, w
  #   cos(w(Rn=0)) = (-w*u + Rl)/(w*v)
  #   Hn = w(Rn=0) = arccos((-w*u + Rl)/(w*v))
  #
  # Variable substitution:
  u <- dsin(user.lat)*dsin(delta)                # unitless
  v <- dcos(user.lat)*dcos(delta)                # unitless
  w <- (1 - kalb_sw)*tau*kGsc*dr                    # W/m^2
  #
  # Check for exceedance similar to sunset angle:
  # NOTE: tan(delta)*tan(phi) == u/v
  if ((Rl - w*u)/(w*v) >= 1.0){
    Hn <- 0  # Net radiation is negative all day
  } else if ((Rl - w*u)/(w*v) <= -1.0){ 
    Hn <- 180 # Net radiation is positive all day
  } else {
    Hn <- (180/pi)*acos((Rl - w*u)/(w*v))
  }
  et.srad$hn_deg <- Hn
  et.srad$hn_lct <- (Hn/15.0) - (EOT/60.0) + LC + 12   # local clock time (hr)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 19. Calculate hours of positive Rn (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Use methodology similar to daylight hours
  dn <- (24.0/180*Hn)
  et.srad$dn_hr <- dn
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 20. Integrate daily Rn (J/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Integrate Rn = w*u + w*v*cos(w) - Rl
  #   Rn_daily = 2*(w*u*Hn + w*v*sin(Hn) - Rl*Hn)
  # - multiply by two (for whole day)
  # - include unit conversion factor: 86400/(2*pi)
  et.srad$rn_j.m2 <- (86400.0/pi)*((w*u - Rl)*(pi/180)*Hn + w*v*dsin(Hn))
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 21. Integrate nighttime Rn (J/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # There are two parts to the nighttime Rn curve:
  # - Rn curve from Hn to Hs
  # - Rl curve from Hs to pi
  # Integrate the piecewise nighttime Rn curve, Rnn:
  #   Rnn = [(sin(Hs)-sin(Hn))*v + (Hs-Hn)*u]*w + (pi+Hn-2*Hs)*Rl
  # - multiply integral by 2 (for whole day)
  # - multiply integral by unit factor: 86400/(2*pi)
  et.srad$rnn_j.m2 <- (86400.0/pi)*(
    ((dsin(Hs) - dsin(Hn))*v + (Hs - Hn)*(pi/180)*u)*w + 
      (pi + Hn*(pi/180) - 2*Hs*(pi/180))*Rl
  )
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 22. Integrate daily PPFD (mol/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Convert PAR (visible light) to PPFD:
  # ref: Ge et al. (2011)
  PPFD.hh <- kfFEC*(1 - kalb_vis)*(tau*Ra.hh)         # umol/s/m^2
  PPFD.sec <- kfFEC*(1 - kalb_vis)*(tau*Ra.sec)       # umol/s/m^2
  et.srad$ppfd_umol.s.m2 <- PPFD.hh
  #
  # To integrate, use methodology for daily Ra:
  # * 1.0e-6 (converts umol to mol), and unit conversion for seconds
  PPFD.day <- (86400.0/pi)*kfFEC*(1 - kalb_vis)*tau*kGsc*dr*cosz.day*1.0e-6
  et.srad$ppfd_mol.m2 <- PPFD.day
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 23. Calculate energy conversion factor of water (m^3/J)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate temperature & pressure dependent variables
  s <- sat_slope(user.tair)                         # Pa/K
  lv <- enthalpy_vap(user.tair)                     # J/kg
  pw <- density_h2o(user.tair, elv2pres(user.elv))  # kg/m^3
  g <- psychro(user.tair, elv2pres(user.elv))       # Pa/K
  #
  # Water to energy conversion factor:
  econ <- s/(lv*pw*(s + g))                         # m^3/J
  et.srad$econ_m3.j <- econ
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 24. Calculate condensation (mm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Assume condensation is equivalent to negative net radiation flux
  # - convert from energy to water, econ
  # - convert from m to mm (x1000.0)
  et.srad$cond_mm <- abs(et.srad$rnn_j.m2)*econ*1000
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 25. Calculate equilibrium evapotranspiration (mm/hr)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 5, Prentice et al. (1993)
  EET.hh <- econ*Rn.hh         # m/s
  EET.sec <- econ*Rn.sec       # m/s
  #
  # ET cannot be negative
  EET.hh[which(EET.hh < 0)] <- 0
  EET.sec[which(EET.sec < 0)] <- 0
  et.srad$eet_mm.hr <- (EET.hh*1000.0*3600.0)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 26. Integrate daily equilibrium evapotranspiration (mm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Apply unit conversion to the integrated daily Rn
  # - convert to units of water, econ
  # - convert from m to mm, x1000
  et.srad$eet_mm <- (1000.0*econ)*et.srad$rn_j.m2
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 27. Calculate potential evapotranspiration (mm/hr)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 14, Priestley & Taylor (1972)
  PET.hh <- (1.0 + kw)*EET.hh    # m/s
  PET.sec <- (1.0 + kw)*EET.sec  # m/s
  et.srad$pet_mm.hr <- (PET.hh*1000.0*3600.0)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 28. Calculate daily demand (mm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Apply scaling factor to integrated daily EET:
  et.srad$di_mm <- (1.0 + kw)*et.srad$eet_mm
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 29. Calculate daily supply (mm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Make certain supply is not negative:
  if (user.supply < 0){
    Sw <- 0.0
  } else {
    Sw <- user.supply
  }
  #
  # Daily supply assuming full sine curve
  # ref: Eq. 18b, Federer (1982)
  et.srad$si_mm <- (ds*Sw)
  #  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 30. Estimate actual evapotranspiration (mm/hr)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The instantaneous minimum of supply and demand: AET=min{Sw,Dp}
  # ref: Eq. 3, Prentice et al. (1993)
  # 
  # Convert from mm/hr to m/s:
  Sw <- user.supply*(1/1000.0)*(1/3600.0)
  #
  AET.sec <- PET.sec
  AET.sec[which(AET.sec > Sw)] <- Sw
  #
  AET.hh <- PET.hh
  AET.hh[which(AET.hh > Sw)] <- Sw
  et.srad$aet_mm.hr <- (AET.hh*1000.0*3600.0)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 31. Calculate the PET--Sw cross-over angle, degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The angle where PET == Sw is found by setting the two equations equal to
  # one another and solving for the hour angle
  # NOTE: Sw and PET are both in m/s
  #   PET = (1 + kw)*econ*(w*(u + v*cos(H)) - Rl) = Sw
  # Such that:
  #   cos(H=Hi) = Sw/(w*v*econ*(1 + kw)) + Rl/(w*v) - u/v
  #   Hi = arccos(Sw/(w*v*econ*(1 + kw)) + Rl/(w*v) - u/v)
  # As Sw -> 0, Hi -> Hn
  cosHi <- Sw/(w*v*econ*(1 + kw)) + Rl/(w*v) - u/v
  if (cosHi >= 1.0){
    Hi <- 0.0       # supply exceeds demand
  } else if (cosHi <= -1.0){
    Hi <- 180.0     # supply limits demand everywhere
  } else {
    Hi <- acos(cosHi)*(180/pi)
  }
  et.srad$hi_deg <- (Hi)
  et.srad$hi_lct <- (Hi/15) - (EOT/60) + LC + 12
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 32. Integrate daily actual evapotranspiration (mm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The AET curve consists of two parts:
  # - Sw from 0 to Hi
  # - PET from Hi to Hn
  # Integrate the piece-wise curve:
  #   AET = Sw*Hi + (1+kw)*econ[
  #           ((sin(Hn)-sin(Hi))*v + (Hn-Hi)*u)*w + (Hi-Hn)*Rl]
  # - multiply integral by 2 (for whole day)
  # - multiply by unit conversion: 86400/(2*pi)
  # - convert m to mm: x1000
  et.srad$aet_mm <- (86400/pi)*1000.0*(
    Sw*Hi*(pi/180) + (1 + kw)*econ*(
      ((dsin(Hn) - dsin(Hi))*v + (Hn-Hi)*(pi/180)*u)*w +
        (Hi - Hn)*(pi/180)*Rl
    )
  )
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  et.srad
}

# ************************************************************************
# * Name: earth_velocity
# *
# * Input: double or numeric list (lon), orbital lon. w.r.t. perihelion
# *        double (P), days in a year
# *
# * Return: double or numeric list, rad/day
# *
# * Features: This function calculates the angular velocity of earth
# *           on its orbit around the sun, rad/day.
# *
# *           Depends on:
# *           - ke ............. eccentricity of earth's orbit, unitless
# *
# * Ref: Kepler's Second Law of Planetary Motion
# ************************************************************************
earth_velocity <- function(lon, P=kP){
  # Kepler's Second Law:
  2*pi/P*(1+ke*dcos(lon))^2*(1-ke^2)^(-1.5)
}

# ************************************************************************
# * Name: julian_day
# *
# * Input: int (year) 
# *        int (month) 
# *        int (day)
# *
# * Return: int, Julian day
# *
# * Features: This function converts a date in the Gregorian calendar
# *           to a Julian day number (i.e., a method of consecutative 
# *           numbering of days---does not have anything to do with 
# *           the Julian calendar!)
# *
# * Ref: Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", 
# *      Astronomical Algorithms
# ************************************************************************
julian_day <- function(Y, M, D){
  if(M<=2){
    Y <- Y-1
    M <- M+12
  }
  A <- floor(Y/100)
  B <- 2 - A + floor(A/4)  # modified leap year def. for Gegorian calendar
                           # for Julian calendar, B = 0
  # Eq. 7.1, J. Meeus (1991), Astronomical Algorithms
  JDE <- 1.0*floor(365.25*(Y+4716)) + 
    1.0*floor(30.6001*(M+1)) + 1.0*D + 1.0*B - 1524.5
  JDE
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
# * Name: dtan
# *
# * Input: double (d), angle in degrees
# *
# * Return: double, tangent of angle
# *
# * Features: This function calculates the tangent of an angle (d) given
# *           in degrees.
# ************************************************************************
dtan <- function(d) {
  tan(d*pi/180)
}

# ************************************************************************
# * Name: equinox
# *
# * Input: double, year (year) 
# *        double (opt), flag, see below for description
# *
# * Return: double, day of the year
# *
# * Features: This function calculates the day of the year on which one 
# *           of the four season dates falls for a given year (-1000 to 
# *           +3000)---centuries outside this range should only have 
# *           slight errors. For season selection:
# *            opt=1 :: vernal equinox
# *            opt=2 :: summer solstice
# *            opt=3 :: autumnal equinox
# *            opt=4 :: winter solstice
# *
# * Depends: - dcos()
# *          - julian_day()
# *
# * Note: This script is based on the Javascript function written by
# *       C Johnson, Theoretical Physicist, Univ of Chicago
# *       URLs:
# *       1. 'Equation of Time': http://mb-soft.com/public3/equatime.html
# *       2. Javascript: http://mb-soft.com/believe/txx/astro22.js
# *
# * Ref: J. Meeus (1991), Chapter 26 "Equinoxes and Solstices", 
# *      Astronomical Algorithms
# ************************************************************************
equinox <- function(year, opt=1) {
  # Define variables:
  #   JDE0 :: instant of the 'mean' equinox or solstice
  #   JDE  :: Julian Ephemeris Day of equinox or solstice
  #   Y    :: year fraction (from Table 26.A / 26.B)
  #   JDE0tab :: Table 26.A or 26.B depending on year
  #   S       :: sum of 24 periodic terms in Table 26.C
  #   Tt      :: \
  #   W       :: -   Values to solve in J. Meeus, Ch. 26
  #   deltaL  :: / to solve for the equinox/solstice dates
  
  # Table 26.C, J. Meeus (1991), Astronomical Algorithms
  EquinoxpTerms <- c(
    # A    B          C
    485, 324.96,   1934.136,
    203, 337.23,  32964.467,
    199, 342.08,     20.186,
    182,  27.85, 445267.112,
    156,  73.14,  45036.886,
    136, 171.52,  22518.443,
    77, 222.54,  65928.934,
    74, 296.72,   3034.906,
    70, 243.58,   9037.513,
    58, 119.81,  33718.147,
    52, 297.17,    150.678,
    50,  21.02,   2281.226,
    45, 247.54,  29929.562,
    44, 325.15,  31555.956,
    29,  60.93,   4443.417,
    18, 155.12,  67555.328,
    17, 288.79,   4562.452,
    16, 198.04,  62894.029,
    14, 199.76,  31436.921,
    12,  95.39,  14577.848,
    12, 287.11,  31931.756,
    12, 320.81,  34777.259,
    9, 227.73,   1222.114,
    8,  15.45,  16859.074
  )
  
  # Table 26.A, J. Meeus (1991), Astronomical Algorithms
  JDE0tab1000 <- matrix(data=c(
    1721139.29189, 365242.13740,  0.06134,  0.00111, -0.00071, # March equinox
    1721233.25401, 365241.72562, -0.05323,  0.00907,  0.00025, # June solstice
    1721325.70455, 365242.49558, -0.11677, -0.00297,  0.00074, # September equinox
    1721414.39987, 365242.88257, -0.00769, -0.00933, -0.00006  # December solstice
  ),
  nrow=4,
  ncol=5,
  byrow=T
  )
  
  # Table 26.B, J. Meeus (1991), Astronomical Algorithms
  JDE0tab2000 <- matrix(data=c(
    2451623.80984, 365242.37404,  0.05169, -0.00411, -0.00057, # March equinox
    2451716.56767, 365241.62603,  0.00325,  0.00888, -0.00030, # June solstice
    2451810.21715, 365242.01767, -0.11575,  0.00337,  0.00078, # September equinox
    2451900.05952, 365242.74049, -0.06223, -0.00823,  0.00032  # December solstice
  ),
  nrow=4,
  ncol=5,
  byrow=T
  )
  
  if (year < 1000) {
    # Use Table 26.A for years -1000 to 1000
    JDE0tab = JDE0tab1000
    Y = year / 1000
  } else {
    # Use Table 26.B for years 1000 to 3000
    JDE0tab = JDE0tab2000
    Y = (year - 2000) / 1000
  }
  
  # Calculate the JDE0 term (from Table 26.A/26.B):
  JDE0 =  JDE0tab[opt,1] +
    (JDE0tab[opt,2] * Y) +
    (JDE0tab[opt,3] * Y * Y) +
    (JDE0tab[opt,4] * Y * Y * Y) +
    (JDE0tab[opt,5] * Y * Y * Y * Y)
  
  # Calculate the other three terms:
  Tt = (JDE0 - 2451545.0) / 36525
  W = (35999.373 * Tt) - 2.47
  deltaL = 1 + (0.0334 * dcos(W)) + (0.0007 * dcos(2*W))
  
  # Calculate the sum of the 24 periodic terms:
  S <- 0
  j <- 1
  for (i in seq(24)) {
    S <- S + EquinoxpTerms[j] * dcos(EquinoxpTerms[j+1] + (EquinoxpTerms[j+2] * Tt))
    j <- j + 3
  }
  
  # Calculate the JDE of the equinox/solstice:
  JDE = JDE0 + ((S * 0.00001) / deltaL)
  
  # Return vernal equinox day of year, i.e., DoY
  #(JDE - gregorian_to_jd(year, 1, 1) + 1)
  (JDE - julian_day(year, 1, 1) + 1)
}

# ************************************************************************
# * Name: simplified_kepler
# *
# * Input: double (lon), degrees
# *        double (P), days in year
# *
# * Return: double (t), days
# *
# * Features: This function calculates the amount of time (days) it 
# *           takes earth to travel a set amount of longitude with
# *           respect to the perihelion
# *
# *           "SIMPLIFIED KEPLER METHOD"
# *
# *           Depends on:
# *           - dsin() ......... dsin(x) = sin(x*pi/180)
# *           - dcos() ......... dcos(x) = cos(x*pi/180)
# *           - ke ............. eccentricity of earth's orbit
# *
# * Note: The original derivation of d(t)/d(lon) was computed in wxmaxima
# *
# *         d(t)/d(lon)= P/(2*pi)*(1-e^2)^1.5/(1+e*cos(lon))^2
# *
# *       The simplification assumes that eccentricity is small and so
# *       neglects the term in the numer and denom terms in the 
# *       full_kepler function.
# *
# *       The correction factor (2sin x) is due to periodic differences
# *       between the simple approach and the full expression presented
# *       in the full_kepler function
# *
# ************************************************************************
simplified_kepler <- function(lon, P=kP){
  P/(pi)*(atan(dsin(lon)/(dcos(lon)+1))-ke*dsin(lon))
}

# ************************************************************************
# * Name: correct_time
# *
# * Input: numeric list (x), longitude degrees
# *        numeric list (y), days
# *
# * Return: numeric list (y), days
# *
# * Features: This function corrects the series of days (y) for a given
# *           set of longitudes (x) calcuated by either full_kepler or 
# *           simplified_kepler functions. 
# *
# *           The problem with these two functions is that the atan()
# *           has an asymptote at pi/2 (180 degrees). In order to 
# *           correctly stack the series the slope of the last two 
# *           points before the asymptote is used to approximate the
# *           value of the first negative term (where it should be along
# *           the line, if it were unbroken). The negative terms are then
# *           shifted to this value to create a complete unbroken line,
# *           see figure below for a visual explanation:
# *
# *          +y                     +y
# *           |     o               |           o
# *           |   o                 |         o  
# *           | o                   |       o    
# *           |------------- +x  => |     o      
# *           |           o         |   o        
# *           |         o           | o          
# *           |       o             |------------ +x
# *
# *           (a) Original Data     (b) Corrected
# ************************************************************************
correct_time <- function(x,y){
  # Find indices for negative and positive values of time:
  neg_points <- which(y<0)
  pos_points <- which(y>=0)
  #
  # Extrapolation values for longitude:
  intrp_a <- pos_points[length(pos_points)]
  intrp_x <- neg_points[1]
  intrp_b <- pos_points[length(pos_points)-1]
  #
  intrp_A <- y[intrp_a]
  intrp_B <- y[intrp_b]
  #
  intrp_y <- intrp_A - (intrp_B - intrp_A)*(intrp_a-intrp_x)/(intrp_b-intrp_a)
  #
  # Calculate the difference between the interpolated value and actual
  diff_y <- intrp_y - y[intrp_x]
  #
  # Add the difference to the negative terms:
  y[neg_points]<- y[neg_points]+diff_y
  y
}

# ************************************************************************
# * Name: map_days
# *
# * Input: double (doy), day of year
# *        double (year)
# *        double (P), days in year
# *
# * Output: double, longitude w.r.t day of year
# *
# * Features: This script maps the earth's heliocentric longitude with
# *           the day of the year. The angle (komega) and calendar day of 
# *           the vernal equinox is used to map the time of earth's orbit to
# *           the correct calendar day. To interpolate the angle for a given
# *           day of the year, the closest known day, t, that precedes 
# *           the requested day, t+1, is moved forward through space, b,
# *           according to the angular velocity at that time w(a) 
# *           and the length of travel time:
# *
# *                             SUN
# *                      ------- O ------- +x
# *                       a    /   \
# *                           /  b  \
# *                          A       B          
# *                     EARTH(t)   EARTH(t+1)
# *
# *           Longitude of Earth(t) = a      <- known
# *           Longitude of Earth(t+1) = a+b  <- unknown
# *
# *           The angle traveled, b = w(a)*dt
# *              where w(a) is the angular velocity at longitude a
# *              and dt is the difference between t and t+1 (days)
# *
# *           Solving for longitude a+b = a + b => a + w(a)*dt
# *
# *           Depends on:
# *           - simplified_kepler() .... Simple Kepler Method, t(nu)
# *           - correct_time() ......... fixes sign change due to tan()
# *           - equinox() .............. calculates day of vernal equinox
# *           - earth_velocity() ....... calculates tangential velocity
# *           - komega ................. angle of the perihelion, degrees
# ************************************************************************
map_days <- function(doy, yr, P=kP){
  # Check doy is correct:
  if (doy < 1 | doy > 366){
    warning("DOY should be from 1 to 366")
  } else{
    # Create longitude field (nu, w.r.t. to perihelion)
    lon <- seq(0,360)
    # 
    # Calculate the time to circumnavigate the sun
    # * ranges between 0 (lon=0) and P (lon=360)
    kepler_days <- correct_time(lon,simplified_kepler(lon, P))
    #
    # Get angle between perihelion and the vernal equinox:
    wp <- (360-komega)
    #
    # Find how long it takes earth to travel to vernal equinox:
    if(wp<180){
      daysToVE <- simplified_kepler(wp, P)
    } else{
      # Set up longitudes before and after the jump:
      my_lons <- c(179,180,181,wp)
      daysToVE <- correct_time(my_lons,simplified_kepler(my_lons, P))[4]
    }
    #
    # Find which day the vernal equinox falls on:
    if (yr == 0){
      dayOfVE <- 80
    } else {
      dayOfVE <- equinox(yr)
    }
    #
    # Calculate the offset:
    dOffset <- (dayOfVE-daysToVE)
    #
    # Compute calendar for this years' orbital period:
    calendar <- (kepler_days + dOffset)
    calendar[which(calendar>=P)]<-calendar[which(calendar>=P)]-P
    #
    # Check to see if doy is in the calendar:
    if (any(doy == calendar)){
      # Get calendar index:
      icalendar <- which(calendar == doy)[1]
      #
      # Get nu directly:
      nu_doy <- lon[icalendar]
    } else {
      # Find the kepler day that preceeds the doy:
      dbefore <- sort(c(calendar,doy))[which(sort(c(calendar,doy))==doy)-1]
      #
      # Get index of this value:
      ibefore <- which(calendar==dbefore)[1]
      #
      # Get the angular velocity at this location, rad/day:
      vbefore <- earth_velocity(lon, P)[ibefore]
      #
      # Calculate the time between doy and the kepler day:
      tbetween <- (doy-dbefore)
      #
      # Calculate the orbit traveled during this time & convert to degrees:
      atraveled <- (vbefore*tbetween)*180/pi
      #
      # Longitude of doy (w.r.t perihelion)
      nu_doy <- lon[ibefore]+atraveled
    }
    # Longitude of doy (w.r.t. vernal equinox)
    lamda_doy <- (nu_doy + komega)
    if(lamda_doy >= 360){
      lamda_doy = lamda_doy-360
    }
    #
    c(nu_doy,lamda_doy)
  }
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
  kPo*(1-kL*z/kTo)^(kG*kMa/(kR*kL))
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
  (17.269)*(237.3)*(610.78)*exp(17.269*tc/(237.3+tc))/(237.3+tc)^2
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
  1.91846e6*((tc+273.15)/(tc+273.15-33.91))^2
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
  # Vo can be found by using do=(1/Vo) by Kell, (1975). 
  # The equation for do, the temperature-dependent density of water at 
  # 1 atm, g/cm^3
  do <- function(tc){
    C1 <- 0.99983952
    C2 <- 6.78826 * 10^(-5)
    C3 <- -9.08659 * 10^(-6)
    C4 <- 1.02213 * 10^(-7)
    C5 <- -1.35439 * 10^(-9)
    C6 <- 1.47115 * 10^(-11)
    C7 <- -1.11663 * 10^(-13)
    C8 <- 5.04407 * 10^(-16)
    C9 <- -1.00659 * 10^(-18)
    #
    d <- (C1 +
      C2 * tc +
      C3 * tc * tc +
      C4 * tc * tc * tc +
      C5 * tc * tc * tc * tc +
      C6 * tc * tc * tc * tc * tc +
      C7 * tc * tc * tc * tc * tc * tc +
      C8 * tc * tc * tc * tc * tc * tc * tc +
      C9 * tc * tc * tc * tc * tc * tc * tc * tc)
    d
  }
  #
  # Equation for the temperature-dependent bulk modulus of water at 1 atm, atm
  Ko <- function(tc){
    C1 <- 19652.17
    C2 <- 148.183
    C3 <- -2.29995
    C4 <- 0.01281
    C5 <- -4.91564 * 10^(-5)
    C6 <- 1.03553 * 10^(-7)
    # 
    k <- C1 +
      C2 * tc +
      C3 * tc * tc +
      C4 * tc * tc * tc +
      C5 * tc * tc * tc * tc +
      C6 * tc * tc * tc * tc * tc
    k
  }
  #
  # Equations for the temperature-dependend coefficients, A & B
  A <- function(tc){
    C1 <- 3.26138
    C2 <- 5.223 * 10^(-4)
    C3 <- 1.324 * 10^(-4)
    C4 <- -7.655 * 10^(-7)
    C5 <- 8.584 * 10^(-10)
    #
    a <- C1 +
      C2 * tc +
      C3 * tc * tc +
      C4 * tc * tc * tc +
      C5 * tc * tc * tc * tc 
    a
  }
  B <- function(tc){
    # T :: temperature, deg. C
    C1 <- 7.2061 * 10^(-5)
    C2 <- -5.8948 * 10^(-6)
    C3 <- 8.699 * 10^(-8)
    C4 <- -1.01 * 10^(-9)
    C5 <- 4.322 * 10^(-12)
    #
    b <- C1 +
      C2 * tc +
      C3 * tc * tc +
      C4 * tc * tc * tc +
      C5 * tc * tc * tc * tc 
    b
  }
  #
  # Convert pressure to bar (1 bar = 100000 Pa)
  Pbar <- (pa/100000.0)
  #
  rho <- (
    1000*do(tc)*(Ko(tc) + A(tc)*Pbar + B(tc)*Pbar^2)/(
      Ko(tc) + A(tc)*Pbar + B(tc)*Pbar^2 - Pbar)
  )
  rho
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
  # Equation for temperature-dependent specific heat capacity of water, J/kg/K
  Cp <- function(tc){
    1.0*10^3*(
      1.0045714270+
        2.050632750*10^(-3)*tc -
        1.631537093*10^(-4)*tc^(2) +
        6.212300300*10^(-6)*tc^(3) -
        8.830478888*10^(-8)*tc^(4) +
        5.071307038*10^(-10)*tc^(5)
    )
  }
  #
  Cp(tc)*kMa*pa/(kMv*enthalpy_vap(tc))
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
ka <- 1.4960e8        # length of earth's semi-major axis, km
kalb_sw <- 0.17       # shortwave albedo (Federer, 1968)
kalb_vis <- 0.03      # visible light albedo (Sellers, 1985)
kb <- 0.20            # constant for Rl (Linacre, 1968; Kramer, 1957)
kc <- 0.25            # constant for Rs (Linacre, 1968)
kCw <- 1.05           # supply constant, mm/hr (Federer, 1982)
kd <- 0.50            # constant for Rs (Linacre, 1968)
ke <- 0.01670         # eccentricity of earth's orbit, 2000CE (Berger 1978)
keps <- 23.44         # obliquity of earth's elliptic, 2000CE (Berger 1978)
kfFEC <- 2.04         # from-flux-to-energy, umol/J (Meek et al., 1984)
kG <- 9.80665         # gravitational acceleration, m/s^2
kGM <- 1.32712e11     # sun's standard gravity, km3/s2 (NASA's Sun Fact Sheet)
kGsc <- 1360.8        # solar constant, W/m^2 (Kopp & Lean, 2011)
kL <- 0.0065          # adiabatic lapse rate, K/m (Cavcar, 2000)
kMa <- 0.028963       # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kMv <- 0.01802        # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
komega <- (103+180)   # lon. of perihelion, degrees, 2000CE (Berger, 1978)
kSecInDay <- 86400    # number of seconds in a day
kPo <- 101325         # mean sea-level pressure, Pa (Cavcar, 2000)
kR <- 8.314           # universal gas constant, J/mol/K (Farquhar et al., 1980)
ksb <- 5.670373e-8    # Stefan-Boltzman constant, W/m^2/K^4
kTo <- 298.15         # base temperature (25 deg C), K (I.C. Prentice)
kWm <- 150            # soil moisture capacity, mm (Cramer-Prentice, 1988)
kw <- 0.26            # PET entrainment, (1+kw)*EET (Priestley-Taylor, 1972)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### STASH 2.0 ################################################################
# /////////////////////////////////////////////////////////////////////////////
all_years <- seq(from=2000, to=2010, by=1)
all_months <- seq(from=1, to=12, by=1)
my_lon <- -0.641
my_lat <- 51.4
my_elv <- 74

# Average Climate Data, Met Office (Ascot Racecourse)
# [http://www.metoffice.gov.uk/public/weather/climate/gcps2xe58]
# * total sunshine hours from Doorenbos & Pruitt (1977), Table 3
#   $fsun, unitless
#   $tair, deg C
#   $pre, mm
DATA <- matrix(nrow=12, ncol=3)
DATA <- as.data.frame(DATA)
names(DATA) <- c('fsun', 'tair', 'pre')
DATA$fsun <- c(0.21, 0.27, 0.30, 0.40, 0.39, 0.39, 0.40, 0.43, 0.36, 0.32, 0.23, 0.19)
DATA$tair <- c(4.80, 4.85, 7.10, 9.10, 12.4, 15.3, 17.6, 17.3, 14.6, 11.2, 7.55, 5.05)
DATA$pre  <- c(61.0, 41.2, 44.5, 48.0, 46.4, 44.6, 46.0, 52.3, 50.3, 71.8, 66.3, 62.9)

# Initial soil moisture:
W <- (0*seq(366))

# Initialize daily results:
daily_results <- matrix(data=rep(0,2928),nrow=366, ncol=8)
daily_results <- as.data.frame(daily_results)
names(daily_results) <- c('Ra_J.m2', 'Rn_J.m2', 'PPFD_mol.m2',
                          'Cond_mm', 'EET_mm', 'D_mm', 'S_mm', 'AET_mm')

# Initialize monthly results:
monthly_results <- matrix(data=rep(0,84),nrow=12, ncol=7)
monthly_results <- as.data.frame(monthly_results)
names(monthly_results) <- c("AET_mm", "EET_mm", "PET_mm", "CP_alpha", 
                            "CWD_mm", "PPFD_mol.m2", "RO_mm")

# yearly:
for (year in all_years){
  # Calculate days in the year
  DYR <- (julian_day(year+1, 1, 1) - julian_day(year, 1, 1))
  #
  # monthly:
  for (month in all_months){
    # Get DOM (days of current month):
    DOM <- (julian_day(year, month+1, 1) - julian_day(year, month, 1))
    #
    # daily:
    for (day in seq(from=1, to=DOM, by=1)){
      # Calculate day of year:
      DOY <- (julian_day(year, month, day) - julian_day(year, 1 , 1) + 1)
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 1. Evaporative supply (mm/hr)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Based on yesterday's available soil moisture
      # ref: Federer (1982); Eq. 4, Prentice et al. (1993)
      idx <- (DOY-1)
      if (idx < 1){
        idx <- DYR
      }
      Sw <- kCw*(W[idx]/kWm)
      #
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 2. Actual evapotranspiration
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Computes dailylight hours, Rs, Rl, Rn, daily_Rn, nighttime Rn (cond.?),
      # PET, and AET
      ET <- evap(
        user.lon=my_lon, 
        user.lat=my_lat, 
        user.day=DOY, 
        user.year=year, 
        user.fsun=DATA$fsun[month], 
        user.tair=DATA$tair[month], 
        user.elv=my_elv,
        user.supply=Sw,
        user.dr='loutre',            # 'loutre' or 'klein'
        user.delta='loutre',         # 'loutre', 'cooper' or 'circle'
        user.lambda='kepler')        # 'kepler', 'woolf' or 'berger'
      #
      # Save daily results:
      daily_results$Ra_J.m2[DOY] <- ET$ra_j.m2
      daily_results$Rn_J.m2[DOY] <- ET$rn_j.m2
      daily_results$PPFD_mol.m2[DOY] <- ET$ppfd_mol.m2
      daily_results$Cond_mm[DOY] <- ET$cond_mm
      daily_results$EET_mm[DOY] <- ET$eet_mm
      daily_results$D_mm[DOY] <- ET$di_mm
      daily_results$S_mm[DOY] <- ET$si_mm
      daily_results$AET_mm[DOY] <- ET$aet_mm
      #
      # Update monthly totals:
      monthly_results$AET_mm[month] <- (
        monthly_results$AET_mm[month] + ET$aet_mm
      )
      monthly_results$EET_mm[month] <- (
        monthly_results$EET_mm[month] + ET$eet_mm
      )
      monthly_results$PET_mm[month] <- (
        monthly_results$PET_mm[month] + ET$di_mm
      )
      monthly_results$PPFD_mol.m2[month] <- (
        monthly_results$PPFD_mol.m2[month] + ET$ppfd_mol.m2
      )
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 3. Update daily soil moisture (mm)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Today's soil moisture equals yesterday's soil moisture plus this 
      # month's fraction of precipitation plus condensation minus AET and
      # any excess is considered runoff
      # ref: Cramer-Prentice (1988)
      RO <- (W[idx] + (DATA$pre[month]/DOM) + ET$cond_mm - ET$aet_mm)
      #
      if (RO > kWm){
        # Set soil moisture to capacity
        W[DOY] <- kWm
        #
        # Add excess to monthly runoff
        monthly_results$RO_mm[month] <- (
          monthly_results$RO_mm[month] + RO - kWm 
        )
      } else if (RO < 0) {
        # Demand exceeds supply, soil is left dry:
        W[DOY] <- 0
      } else {
        # No runoff
        W[DOY] <- RO
      }
      #
    } # end daily
    monthly_results$CP_alpha[month] <- (
      monthly_results$AET_mm[month]/monthly_results$PET_mm[month]
    )
    monthly_results$CWD_mm[month] <- (
      monthly_results$PET_mm[month] - monthly_results$AET_mm[month]
    )
  } # end monthly
} # end yearly

#### View Results ####
##
## Plot monthly results
##
par(mar=c(4,4.5,1,1))
plot(monthly_results$PET_mm,type='l',lwd=2, col='black',
     ylim=c(0,max(monthly_results$PET_mm)), xlab=NA,ylab=NA,axes=F)
lines(monthly_results$EET_mm,lty=1,lwd=2,col='orange')
lines(monthly_results$AET_mm,lty=2,lwd=2,col='blue')
lines(monthly_results$CWD_mm,lty=1,lwd=2,col='red')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, 'Month', line=2)
mtext(side=2, as.expression(ET~(mm)), line=3)
legend('topleft',legend=c("PET","EET","AET","CWD"),
       col=c('black','orange','blue','red'),
       lty=c(1,1,2,1),lwd=c(2,2,2,2), inset=0.01,
       y.intersp=1.2,horiz=FALSE, bty='n')

##
## Plot daily ET results
##
par(mar=c(4,4.5,1,1))
plot(daily_results$D_mm, type='l', lwd=2, col='black', 
     ylim=c(0,15),
     xlab=NA, ylab=NA, axes=F)
lines(daily_results$S_mm, lty=1, lwd=2, col='red')
lines(daily_results$EET_mm, lty=1, lwd=2, col='orange')
lines(daily_results$AET_mm, lty=2, lwd=2, col='blue')
lines(daily_results$Cond_mm, lty=1, lwd=2, col='green')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, 'Day', line=2)
mtext(side=2, as.expression(ET~(mm)), line=3)
legend('top',legend=c("Demand", "Supply", "EET", "AET", "Cond"),
       col=c('black','red','orange','blue', 'green'),
       lty=c(1,1,1,2,1),lwd=c(2,2,2,2,2), inset=0.01,
       y.intersp=1.0,horiz=TRUE, bty='n')

##
## Plot single day ET results
##
par(mar=c(4,4.5,1,1))
plot(ET$time_hr, ET$pet_mm.hr, type='l', lwd=2, xlab=NA,ylab=NA,axes=F)
lines(ET$time_hr, ET$eet_mm.hr, lty=1, lwd=2, col='orange')
lines(ET$time_hr, ET$aet_mm.hr, lty=2, lwd=2, col='blue')
lines(ET$time_hr, (ET$pet_mm.hr - ET$aet_mm.hr), lwd=2, col='red')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, as.expression(Time~(hr)), line=2)
mtext(side=2, as.expression(ET~(mm%.%hr^-1)), line=3)
legend('topleft',legend=c("PET","EET","AET","CWD"),
       col=c('black','orange','blue','red'),
       lty=c(1,1,2,1),lwd=c(2,2,2,2), inset=0.01,
       y.intersp=1.2,horiz=FALSE, bty='n')

##
## Figure for AET
##
par(mar=c(4,4.5,1,1))
plot(ET$time_hr,ET$pet_mm.hr,type='l',lwd=3,lty=2,xlab=NA,ylab=NA,axes=F)
lines(ET$time_hr,ET$aet_mm.hr,lty=1,lwd=3,col='black')
abline(h=Sw, lty=3, lwd=3)
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, as.expression(Time~(hr)), line=2)
mtext(side=2, as.expression(ET~(mm%.%hr^-1)), line=2)
legend('topleft',legend=c("AET",
                          expression(italic(D[p])),
                          expression(italic(S[w]))),
       col=c('black','black','black'),
       lty=c(1,2,3),lwd=c(3,3,3), inset=0.02,
       y.intersp=1.5,horiz=FALSE, bty='n')

##
## Plot single day ET radiation results
##
par(mar=c(4,4.5,1,1))
plot(ET$time_hr, ET$rn_w.m2, type='l', lwd=2, xlab=NA, ylab=NA, 
     ylim=c(min(ET$rn_w.m2),max(ET$ra_w.m2)), axes=F)
lines(ET$time_hr, ET$ra_w.m2, col='red', lty=1, lwd=2)
lines(ET$time_hr, ET$rs_w.m2, col='blue', lty=2, lwd=2)
abline(a=ET$rl_w.m2, b=0, col='green', lty=4, lwd=2)
abline(h=0, lty=2)            # Zero line
abline(v=ET$hs_lct, lty=2)    # Sunset hour
abline(v=ET$hn_lct, lty=2)    # Rn cross-over hour
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, as.expression(Time~(hr)), line=2)
mtext(side=2, as.expression(Radiative~Flux~(W%.%m^-2)), line=3)
legend('topleft',legend=c(expression(italic(R[a])),
                          expression(italic(R[s])),
                          expression(italic(R[l])),
                          expression(italic(R[n]))),
       col=c('red','blue','green','black'),
       lty=c(1,2,4,1),lwd=c(2,2,2,2), inset=0.01,
       y.intersp=1.2,horiz=FALSE, bty='n')

##
## Daily Ra Figure
##
sol_noon <- 12-ET$eot_min/60+ET$lc_hr
par(mar=c(4,4.5,1,1))
plot(ET$time_hr, ET$ra_w.m2, type='l', lwd=2, xlab=NA, ylab=NA, 
     ylim=c(min(ET$rn_w.m2),max(ET$ra_w.m2)), axes=F)
abline(h=0, lty=2)            # Zero line
abline(v=ET$hs_lct, lty=2)    # Sunset hour
abline(v=sol_noon, lty=2)     # Solar noon
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, as.expression(Time~(hr)), line=2)
mtext(side=2, as.expression(Radiative~Flux~(W%.%m^-2)), line=3)

##
## Daily Rn Figure
##
sol_noon <- 12-ET$eot_min/60+ET$lc_hr
x <- c(ET$time_hr,24)
y <- c(ET$rn_w.m2,ET$rn_w.m2[length(ET$rn_w.m2)])
plot(x, y, type='l', lwd=2, xlab=NA, ylab=NA, axes=F)
abline(h=0, lty=2)            # Zero line
abline(v=ET$hn_lct, lty=2)    # Rn cross-over hour
abline(v=ET$hs_lct, lty=2)    # Sunset hour
#abline(v=sol_noon, lty=2)     # Solar noon
abline(v=24, lty=2)           # Pi
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, as.expression(Time~(hr)), line=2)
mtext(side=2, as.expression(Radiative~Flux~(W%.%m^-2)), line=3)
