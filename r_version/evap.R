# R version 3.2.2 (2015-08-14) -- "Fire Safety"
#
# evap.R
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
# This script contains functions to calculate daily radiation, condensation,
# and evapotranspiration, i.e.:
#   berger_tls(double n, double N)
#   density_h2o(double tc, double pa)
#   dcos(double d)
#   dsin(double d)
#   elv2pres(double z)
#   enthalpy_vap(double tc)
#   evap(double lat, double n, double elv=0, double y=0, double sf=1,
#        double tc=23.0, double sw=1.0)
#   julian_day(double y, double m, double i)
#   psychro(double tc, double pa)
#   sat_slope(double tc)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - fixed Cooper's and Spencer's declination angle equations [14.11.25]
# - replaced simplified_kepler with full_kepler method [14.11.25]
# - added berger_tls function [15.01.13]
# - updated evap function (similar to stash.py EVAP class) [15.01.13]
#
#### IMPORT SOURCES ##########################################################
source("const.R")


#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     berger_tls
# Inputs:   - double, day of year (n)
#           - double, days in year (N)
# Returns:  numeric list, true anomaly and true longitude
# Features: Returns true anomaly and true longitude for a given day.
# Depends:  - ke ............. eccentricity of earth's orbit, unitless
#           - komega ......... longitude of perihelion
#  Ref:     Berger, A. L. (1978), Long term variations of daily insolation
#             and quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
# ************************************************************************
berger_tls <- function(n, N) {
    # Variable substitutes:
    xee <- ke^2
    xec <- ke^3
    xse <- sqrt(1 - ke^2)

    # Mean longitude for vernal equinox:
    xlam <- (ke/2.0 + xec/8.0)*(1 + xse)*sin(komega*pir) -
        xee/4.0*(0.5 + xse)*sin(2.0*komega*pir) +
        xec/8.0*(1.0/3.0 + xse)*sin(3.0*komega*pir)
    xlam <- 2.0*xlam/pir

    # Mean longitude for day of year:
    dlamm <- xlam + (n - 80.0)*(360.0/N)

    # Mean anomaly:
    anm <- dlamm - komega
    ranm <- anm*pir

    # True anomaly (uncorrected):
    ranv <- ranm + (2.0*ke - xec/4.0)*sin(ranm) +
        5.0/4.0*xee*sin(2.0*ranm) +
        13.0/12.0*xec*sin(3.0*ranm)
    anv <- ranv/pir

    # True longitude:
    my_tls <- anv + komega
    if (my_tls < 0){
        my_tls <- my_tls + 360
    } else if (my_tls > 360) {
        my_tls <- my_tls - 360
    }

    # True anomaly:
    my_nu <- my_tls - komega
    if (my_nu < 0){
        my_nu <- my_nu + 360
    }
    return (c(my_nu, my_tls))
}


# ************************************************************************
# Name:     density_h2o
# Inputs:   - double (tc), air temperature, degrees C
#           - double (pa), atm pressure, Pa
# Returns:  double, kg/m^3
# Features: This function calculates the temperature and pressure
#           dependent density of pure water
# * Ref:    Chen, C.T., R.A. Fine, and F.J. Millero (1977), The equation
#             of state of pure water determined from sound speeds, The
#             Journal of Chemical Physics 66, 2142;
#             doi:10.1063/1.434179
# ************************************************************************
density_h2o <- function(tc, pa) {
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

    # Calculate the bulk modulus of water at 1 atm, atm
    ko <- 19652.17 +
        148.1830*tc +
        -2.29995*tc*tc +
        0.01281*tc*tc*tc +
        -(4.91564e-5)*tc*tc*tc*tc +
        (1.035530e-7)*tc*tc*tc*tc*tc

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

    # Convert pressure to bar (1 bar = 100000 Pa)
    pbar <- (1e-5)*pa

    pw <- (1e3)*po*(ko + ca*pbar + cb*pbar^2)/(ko + ca*pbar + cb*pbar^2 - pbar)
    return(pw)
}


# ************************************************************************
# Name:     dcos
# Inputs:   double (d), angle in degrees
# Returns:  double, cosine of angle
# Features: This function calculates the cosine of an angle (d) given
#           in degrees.
# Ref:      This script is based on the Javascript function written by
#           C Johnson, Theoretical Physicist, Univ of Chicago
#           - 'Equation of Time' URL: http://mb-soft.com/public3/equatime.html
#           - Javascript URL: http://mb-soft.com/believe/txx/astro22.js
# ************************************************************************
dcos <- function(d) {
  cos(d*pi/180)
}


# ************************************************************************
# Name:     dsin
# Inputs:   double (d), angle in degrees
# Returns:  double, sine of angle
# Features: This function calculates the sine of an angle (d) given
#           in degrees.
# ************************************************************************
dsin <- function(d) {
  sin(d*pi/180)
}


# ************************************************************************
# Name:     elv2pres
# Inputs:   double (z), meters
# Returns:  double, Pa
# Features: Calculates atmospheric pressure for a given elevation
# Depends:  - kPo ............ base pressure, Pa
#           - kTo ............ base temperature, C
#           - kL ............. temperature lapse rate, K/m
#           - kR ............. universal gas constant, J/mol/K
#           - kMa ............ molecular weight of dry air, kg/mol
#           - kG ............. gravity, m/s^2
# Ref:      Allen et al. (1998)
# ************************************************************************
elv2pres <- function(z) {
    kPo*(1 - kL*z/kTo)^(kG*kMa/(kR*kL))
}


# ************************************************************************
# Name:     enthalpy_vap
# Inputs:   double (tc), degrees C
# Returns:  double, J/kg
# Features: This function calculates the temperature-dependent enthalpy
#           of vaporization (latent heat of vaporization)
# Ref:      Eq. 8, Henderson-Sellers (1984), A new formula for latent heat
#             of vaporization of water as a function of temperature, Quarterly
#             Journal of the Royal Meteorological Society, vol. 110, pp. 1186--
#             1190.
# ************************************************************************
enthalpy_vap <- function(tc) {
    1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))^2
}


# ************************************************************************
# Name:     evap
# Inputs:   - double, latitude, degrees (lat)
#           - double, day of year (n)
#           - double, elevation (elv)  *optional
#           - double, year (y)         *optional
#           - double, fraction of sunshine hours (sf)        *optional
#           - double, mean daily air temperature, deg C (tc) *optional
#           - double, evaporative supply rate, mm/hr (sw)    *optional
# Returns:  list object (et.srad)
#             $nu_deg ............ true anomaly, degrees
#             $lambda_deg ........ true longitude, degrees
#             $dr ................ distance factor, unitless
#             $delta_deg ......... declination angle, degrees
#             $hs_deg ............ sunset angle, degrees
#             $ra_j.m2 ........... daily extraterrestrial radiation, J/m^2
#             $tau ............... atmospheric transmittivity, unitless
#             $ppfd_mol.m2 ....... daily PPFD, mol/m^2
#             $hn_deg ............ net radiation hour angle, degrees
#             $rn_j.m2 ........... daily net radiation, J/m^2
#             $rnn_j.m2 .......... daily nighttime net radiation, J/m^2
#             $econ_m3.j ......... water to energy conversion, m^3/J
#             $cond_mm ........... daily condensation, mm
#             $eet_mm ............ daily equilibrium ET, mm
#             $pet_mm ............ daily potential ET, mm
#             $hi_deg ............ intersection hour angle, degrees
#             $aet_mm ............ daily actual ET, mm
# Features: This function calculates daily radiation, condenstaion, and
#           evaporation fluxes.
# Depends:  - kalb_sw ........ shortwave albedo
#           - kalb_vis ....... visible light albedo
#           - kb ............. empirical constant for longwave rad
#           - kc ............. empirical constant for shortwave rad
#           - kd ............. empirical constant for shortwave rad
#           - ke ............. eccentricity
#           - keps ........... obliquity
#           - kfFEC .......... from-flux-to-energy conversion, umol/J
#           - kGsc ........... solar constant
#           - kw ............. entrainment factor for PET
#           - berger_tls() ... calc true anomaly and longitude
#           - dcos() ......... cos(x*pi/180), where x is in degrees
#           - dsin() ......... sin(x*pi/180), where x is in degrees
#           - density_h2o() .. density of water
#           - elv2pres() ..... elevation dependent atm. pressure
#           - enthalpy_vap() . latent heat of vaporization
#           - julian_day() ... date to julian day
#           - psychro() ...... psychrometric constant
#           - sat_slope() .... slope of sat. pressure temp curve
# ************************************************************************
evap <- function(lat, n, elv=0, y=0, sf=1, tc=23.0, sw=1.0) {
    # ~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if (lat > 90 || lat < -90) {
        stop("Warning: Latitude outside range of validity (-90 to 90)!")
    }
    if (n < 1 || n > 366) {
        stop("Warning: Day outside range of validity (1 to 366)!")
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    et.srad <- list()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 01. Calculate the number of days in yeark (kN), days
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (y == 0) {
        kN <- 365
    } else {
        kN <- (julian_day(y + 1, 1, 1) - julian_day(y, 1, 1))
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 02. Calculate heliocentric longitudes (nu and lambda), degrees
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    my_helio <- berger_tls(n, kN)
    my_nu <- my_helio[1]
    my_lambda <- my_helio[2]

    et.srad$nu_deg <- my_nu
    et.srad$lambda_deg <- my_lambda

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03. Calculate distance factor (dr), unitless
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Berger et al. (1993)
    kee <- ke^2
    my_rho <- (1 - kee)/(1 + ke*dcos(my_nu))
    dr <- (1/my_rho)^2

    et.srad$dr <- dr

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 04. Calculate the declination angle (delta), degrees
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Woolf (1968)
    delta <- asin(dsin(my_lambda)*dsin(keps))
    delta <- delta/pir

    et.srad$delta_deg <- delta

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 05. Calculate variable substitutes (u and v), unitless
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ru <- dsin(delta)*dsin(lat)
    rv <- dcos(delta)*dcos(lat)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 06. Calculate the sunset hour angle (hs), degrees
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Note: u/v equals tan(delta) * tan(lat)
    if (ru/rv >= 1.0) {
        hs <- 180  # Polar day (no sunset)
    } else if (ru/rv <= -1.0) {
        hs <- 0 # Polar night (no sunrise)
    } else {
        hs <- acos(-1.0*ru/rv)
        hs <- hs / pir
    }
    et.srad$hs_deg <- hs

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 07. Calculate daily extraterrestrial radiation (ra_d), J/m^2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq. 1.10.3, Duffy & Beckman (1993)
    ra_d <- (86400/pi)*kGsc*dr*(ru*pir*hs + rv*dsin(hs))
    et.srad$ra_j.m2 <- ra_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 08. Calculate transmittivity (tau), unitless
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref:  Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
    tau_o <- (kc + kd*sf)
    tau <- tau_o*(1 + (2.67e-5)*elv)

    et.srad$tau <- tau

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 09. Calculate daily PPFD (ppfd_d), mol/m^2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ppfd_d <- (1e-6)*kfFEC*(1 - kalb_vis)*tau*ra_d
    et.srad$ppfd_mol.m2 <- ppfd_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 10. Estimate net longwave radiation (rnl), W/m^2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnl <- (kb + (1 - kb)*sf)*(kA - tc)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 11. Calculate variable substitue (rw), W/m^2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw <- (1 - kalb_sw)*tau*kGsc*dr

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 12. Calculate net radiation cross-over angle (hn), degrees
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((rnl - rw*ru)/(rw*rv) >= 1.0) {
        hn <- 0  # Net radiation is negative all day
    } else if ((rnl - rw*ru)/(rw*rv) <= -1.0) {
        hn <- 180 # Net radiation is positive all day
    } else {
        hn <- acos((rnl - rw*ru)/(rw*rv))
        hn <- hn/pir
    }
    et.srad$hn_deg <- hn

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 13. Calculate daytime net radiation (rn_d), J/m^2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rn_d <- (86400/pi)*(hn*pir*(rw*ru - rnl) + rw*rv*dsin(hn))
    et.srad$rn_j.m2 <- rn_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 14. Calculate nighttime net radiation (rnn_d), J/m^2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnn_d <- (86400/pi)*(
        rw*ru*(hs - hn)*pir + rw*rv*(dsin(hs) - dsin(hn)) +
            rnl*(pi - 2*hs*pir + hn*pir)
    )
    et.srad$rnn_j.m2 <- rnn_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 15. Calculate water-to-energy conversion (econ), m^3/J
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Slope of saturation vap press temp curve, Pa/K
    s <- sat_slope(tc)

    # Enthalpy of vaporization, J/kg
    lv <- enthalpy_vap(tc)

    # Density of water, kg/m^3
    pw <- density_h2o(tc, elv2pres(elv))

    # Psychrometric constant, Pa/K
    g <- psychro(tc, elv2pres(elv))

    econ <- s/(lv*pw*(s + g))
    et.srad$econ_m3.j <- econ

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 16. Calculate daily condensation (cn), mm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cn <- (1e3)*econ*abs(rnn_d)
    et.srad$cond_mm <- cn

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 17. Estimate daily equilibrium evapotranspiration (eet_d), mm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eet_d <- (1e3)*econ*rn_d
    et.srad$eet_mm <- eet_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 18. Estimate daily potential evapotranspiration (pet_d), mm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pet_d <- (1 + kw)*eet_d
    et.srad$pet_mm <- pet_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx <- (3.6e6)*(1 + kw)*econ

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 20. Calculate the intersection hour angle (hi), degrees
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cos_hi <- sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
    if (cos_hi >= 1.0) {
        hi <- 0.0       # supply exceeds demand
    } else if (cos_hi <= -1.0) {
        hi <- 180.0     # supply limits demand everywhere
    } else {
        hi <- acos(cos_hi)
        hi <- hi/pir
    }
    et.srad$hi_deg <- hi

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 21. Estimate daily actual evapotranspiration (aet_d), mm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aet_d <- (24/pi)*(
        sw*hi*pir +
            rx*rw*rv*(dsin(hn) - dsin(hi)) +
                (rx*rw*ru - rx*rnl)*(hn - hi)*pir
    )
    et.srad$aet_mm <- aet_d

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    et.srad
}


# ************************************************************************
# Name:     julian_day
# Inputs:   - double, year (y)
#           - double, month (m)
#           - double, day of month (i)
# Returns:  double, Julian day
# Features: This function converts a date in the Gregorian calendar
#           to a Julian day number (i.e., a method of consecutative
#           numbering of days---does not have anything to do with
#           the Julian calendar!)
# Ref:      Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", Astronomical
#             Algorithms
# ************************************************************************
julian_day <- function(y, m, i) {
    if (m <= 2) {
        y <- y - 1
        m <- m + 12
    }
    a <- floor(y/100)
    b <- 2 - a + floor(a/4)

    jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
    return(jde)
}


# ************************************************************************
# Name:     psychro
# Inputs:   - double (tc), air temperature, degrees C
#           - double (pa), atm pressure, Pa
# Returns:  double, Pa/K
# Features: This function calculates the temperature and pressure
#           dependent psychrometric constant
# Depends:  enthalpy_vap
# Ref:      - Temperature dependency of psychrometric constant (Eq. 8):
#               Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998),
#                 'Meteorological data,' Crop evapotranspiration - Guidelines
#                 for computing crop water requirements - FAO Irrigation and
#                 drainage paper 56, Food and Agriculture Organization of the
#                 United Nations, Available:
#                 http://www.fao.org/docrep/x0490e/x0490e07.htm
#           - Pressure dependency of psychrometric constant (i.e., Cp):
#               Tsilingris (2008), Thermophysical and transport properties of
#                 humid air at temperature range between 0 and 100 °C, Energy
#                 Conversion and Management, vol. 49, pp. 1098--1110.
# ************************************************************************
psychro <- function(tc, pa) {
    # Calculate the specific heat capacity of water, J/kg/K
    cp <- 1.0045714270 +
        (2.050632750e-3)*tc -
        (1.631537093e-4)*tc*tc +
        (6.212300300e-6)*tc*tc*tc -
        (8.830478888e-8)*tc*tc*tc*tc +
        (5.071307038e-10)*tc*tc*tc*tc*tc
    cp <- (1e3)*cp

    # Calculate latent heat of vaporization, J/kg
    lv <- enthalpy_vap(tc)

    # Calculate psychrometric constant, Pa/K
    return(cp*kMa*pa/(kMv*lv))
}


# ************************************************************************
# Name:     sat_slope
# Inputs:   double (tc), degrees C
# Returns:  double, Pa/K
# Features: This function calculates the temperature-dependent slope of
#           the saturation pressure temperature curve using the
#           methodology presented in the eMast energy.cpp script
# Ref:      - Eq. 6, Prentice et al. (1993);
#           - Eq. 13, Allen et al. (1998)
# ************************************************************************
sat_slope <- function(tc) {
    (17.269)*(237.3)*(610.78)*exp(17.269*tc/(237.3 + tc))/(237.3 + tc)^2
}
