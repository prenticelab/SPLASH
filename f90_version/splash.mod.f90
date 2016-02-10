module _splash
  !////////////////////////////////////////////////////////////////
  ! SPLASH 
  ! This module contains all functions and parameter values to run
  ! the SPLASH model, given inputs 
  ! - temperature (deg C, daily values)
  ! - precipitation (mm/day, daily values)
  ! - sunshine fraction (unitless, daily values)
  ! - latitude (deg N)
  ! - elevation (m.a.s.l.)
  ! 
  ! author: T.W. Davis, Fortran version by B. Stocker
  !
  ! last updated: 2016-02-05
  !
  ! citation:
  ! T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
  ! Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
  ! led algorithms for simulating habitats (SPLASH): Robust indices of radiation,
  ! evapotranspiration and plant-available moisture, Geoscientific Model
  ! Development, 2016 (in progress)
  !----------------------------------------------------------------   
  implicit none

  ! soil water balance variables as derived type
  type waterbaltype
    real :: sm     ! soil moisture (=water content, mm)
    real :: ro     ! daily runoff (mm)
  end type

  ! function return value of evap() as derived type
  type outtype_evap
    real :: ra         ! daily TOA solar irradiation (J/m2)
    real :: rn         ! daily net radiation (J/m2)
    real :: ppfd       ! daily PPFD (mol/m^2)
    real :: eet        ! daily out_evap%EET (mm)
    real :: pet        ! daily PET (mm)
    real :: aet        ! daily AET (mm)
    real :: cn         ! daily condensation (mm)
  end type

  ! function return value of berger_tls() as derived type
  type outtype_berger
    real :: nu
    real :: lambda
  end type

  ! SPLASH state variables, daily updated
  type( waterbaltype ) :: waterbal
  type( outtype_evap ) :: out_evap

  ! verbose prints various state variables to screen
  logical :: verbose

  ! Chose whether to use monthly or daily input data
  logical :: use_daily_input

  !----------------------------------------------------------------   
  ! model parameters
  !----------------------------------------------------------------   
  real, parameter :: kA = 107             ! constant for Rnl (Monteith & Unsworth, 1990)
  real, parameter :: kalb_sw = 0.17       ! shortwave albedo (Federer, 1968)
  real, parameter :: kalb_vis = 0.03      ! visible light albedo (Sellers, 1985)
  real, parameter :: kb = 0.20            ! constant for Rnl (Linacre, 1968)
  real, parameter :: kc = 0.25            ! cloudy transmittivity (Linacre, 1968)
  real, parameter :: kCw = 1.05           ! supply constant, mm/hr (Federer, 1982)
  real, parameter :: kd = 0.50            ! angular coefficient of transmittivity (Linacre, 1968)
  real, parameter :: ke = 0.0167          ! eccentricity for 2000 CE (Berger, 1978)
  real, parameter :: keps = 23.44         ! obliquity for 2000 CE, degrees (Berger, 1978)
  real, parameter :: kfFEC = 2.04         ! from flux to energy conversion, umol/J (Meek et al., 1984)
  real, parameter :: kG = 9.80665         ! gravitational acceleration, m/s^2 (Allen, 1973)
  real, parameter :: kGsc = 1360.8        ! solar constant, W/m^2 (Kopp & Lean, 2011)
  real, parameter :: kL = 0.0065          ! temperature lapse rate, K/m (Cavcar, 2000)
  real, parameter :: kMa = 0.028963       ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  real, parameter :: kMv = 0.01802        ! molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
  real, parameter :: kPo = 101325         ! standard atmosphere, Pa (Allen, 1973)
  real, parameter :: kR = 8.31447         ! universal gas constant, J/mol/K (Moldover et al., 1988)
  real, parameter :: kTo = 288.15         ! base temperature, K (Berberan-Santos et al., 1997)
  real, parameter :: kWm = 150            ! soil moisture capacity, mm (Cramer & Prentice, 1988)
  real, parameter :: kw = 0.26            ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
  real, parameter :: komega = 283.0       ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)
  real, parameter :: pi = 3.14159

  integer, parameter :: nmonth = 12       ! number of months in year

  !----------------------------------------------------------------     
  ! output variables
  !----------------------------------------------------------------   
  ! daily totals
  real, dimension(366) :: outdra            ! daily solar irradiation, J/m2
  real, dimension(366) :: outdrn            ! daily net radiation, J/m2
  real, dimension(366) :: outdppfd          ! daily PPFD, mol/m2
  real, dimension(366) :: outdcn            ! daily condensation water, mm
  real, dimension(366) :: outdsm            ! daily soil moisture, mm
  real, dimension(366) :: outdro            ! daily runoff, mm
  real, dimension(366) :: outdeet           ! daily equilibrium ET, mm
  real, dimension(366) :: outdpet           ! daily potential ET, mm
  real, dimension(366) :: outdaet           ! daily actual ET, mm

  ! monthly totals
  real, dimension(12) :: outmeet
  real, dimension(12) :: outmpet
  real, dimension(12) :: outmaet
  real, dimension(12) :: outmcpa
  real, dimension(12) :: outmcwd
  real, dimension(12) :: outmppfd


contains 


  subroutine spin_up_sm( yr, lat, elv, ppt, tc, sf, inlen )
    !----------------------------------------------------------------   
    ! Spins up the daily soil moisture
    !----------------------------------------------------------------   
    ! arguments
    integer, intent(in)                :: yr    ! year AD  
    real, intent(in)                   :: lat   ! latitude (degrees)
    real, intent(in)                   :: elv   ! altitude (m)
    real, dimension(inlen), intent(in) :: ppt   ! monthly precip (mm) 
    real, dimension(inlen), intent(in) :: tc    ! mean monthly temperature (deg C)
    real, dimension(inlen), intent(in) :: sf    ! mean monthly sunshine fraction (unitless)
    integer, intent(in)                :: inlen ! =12 if monthly input data, =ndayyear if daily input data

    ! local variables
    integer :: spin_count                            ! counter variable
    real :: start_sm                                 ! soil moisture in first day of year
    real :: end_sm                                   ! soil moisture in last day of year
    real :: diff_sm                                  ! difference in soil moisture between first and last day of year

    ! control variables
    spin_count = 1
    diff_sm = 9999.0

    do while (diff_sm>1.0e-4)

      ! Run one year
      print*,'********************************************'
      print*,'spinup year', spin_count
      print*,'--------------------------------------------'
      start_sm = waterbal%sm
      call run_one_year( lat, elv, yr, ppt, tc, sf, inlen )
      end_sm   = waterbal%sm

      ! Check to see if 1 Jan soil moisture matches 31 Dec:
      diff_sm  = abs(end_sm - start_sm)

      print*,'soil moisture in equil. when diff < 0.0001'
      print*,'diff ',diff_sm

      ! Increase counter
      spin_count = spin_count + 1

    enddo

  end subroutine spin_up_sm


  subroutine run_one_year( lat, elv, yr, ppt, tc, sf, inlen )
    !----------------------------------------------------------------   
    ! Calculates daily and monthly quantities for one year
    !----------------------------------------------------------------  
    ! arguments
    real, intent(in)                   :: lat   ! latitude (degrees)
    real, intent(in)                   :: elv   ! altitude (m)
    integer, intent(in)                :: yr    ! year AD
    real, dimension(inlen), intent(in) :: ppt   ! monthly precip (mm) 
    real, dimension(inlen), intent(in) :: tc    ! mean monthly temperature (deg C)
    real, dimension(inlen), intent(in) :: sf    ! mean monthly sunshine fraction (unitless)

    integer, intent(in) :: inlen                         ! =12 if monthly input data, =ndayyear if daily input data

    ! local variables
    real :: use_tc                       ! mean monthly or daily temperature (deg C)
    real :: use_sf                       ! mean monthly or daily sunshine fraction (unitless)
    real :: use_pr                       ! monthly or daily precipitation (mm)

    real :: sw                           ! evaporative supply rate (mm/h)

    integer :: doy                       ! day of year
    integer :: moy                       ! month of year
    integer :: ndaymonth                 ! number of days in month (formerly 'nm')
    integer :: idx                       ! day of year corresponding to yesterday
    integer :: dm                        ! day of month

    ! Reset monthly totals (output variables)
    call initmonthly

    ! Iterate through months
    do moy=1,nmonth

      ndaymonth = get_julian_day(yr, moy+1, 1) - get_julian_day(yr, moy, 1)

      ! Iterate through days in this month
      do dm=1,ndaymonth

        ! Calculate the day of the year
        ! xxx is this '+1' really necessary here?
        doy = get_julian_day(yr, moy, dm) - get_julian_day(yr, 1, 1) + 1

        ! Select daily or monthly input data
        if (use_daily_input) then
          use_tc   = tc(doy)
          use_sf   = sf(doy)
          use_pr   = ppt(doy)
        else
          use_tc   = tc(moy)
          use_sf   = sf(moy)
          use_pr   = ppt(moy)/ndaymonth
        end if

        call run_one_day( lat, elv, doy, yr, use_pr, use_tc, use_sf )

        ! Collect daily output variables
        call getout_daily( doy, moy ) 

      enddo

      call getout_monthly( moy )

    enddo

  end subroutine run_one_year


  subroutine run_one_day( lat, elv, doy, yr, pr, tc, sf )
    !----------------------------------------------------------------   
    ! Calculates daily and monthly quantities for one year
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in)    :: lat   ! latitude (degrees)
    real, intent(in)    :: elv   ! altitude (m)
    integer, intent(in) :: doy   ! altitude (m)
    integer, intent(in) :: yr    ! year AD
    real, intent(in)    :: pr    ! monthly or daily precip (mm) 
    real, intent(in)    :: tc    ! mean monthly or daily temperature (deg C)
    real, intent(in)    :: sf    ! mean monthly or daily sunshine fraction (unitless)

    ! local variables
    real  :: sw            ! evaporative supply rate (mm/h)

    ! Calculate evaporative supply rate, mm/h
    sw = kCw * waterbal%sm / kWm

    ! Calculate radiation and evaporation quantities
    out_evap = evap( lat, doy, elv, yr, sf, tc, sw )

    ! Update soil moisture
    waterbal%sm = waterbal%sm + pr + out_evap%cn - out_evap%aet

    if (waterbal%sm>kWm) then
      ! Bucket is full 
      ! * set soil moisture to capacity
      ! * add remaining water to monthly runoff total
      waterbal%ro = waterbal%sm - kWm  ! xxx add in sofun
      waterbal%sm = kWm

    elseif (waterbal%sm<0) then
      ! Bucket is empty
      ! * set soil moisture to zero
      out_evap%aet = out_evap%aet + waterbal%sm
      waterbal%ro = 0.0
      waterbal%sm = 0.0 ! xxx check in sofun

    else
      ! No runoff occurrs
      waterbal%ro = 0.0

    endif

  end subroutine run_one_day


  function evap( lat, doy, elv, yr, sf, tc, sw ) result( out_evap )
    !----------------------------------------------------------------   
    ! This subroutine calculates daily radiation and evapotranspiration
    ! quantities
    ! - daily PPFD (ppfd), mol/m^2
    ! - daily EET (eet), mm
    ! - daily PET (pet), mm
    ! - daily AET (aet), mm
    ! - daily condensation (wc), mm
    !-------------------------------------------------------------  
    ! arguments
    real,    intent(in) :: lat           ! latitude, degrees
    integer, intent(in) :: doy           ! day of the year (formerly 'n')
    real,    intent(in) :: elv ! elevation, metres
    integer, intent(in) :: yr  ! year
    real,    intent(in) :: sf  ! fraction of sunshine hours 
    real,    intent(in) :: tc  ! mean daily air temperature, C
    real,    intent(in) :: sw  ! evaporative supply rate, mm/hr

    ! function return variable
    type( outtype_evap )  :: out_evap
    type( outtype_berger) :: out_berger

    integer :: ndayyear

    real :: my_rho
    real :: dr                 ! distance factor
    real :: delta              ! declination angle 
    real :: ru                 ! variable substitute for u
    real :: rv                 ! variable substitute for v
    real :: hs                 ! sunset hour angle
    real :: tau_o, tau         ! transmittivity (unitless)
    real :: rnl                ! net longwave radiation (W/m^2)
    real :: rnn_d              ! nighttime net radiation (J/m^2)
    real :: rw                 ! variable substitute (W/m^2)
    real :: hn                 ! net radiation cross-over hour angle
    real :: s                  ! slope of saturation vap press temp curve, Pa/K
    real :: pw                 ! density of water, kg/m^3
    real :: lv                 ! enthalpy of vaporization, J/kg
    real :: g                  ! psychrometric constant, Pa/K
    real :: econ               ! Eq. 58, Documentation
    real :: rx                 ! variable substitute (mm/hr)/(W/m^2)
    real :: hi, cos_hi         ! intersection hour angle, degrees

    ! Error handle and assign required public variables:
    if (lat>90.0 .or. lat<(-90.0)) then
      print*,"Latitude outside range of validity (-90 to 90)!"
      stop
    endif

    if (doy<1 .or. doy>366) then
      print*,"Day of year outside range of validity (1 to 366)!"
      print*,"doy:",doy
      stop
    endif

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1. Calculate number of days in year (kN), days (take from argument)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ndayyear = get_julian_day(yr+1, 1, 1) - get_julian_day(yr, 1, 1)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger (1978)
    out_berger = berger_tls( doy, ndayyear )
    
    ! consistency check
    if (verbose) print*,'out_berger', out_berger

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3. Calculate distance factor (dr), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger et al. (1993)
    my_rho = (1.0 - ke**2)/(1.0 + ke * dgcos( out_berger%nu ))

    ! consistency check
    if (verbose) print*,'my_rho', my_rho

    dr = (1.0/my_rho)**2

    ! consistency check
    if (verbose) print*,'dr', dr

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 4. Calculate declination angle (delta), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Woolf (1968)
    delta = dasin( dgsin( out_berger%lambda ) * dgsin( keps ) )   ! xxx arcsin is asin in Fortran?
    delta = degrees( delta )

    ! consistency check
    if (verbose) print*,'delta', delta

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 5. Calculate variable substitutes (u and v), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ru = dgsin(delta) * dgsin(lat)
    rv = dgcos(delta) * dgcos(lat)

    ! consistency check
    if (verbose) print*,'ru, rv ',ru,rv

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 6. Calculate the sunset hour angle (hs), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Note: u/v == tan(delta)*tan(lat)
    ! Eq. 3.22, Stine & Geyer (2001)
    if ((ru/rv) >= 1.0) then
      ! Polar day (no sunset)
      hs = 180.0 
    elseif ((ru/rv) <= -1.0) then
      ! Polar night (no sunrise)
      hs = 0.0
    else
      hs = degrees(acos(-1.0*ru/rv))
    endif

    ! consistency check
    if (verbose) print*,'hs',hs

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 7. Calculate daily extraterrestrial solar radiation / irradiation (out_evap%ra), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 1.10.3, Duffy & Beckman (1993)
    out_evap%ra = (86400.0/pi)*kGsc*dr*(ru*radians(hs) + rv * dgsin(hs))

    ! consistency check
    if (verbose) print*,'daily extraterrestrial solar radiation / irradiation (out_evap%ra), J/m^2', out_evap%ra

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 8. Calculate transmittivity (tau), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Linacre (1968)
    tau_o = (kc + kd*sf)

    ! Eq. 2, Allen (1996)
    tau = tau_o*(1.0 + (2.67e-5)*elv)

    ! consistency check
    if (verbose) print*,'tau_o ',tau_o
    if (verbose) print*,'tau   ',tau

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 9. Calculate daily PPFD (ppfd_d), mol/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 57, Documentation
    out_evap%ppfd = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*out_evap%ra

    ! consistency check
    if (verbose) print*,'ppfd_d',out_evap%ppfd

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 10. Estimate net longwave radiation (rnl), W/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
    rnl = (kb + (1.0 - kb)*sf)*(kA - tc)

    ! consistency check
    if (verbose) print*,'net longwave radiation, rnl',rnl

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 11. Calculate variable substitute (rw), W/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw = (1.0-kalb_sw)*tau*kGsc*dr

    ! consistency check
    if (verbose) print*,'rw ', rw

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 12. Calculate net radiation cross-over hour angle (hn), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((rnl - rw*ru)/(rw*rv) >= 1.0) then
      ! Net radiation negative all day
      hn = 0.0
    else if ((rnl - rw*ru)/(rw*rv) <= -1.0) then
      ! Net radiation positive all day
      hn = 180.0
    else
      hn = degrees( acos((rnl - rw*ru)/(rw*rv)) )
    endif

    ! consistency check
    if (verbose) print*,'cross-over hour angle, hn ',hn

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 13. Calculate daytime net radiation (out_evap%rn), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 53, Documentation
    out_evap%rn = (86400.0/pi) * (hn*(pi/180.0)*(rw*ru - rnl) + rw*rv*dgsin(hn))

    ! consistency check
    if (verbose) print*,'daytime net radiation, out_evap%rn', out_evap%rn

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 14. Calculate nighttime net radiation (rnn_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 56, Documentation
    rnn_d = (86400.0/pi)*(radians(rw*ru*(hs-hn)) + rw*rv*(dgsin(hs)-dgsin(hn)) + rnl*(pi - 2.0*radians(hs) + radians(hn)))

    ! consistency check
    if (verbose) print*,'nighttime net radiation, rnn_d ',rnn_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 15. Calculate water-to-energy conversion (econ), m^3/J
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Slope of saturation vap press temp curve, Pa/K
    s = get_sat_slope( tc )

    ! consistency check
    if (verbose) print*,'Slope of saturation vap press temp curve, Pa/K, hn ', s

    ! Enthalpy of vaporization, J/kg
    lv = get_enthalpy_vap( tc )

    ! consistency check
    if (verbose) print*,'Enthalpy of vaporization, J/kg, lv ', lv

    ! Density of water, kg/m^3
    pw = get_density_h2o( tc, elv2pres( elv ) )

    ! consistency check
    if (verbose) print*,'Density of water, kg/m^3, pw ', pw

    ! Psychrometric constant, Pa/K
    g = get_psychro( tc, elv2pres( elv ) )

    ! consistency check
    if (verbose) print*,'Psychrometric constant, Pa/K ', g

    ! Eq. 58, Documentation
    econ = s/(lv*pw*(s + g))

    ! consistency check
    if (verbose) print*,'Econ ',econ

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 16. Calculate daily condensation (wc), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 68, Documentation
    out_evap%cn = 1000.0 * econ * abs(rnn_d)

    ! consistency check
    if (verbose) print*,'daily condensation (mm) ',out_evap%cn

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 17. Estimate daily EET (eet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 70, Documentation
    out_evap%eet = 1000.0*econ*(out_evap%rn)

    ! consistency check
    if (verbose) print*,'out_evap%EET (mm) ', out_evap%eet

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 18. Estimate daily PET (pet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 72, Documentation
    out_evap%pet = (1.0+kw)*out_evap%eet

    ! consistency check
    if (verbose) print*,'daily out_evap%PET (mm) ',out_evap%pet

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx = 1000.0*3600.0*(1.0+kw)*econ

    ! consistency check
    if (verbose) print*,'variable substitute rx ',rx

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 20. Calculate the intersection hour angle (hi), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cos_hi = sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv

    if (cos_hi >= 1.0) then
      ! Supply exceeds demand:
      hi = 0.0
    elseif (cos_hi <= -1.0) then
      ! Supply limits demand everywhere:
      hi = 180.0
    else
      hi = degrees( acos(cos_hi) )
    endif

    ! consistency check - XXX PROBLEM: THIS LEADS TO DIFFERENCE WITH OTHER VERSIONS XXX
    if (verbose) print*,'intersection hour angle ', hi

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    ! 21. Estimate daily out_evap%AET (out_evap%aet), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 81, Documentation
    out_evap%aet = (24.0/pi)*(radians(sw*hi) + rx*rw*rv*(dgsin(hn) - dgsin(hi)) + radians((rx*rw*ru - rx*rnl)*(hn - hi)))

    ! consistency check
    if (verbose) print*,'daily out_evap%AET (mm) ', out_evap%aet

    !-------------------------------------------------------------   
    ! Refs: Allen, R.G. (1996), Assessing integrity of weather data for 
    !         reference evapotranspiration estimation, Journal of Irrigation
    !         and Drainage Engineering, vol. 122, pp. 97--106.
    !       Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998), 
    !         'Meteorological data,' Crop evapotranspiration - Guidelines for 
    !         computing crop water requirements - FAO Irrigation and drainage 
    !         paper 56, Food and Agriculture Organization of the United 
    !         Nations, online: http://www.fao.org/docrep/x0490e/x0490e07.htm
    !       Berger, A.L. (1978), Long-term variations of daily insolation and 
    !         quarternary climatic changes, Journal of Atmospheric Sciences, 
    !         vol. 35, pp. 2362--2367.
    !       Berger, A.L., M.F. Loutre, and C. Tricot (1993), Insolation and 
    !         Earth's orbital periods, J. Geophys. Res., 98, 10341--10362.
    !       Duffie, J. A. and W. A. Beckman (1991). Solar engineering of 
    !         thermal processes. 4th ed. New Jersey: John Wiley and Sons
    !       Federer (1982), Transpirational supply and demand: plant, soil, 
    !         and atmospheric effects evaluated by simulation, Water 
    !         Resources Research, vol. 18, no. 2, pp. 355--362.
    !       Ge, S., R.G. Smith, C.P. Jacovides, M.G. Kramer, R.I. Carruthers 
    !         (2011), Dynamics of photosynthetic photon flux density (PPFD) 
    !         and estimates in coastal northern California, Theoretical and 
    !         Applied Climatology, vol. 105, pp. 107--118.
    !       Henderson-Sellers, B. (1984), A new formula for latent heat of 
    !         vaporization of water as a function of temperature, Quarterly 
    !         Journal of the Royal Meteorological Society 110, pp. 1186–1190
    !       Linacre (1968), Estimating the net-radiation flux, Agricultural 
    !         Meteorology, vol. 5, pp. 49--63.
    !       Prentice, I.C., M.T. Sykes, W. Cramer (1993), A simulation model 
    !         for the transient effects of climate change on forest 
    !         landscapes, Ecological Modelling, vol. 65, pp. 51--70.
    !       Priestley, C.H.B. and R.J. Taylor (1972), On the assessment of 
    !         surface heat flux and evaporation using large-scale parameters, 
    !         Monthly Weather Review, vol. 100 (2), pp. 81--92.
    !       Spencer, J. W. (1971), Fourier series representation of the 
    !         position of the sun, Search, vol. 2, p. 172.
    !       Stine, W. B. and M. Geyer (2001). “Power from the Sun”. 
    !         online: http://www.powerfromthesun.net/Book/chapter03/chapter03
    !       Wetherald, R.T., S. Manabe (1972), Response to joint ocean-
    !         atmosphere model to the seasonal variation of the solar 
    !         radiation, Monthly Weather Review, vol. 100 (1), pp. 42--59.
    !       Woolf, H. M. (1968). On the computation of solar evaluation 
    !         angles and the determination of sunrise and sunset times. 
    !         Tech. rep. NASA-TM-X-164. National Aeronautics and Space 
    !         Administration (NASA).
    !-------------------------------------------------------------   

  end function evap


  function dgcos( x )
    !----------------------------------------------------------------   
    ! Calculates the cosine of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real, intent(out) :: dgcos ! cosine value of x when x is in degrees

    dgcos = dcos(x*pi/180.0)

  end function dgcos


  function dgsin( x )
    !----------------------------------------------------------------   
    ! Calculates the sinus of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real, intent(out) :: dgsin ! sinus value of x when x is in degrees

    dgsin = dsin(x*pi/180.0)

  end function dgsin


  function degrees( x )
    !----------------------------------------------------------------   
    ! Returns corresponding degrees if x is given in radians
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real, intent(out) :: degrees

    degrees = x*180.0/pi

  end function degrees


  function radians( x )
    !----------------------------------------------------------------   
    ! Returns corresponding radians if x is given in degrees
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real, intent(out) :: radians

    radians = x*pi/180.0

  end function radians


  function berger_tls( day, ndayyear ) result( out_berger )
    !----------------------------------------------------------------   
    ! Returns true anomaly and true longitude for a given day
    ! Reference: Berger, A. L. (1978), Long term variations of daily 
    ! insolation and quaternary climatic changes, J. Atmos. Sci., 35, 
    ! 2362-2367.
    !----------------------------------------------------------------   
    ! arguments
    integer, intent(in) :: day   ! day of the year
    integer, intent(in) :: ndayyear

    ! function return variable
    type(outtype_berger) :: out_berger

    ! local variables
    real :: anm, ranm, anv, ranv
    real :: dlamm                ! Mean longitude for day of year
    real :: xee, xec, xse        ! variable substitutes
    real :: xlam                 ! Mean longitude for vernal equinox
    real :: tmp1, tmp2, tmp3     ! variable substitutes

    ! Variable substitutes:
    xee = ke**2 
    xec = ke**3
    xse = dsqrt(1.0 - xee)

    ! Mean longitude for vernal equinox:
    tmp1 = (ke/2.0 + xec/8.0)*(1.0 + xse)*dgsin(komega)
    tmp2 = xee/4.0*(0.5 + xse)*dgsin(2.0*komega)
    tmp3 = xec/8.0*(1.0/3.0 + xse)*dgsin(3.0*komega)
    xlam = tmp1 - tmp2 + tmp3
    xlam = degrees(2.0*xlam)

    ! Mean longitude for day of year:
    dlamm = xlam + (day - 80.0)*(360.0/ndayyear)

    ! Mean anomaly:
    anm = dlamm - komega
    ranm = radians(anm)

    ! True anomaly:
    ranv = (ranm + (2.0*ke - xec/4.0)*dsin(ranm) + 5.0/4.0*xee*dsin(2.0*ranm) + 13.0/12.0*xec*dsin(3.0*ranm))
    anv = degrees(ranv)

    ! True longitude:
    out_berger%lambda = anv + komega
    if (out_berger%lambda < 0.0) then
        out_berger%lambda = out_berger%lambda + 360.0
    else if (out_berger%lambda > 360.0) then
        out_berger%lambda = out_berger%lambda - 360.0
    endif

    ! True anomaly:
    out_berger%nu = (out_berger%lambda - komega)
    if (out_berger%nu < 0.0) then
        out_berger%nu = out_berger%nu + 360.0
    endif

  end function berger_tls


  function get_julian_day( yr, moy, dom ) result( out_julian_day )
    !----------------------------------------------------------------   
    ! Converts Gregorian date (year, month, day) to Julian 
    ! Ephemeris Day
    ! Reference:  Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," 
    ! Astronomical Algorithms
    ! xxx zero indexing in Python => usage of this function ok for Fortran? xxx
    !----------------------------------------------------------------   
    ! arguments
    integer, intent(in) :: yr      ! year
    integer, intent(in) :: moy     ! month of year
    integer, intent(in) :: dom     ! day of month

    ! local variables
    integer :: my_yr
    integer :: my_moy
    integer :: a, b

    ! function return value
    integer :: out_julian_day                   ! Julian Ephemeris Day

    if (moy <= 2) then
      my_yr  = yr - 1
      my_moy = moy + 12
    else
      my_yr  = yr
      my_moy = moy
    endif

    a = int(real(my_yr)/real(100))
    b = 2 - a + int(real(a)/real(4))
    out_julian_day=int(real(int(365.25*real(my_yr+4716)))+real(int(30.6001*real(my_moy+1)))+real(dom)+real(b)-1524.5)  ! xxx use moy instead of (moy+1) without zero-indexing as in Python?

  end function get_julian_day


  function get_sat_slope( tc ) result( out_sat_slope )
    !----------------------------------------------------------------   
    ! Calculates the slope of the sat pressure temp curve, Pa/K
    ! Ref:      Eq. 13, Allen et al. (1998)
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real :: out_sat_slope  ! slope of the sat pressure temp curve, Pa/K

    out_sat_slope = (17.269)*(237.3)*(610.78)*(exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2))

  end function get_sat_slope


  function get_enthalpy_vap( tc ) result( out_enthalpy_vap )
    !----------------------------------------------------------------   
    ! Calculates the enthalpy of vaporization, J/kg
    ! Ref:      Eq. 8, Henderson-Sellers (1984)
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real ::  out_enthalpy_vap ! enthalpy of vaporization, J/kg

    out_enthalpy_vap = 1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2

  end function get_enthalpy_vap


  function elv2pres( alt ) result( press )
    !----------------------------------------------------------------   
    ! Calculates atm. pressure for a given elevation
    ! Ref:      Allen et al. (1998)
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: alt ! elevation above sea level, m

    ! function return value
    real ::  press ! atm. pressure for a given elevation

    press = kPo*(1.0 - kL*alt/kTo)**(kG*kMa/(kR*kL))

  end function elv2pres


  function get_density_h2o( tc, press ) result( density_h2o )
    !----------------------------------------------------------------   
    ! Calculates density of water at a given temperature and pressure
    ! Ref: Chen et al. (1977)
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: tc     ! air temperature (degrees C)
    real, intent(in) :: press  ! atmospheric pressure (Pa)

    ! local variables
    real :: po, ko, ca, cb
    real :: pbar               ! atmospheric pressure (bar)

    ! function return value
    real :: density_h2o  ! density of water, kg/m^3

    ! Calculate density at 1 atm:
    po = (&
             0.99983952&
             + 6.788260e-5  *tc&
             - 9.08659e-6   *tc*tc&
             + 1.022130e-7  *tc*tc*tc  &
             - 1.35439e-9   *tc*tc*tc*tc&
             + 1.471150e-11 *tc*tc*tc*tc*tc&
             - 1.11663e-13  *tc*tc*tc*tc*tc*tc&
             + 5.044070e-16 *tc*tc*tc*tc*tc*tc*tc&
             - 1.00659e-18  *tc*tc*tc*tc*tc*tc*tc*tc&
         )

    ! Calculate bulk modulus at 1 atm:
    ko = (&
             19652.17&
             + 148.1830   *tc&
             - 2.29995    *tc*tc&
             + 0.01281    *tc*tc*tc&
             - 4.91564e-5 *tc*tc*tc*tc&
             + 1.035530e-7*tc*tc*tc*tc*tc&
         )

    ! Calculate temperature dependent coefficients:
    ca = (&
             3.26138&
             + 5.223e-4  *tc&
             + 1.324e-4  *tc*tc&
             - 7.655e-7  *tc*tc*tc&
             + 8.584e-10 *tc*tc*tc*tc&
         )
    cb = (&
             7.2061e-5&
             - 5.8948e-6  *tc&
             + 8.69900e-8 *tc*tc&
             - 1.0100e-9  *tc*tc*tc&
             + 4.3220e-12 *tc*tc*tc*tc&
         )

    ! Convert atmospheric pressure to bar (1 bar = 100000 Pa)
    pbar = (1.0e-5)*press
    
    ! ! consistency check
    ! if (verbose) print*,'atmospheric pressure, bar', pbar

    density_h2o = 1000.0*po*(ko + ca*pbar + cb*pbar**2.0)/(ko + ca*pbar + cb*pbar**2.0 - pbar)

  end function get_density_h2o


  function get_psychro( tc, press ) result( psychro )
    !----------------------------------------------------------------   
    ! Calculates the psychrometric constant for a given temperature and pressure
    ! Ref: Allen et al. (1998); Tsilingiris (2008) 
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C
    real, intent(in) :: press  ! atmospheric pressure, Pa

    ! local variables
    real :: lv  ! latent heat of vaporization (J/kg)
    real :: cp

    ! function return value
    real :: psychro  ! psychrometric constant, Pa/K

    ! Calculate the specific heat capacity of water, J/kg/K
    ! Eq. 47, Tsilingiris (2008)
    cp = 1.0e3*(&
               1.0045714270&
             + 2.050632750e-3  *tc&
             - 1.631537093e-4  *tc*tc&
             + 6.212300300e-6  *tc*tc*tc&
             - 8.830478888e-8  *tc*tc*tc*tc&
             + 5.071307038e-10 *tc*tc*tc*tc*tc&
            )

    ! Calculate latent heat of vaporization, J/kg
    lv = get_enthalpy_vap( tc )

    ! Calculate psychrometric constant, Pa/K
    ! Eq. 8, Allen et al. (1998)
    psychro = cp*kMa*press/(kMv*lv)

  end function get_psychro


  subroutine initdaily
    !----------------------------------------------------------------   
    ! Initialises daily output variables
    !----------------------------------------------------------------   
    outdra(:)   = 0.0
    outdrn(:)   = 0.0
    outdppfd(:) = 0.0
    outdcn(:)   = 0.0
    outdsm(:)   = 0.0
    outdro(:)   = 0.0
    outdeet(:)  = 0.0
    outdpet(:)  = 0.0
    outdaet(:)  = 0.0

  end subroutine initdaily


  subroutine initmonthly
    !----------------------------------------------------------------   
    ! Initialises monthly output variables
    !----------------------------------------------------------------   
    outmeet(:) = 0.0
    outmpet(:) = 0.0
    outmaet(:) = 0.0
    outmcpa(:) = 0.0
    outmcwd(:) = 0.0
    outmppfd(:)  = 0.0

  end subroutine initmonthly 


  subroutine getout_daily( doy, moy )
    !----------------------------------------------------------------   
    ! Collects daily output variables and sums up monthly output vars
    !----------------------------------------------------------------   
    ! argument
    integer, intent(in) :: doy 
    integer, intent(in) :: moy 

    ! Save the daily totals:
    outdra(doy)   = out_evap%ra
    outdrn(doy)   = out_evap%rn
    outdppfd(doy) = out_evap%ppfd
    outdcn(doy)   = out_evap%cn
    outdeet(doy)  = out_evap%eet
    outdpet(doy)  = out_evap%pet
    outdaet(doy)  = out_evap%aet
    
    outdsm(doy)   = waterbal%sm
    outdro(doy)   = waterbal%ro

    ! Update monthly totals:
    outmeet(moy) = outmeet(moy) + out_evap%eet
    outmpet(moy) = outmpet(moy) + out_evap%pet
    outmaet(moy) = outmaet(moy) + out_evap%aet
    outmppfd(moy)= outmppfd(moy)+ out_evap%ppfd

  end subroutine getout_daily


  subroutine getout_monthly( moy )
    !----------------------------------------------------------------   
    ! Initialises monthly output variables
    !----------------------------------------------------------------   
    ! argument
    integer, intent(in) :: moy

    ! Calculate other monthly totals:
    outmcpa(moy) = outmaet(moy) / outmpet(moy)
    outmcwd(moy) = outmpet(moy) - outmaet(moy)

  end subroutine getout_monthly


  subroutine write_to_file()
    !----------------------------------------------------------------   
    ! Writes daily and monthly values to files
    !----------------------------------------------------------------   
    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    ! OPEN ASCII OUTPUT FILES FOR DAILY OUTPUT
    !----------------------------------------------------------------
    prefix = "./output/"

    ! WRTIING DAILY TOTALS
    !----------------------------------------------------------------
    print*,"writing daily totals to files"

    ! HO: daily solar irradiation, J/m2
    filnam=trim(prefix)//'ra.d.out'
    open(101,file=filnam,err=888,status='unknown')
    write(101,999) outdra

    ! HN: daily net radiation, J/m2
    filnam=trim(prefix)//'rn.d.out'
    open(102,file=filnam,err=888,status='unknown')
    write(102,999) outdrn

    ! QN: daily PPFD, mol/m2
    filnam=trim(prefix)//'ppfd.d.out'
    open(103,file=filnam,err=888,status='unknown')
    write(103,999) outdppfd

    ! CN: daily condensation water, mm
    filnam=trim(prefix)//'cn.d.out'
    open(104,file=filnam,err=888,status='unknown')
    write(104,999) outdcn

    ! WN: daily soil moisture, mm
    filnam=trim(prefix)//'sm.d.out'
    open(105,file=filnam,err=888,status='unknown')
    write(105,999) outdsm

    ! RO: daily runoff, mm
    filnam=trim(prefix)//'ro.d.out'
    open(107,file=filnam,err=888,status='unknown')
    write(107,999) outdro

    ! EQ_N: daily equilibrium ET, mm
    filnam=trim(prefix)//'eet.d.out'
    open(108,file=filnam,err=888,status='unknown')
    write(108,999) outdeet

    ! EP_N: daily potential ET, mm
    filnam=trim(prefix)//'pet.d.out'
    open(109,file=filnam,err=888,status='unknown')
    write(109,999) outdpet

    ! EA_N: daily actual ET, mm
    filnam=trim(prefix)//'aet.d.out'
    open(110,file=filnam,err=888,status='unknown')
    write(110,999) outdaet


    ! WRTIING MONTHLY TOTALS
    !----------------------------------------------------------------
    print*,"writing monthly totals to files"

    ! eq_m
    filnam=trim(prefix)//'eet.m.out'
    open(111,file=filnam,err=888,status='unknown')
    write(111,999) outmeet

    ! ep_m
    filnam=trim(prefix)//'pet.m.out'
    open(112,file=filnam,err=888,status='unknown')
    write(112,999) outmpet

    ! ea_m
    filnam=trim(prefix)//'aet.m.out'
    open(113,file=filnam,err=888,status='unknown')
    write(113,999) outmaet

    ! cpa
    filnam=trim(prefix)//'cpa.m.out'
    open(114,file=filnam,err=888,status='unknown')
    write(114,999) outmcpa

    ! cwd
    filnam=trim(prefix)//'cwd.m.out'
    open(115,file=filnam,err=888,status='unknown')
    write(115,999) outmcwd

    ! qm
    filnam=trim(prefix)//'ppfd.m.out'
    open(116,file=filnam,err=888,status='unknown')
    write(116,999) outmppfd

    return

    888    stop 'WRITE_TO_FILE: error opening output file'
    999    format (F20.8)

  end subroutine write_to_file

end module _splash
