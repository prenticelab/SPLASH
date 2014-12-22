      module modelparams
        !////////////////////////////////////////////////////////////////
        ! Module contains model parameters
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
        real, parameter :: kR = 8.3143          ! universal gas constant, J/mol/K (Allen, 1973)
        real, parameter :: kTo = 298.15         ! base temperature, K (Prentice, unpublished)
        real, parameter :: kWm = 150            ! soil moisture capacity, mm (Cramer & Prentice, 1988)
        real, parameter :: kw = 0.26            ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
        real, parameter :: komega = 283.0       ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)
        real, parameter :: pi = 3.14159

      end module modelparams


      module auxfunc
        !////////////////////////////////////////////////////////////////
        ! Module contains auxiliary functions required by stash
        !----------------------------------------------------------------   
        
        function dcos( x )
          !----------------------------------------------------------------   
          ! Calculates the cosine of an angle given in degrees
          !----------------------------------------------------------------   
          use modelparams, only: pi
          implicit none

          ! arguments
          real, intent(in) :: x  ! angle, degrees (0-360)

          ! function return value
          real, intent(out) :: dcos ! cosine value of x when x is in degrees

          dcos = cos(x*pi/180.0)

        end function dcos


        function dsin( x )
          !----------------------------------------------------------------   
          ! Calculates the sinus of an angle given in degrees
          !----------------------------------------------------------------   
          use modelparams, only: pi
          implicit none

          ! arguments
          real, intent(in) :: x  ! angle, degrees (0-360)

          ! function return value
          real, intent(out) :: dsin ! sinus value of x when x is in degrees

          dsin = sin(x*pi/180.0)

        end function dsin

        function degrees( x )
          !----------------------------------------------------------------   
          ! Returns corresponding degrees if x is given in radians
          !----------------------------------------------------------------   
          use modelparams, only: pi
          implicit none

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
          use modelparams, only: pi
          implicit none

          ! arguments
          real, intent(in) :: x  ! angle, radians

          ! function return value
          real, intent(out) :: radians

          radians = x*pi/180.0

        end function radians


        function berger_tls( day )
          !----------------------------------------------------------------   
          ! Returns true anomaly and true longitude for a given day
          ! Reference: Berger, A. L. (1978), Long term variations of daily 
          ! insolation and quaternary climatic changes, J. Atmos. Sci., 35, 
          ! 2362-2367.
          !----------------------------------------------------------------   
          use modelparams, only: ke
          use auxfunc
          implicit none

          ! arguments
          integer, intent(in) :: day   ! day of the year

          ! function return value
          real, dimension(2) :: berger_tls ! return "tuple": (my_nu = true anomaly, degrees; my_tls = true longitude, degrees)

          ! Variable substitutes:
          xee = ke**2 
          xec = ke**3
          xse = sqrt(1.0 - xee)

          ! Mean longitude for vernal equinox:
          xlam = (ke/2.0 + xec/8.0)*(1.0 + xse)*dsin(komega) - xee/4.0*(0.5 + xse)*dsin(2.0*komega) + xec/8.0*(1.0/3.0 + xse)*dsin(3.0*komega)
          xlam = degrees(2.0*xlam)

          ! Mean longitude for day of year:
          dlamm = xlam + (day - 80.0)*(360.0/kN)

          ! Mean anomaly:
          anm = dlamm - komega
          ranm = radians(anm)

          ! True anomaly:
          ranv = (ranm + (2.0*ke - xec/4.0)*sin(ranm) + 5.0/4.0*xee*sin(2.0*ranm) + 13.0/12.0*xec*sin(3.0*ranm))
          anv = degrees(ranv)

          ! True longitude:
          my_tls = anv + komega
          if (my_tls < 0.0) then
              my_tls = my_tls + 360.0
          else if (my_tls > 360.0) then
              my_tls = my_tls - 360.0
          endif

          ! True anomaly:
          my_nu = (my_tls - komega)
          if (my_nu < 0.0) then
              my_nu = my_nu + 360.0
          endif

          berger_tls = (/my_nu, my_tls/)

        end function berger_tls



        function julian_day( y, m, i )
          !----------------------------------------------------------------   
          ! Converts Gregorian date (year, month, day) to Julian 
          ! Ephemeris Day
          ! Reference:  Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," 
          ! Astronomical Algorithms
          !----------------------------------------------------------------   
          implicit none

          ! arguments
          integer, intent(in) :: y     ! year
          integer, intent(in) :: m     ! month
          integer, intent(in) :: i     ! day of month

          ! local variables
          integer :: a, b

          ! function return value
          real :: julian_day                   ! Julian Ephemeris Day

          if (m <= 2) then
            y = y - 1
            m = m + 12
          endif

          a = int(real(y)/real(100))
          b = 2 - a + int(real(a)/real(4))
          julian_day = real(int( 365.25 * real(y+4716) )) + real(int( 30.6001 * real(m+1) )) + real(i) + real(b) - 1524.5

        end function julian_day


        function sat_slope( tc )
          !----------------------------------------------------------------   
          ! Calculates the slope of the sat pressure temp curve, Pa/K
          ! Ref:      Eq. 13, Allen et al. (1998)
          !----------------------------------------------------------------   
          implicit none

          ! arguments
          real, intent(in) :: tc ! air temperature, degrees C

          ! function return value
          real :: sat_slope  ! slope of the sat pressure temp curve, Pa/K

          sat_slope = (17.269)*(237.3)*(610.78)*(exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2))

        end function sat_slope


        function enthalpy_vap( tc )
          !----------------------------------------------------------------   
          ! Calculates the enthalpy of vaporization, J/kg
          ! Ref:      Eq. 8, Henderson-Sellers (1984)
          !----------------------------------------------------------------   
          implicit none

          ! arguments
          real, intent(in) :: tc ! air temperature, degrees C

          ! function return value
          real ::  enthalpy_vap ! enthalpy of vaporization, J/kg

          enthalpy_vap = (1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2)

        end function enthalpy_vap


        function elv2pres( alt )
          !----------------------------------------------------------------   
          ! Calculates atm. pressure for a given elevation
          ! Ref:      Allen et al. (1998)
          !----------------------------------------------------------------   
          use modelparams, only: kPo, kTo, kL, kMa, kG, kR
          implicit none

          ! arguments
          real, intent(in) :: alt ! elevation above sea level, m

          ! function return value
          real ::  elv2pres ! atm. pressure for a given elevation

          elv2pres = kPo*(1.0 - kL*alt/kTo)**(kG*kMa/(kR*kL))

        end function elv2pres


        function density_h2o( tc, p )
          !----------------------------------------------------------------   
          ! Calculates density of water at a given temperature and pressure
          ! Ref: Chen et al. (1977)
          !----------------------------------------------------------------   
          implicit none

          ! arguments
          real, intent(in) :: tc ! air temperature, degrees C
          real, intent(in) :: press  ! atmospheric pressure, Pa

          ! function return value
          real :: density_h2o  ! density of water, kg/m^3

          ! Calculate density at 1 atm:
          po = (
     $         0.99983952 + 
     $         (6.788260e-5)*tc + 
     $         -(9.08659e-6)*tc*tc +
     $         (1.022130e-7)*tc*tc*tc + 
     $         -(1.35439e-9)*tc*tc*tc*tc +
     $         (1.471150e-11)*tc*tc*tc*tc*tc +
     $         -(1.11663e-13)*tc*tc*tc*tc*tc*tc + 
     $         (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc + 
     $         -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
     $     )

          ! Calculate bulk modulus at 1 atm:
          ko = (
     $         19652.17 +
     $         148.1830*tc + 
     $         -2.29995*tc*tc + 
     $         0.01281*tc*tc*tc + 
     $         -(4.91564e-5)*tc*tc*tc*tc + 
     $         (1.035530e-7)*tc*tc*tc*tc*tc
     $     )

          ! Calculate temperature dependent coefficients:
          ca = (
     $         3.26138 + 
     $         (5.223e-4)*tc + 
     $         (1.324e-4)*tc*tc + 
     $         -(7.655e-7)*tc*tc*tc + 
     $         (8.584e-10)*tc*tc*tc*tc
     $     )
          cb = (
     $         (7.2061e-5) +
     $         -(5.8948e-6)*tc + 
     $         (8.69900e-8)*tc*tc + 
     $         -(1.0100e-9)*tc*tc*tc + 
     $         (4.3220e-12)*tc*tc*tc*tc
     $     )

          ! Convert atmospheric pressure to bar (1 bar = 100000 Pa)
          pbar = (1.0e-5)*press

          density_h2o = 1000.0*po*(ko + ca*pbar + cb*pbar**2.0)/(ko + ca*pbar + cb*pbar**2.0 - pbar)

        end function density_h2o


        function psychro( tc, press )
          !----------------------------------------------------------------   
          ! Calculates the psychrometric constant for a given temperature and pressure
          ! Ref: Allen et al. (1998); Tsilingiris (2008) 
          !----------------------------------------------------------------   
          use modelparams, only: kMa, kMv
          implicit none

          ! arguments
          real, intent(in) :: tc ! air temperature, degrees C
          real, intent(in) :: press  ! atmospheric pressure, Pa

          ! local variables
          real :: lv  ! latent heat of vaporization (J/kg)

          ! function return value
          real :: psychro  ! psychrometric constant, Pa/K

          ! Calculate the specific heat capacity of water, J/kg/K
          ! Eq. 47, Tsilingiris (2008)
          cp = (1.0e3)*(
     $         1.0045714270 +
     $         (2.050632750e-3)*tc -
     $         (1.631537093e-4)*tc*tc +
     $         (6.212300300e-6)*tc*tc*tc -
     $         (8.830478888e-8)*tc*tc*tc*tc +
     $         (5.071307038e-10)*tc*tc*tc*tc*tc
     $        )

          ! Calculate latent heat of vaporization, J/kg
          lv = enthalpy_vap(tc)

          ! Calculate psychrometric constant, Pa/K
          ! Eq. 8, Allen et al. (1998)
          psychro = cp*kMa*press/(kMv*lv)

        end function psychro

      end module auxfunc


      program stash
      !////////////////////////////////////////////////////////////////
      ! main program 
      !----------------------------------------------------------------
      use modelparams, only: keps, ke, pi, kGsc, kc, kd, kfFEC, kalb_vis
      use auxfunc 

      implicit none

      lon = 2.0
      lat = 46.5
      elv = 300

      contains

      subroutine evap( lon, lat, n, ppfd_d, eet_d, pet_d, aet_d, wc, elv, year, sf, tc, sw )
        !----------------------------------------------------------------   
        ! This subroutine calculates daily radiation and evapotranspiration
        ! quantities
        ! - daily PPFD (ppfd_d), mol/m^2
        ! - daily EET (eet_d), mm
        ! - daily PET (pet_d), mm
        ! - daily AET (aet_d), mm
        ! - daily condensation (wc), mm
        !-------------------------------------------------------------  
        implicit none

        ! ARGUMENTS
        ! inputs
        real, intent(in) :: lon              ! longitude, degrees
        real, intent(in) :: lat              ! latitude, degrees
        integer, intent(in) :: n             ! day of the year
        real, intent(in) :: elv              ! elevation, metres
        integer, intent(in) :: year          ! year
        real, intent(in), optional :: sf     ! fraction of sunshine hours 
        real, intent(in), optional :: tc     ! mean daily air temperature, C
        real, intent(in), optional :: sw     ! evaporative supply rate, mm/hr
        ! outputs
        real, intent(out) :: ppfd_d          ! daily PPFD (ppfd_d), mol/m^2
        real, intent(out) :: eet_d           ! daily EET (eet_d), mm
        real, intent(out) :: pet_d           ! daily PET (pet_d), mm
        real, intent(out) :: aet_d           ! daily AET (aet_d), mm
        real, intent(out) :: wc              ! daily condensation (wc), mm

        ! LOCAL VARIABLES
        integer :: kN                        ! number of days in year
        real, dimension(2) :: out_berger     ! temporary var to store output of berger_tls
        real :: my_nu                        ! heliocentric longitudes
        real :: my_lambda                    ! heliocentric longitudes
        real :: my_rho
        real :: dr                           ! distance factor
        real :: delta                        ! declination angle 
        real :: ru                           ! variable substitute for u
        real :: rv                           ! variable substitute for v
        real :: hs                           ! sunset hour angle
        real :: ra_d                         ! daily extraterrestrial solar radiation (J/m^2)
        real :: tau_o, tau                   ! transmittivity (unitless)
        real :: rnl                          ! net longwave radiation (W/m^2)
        real :: rw                           ! variable substitute (W/m^2)
        real :: hn                           ! net radiation cross-over hour angle
        real :: rn_d                         ! daytime net radiation, J/m^2
        real :: s                            ! slope of saturation vap press temp curve, Pa/K
        real :: pw                           ! density of water, kg/m^3
        real :: lv                           ! enthalpy of vaporization, J/kg
        real :: g                            ! psychrometric constant, Pa/K
        real :: econ                         ! Eq. 58, STASH 2.0 Documentation
        real :: rx                           ! variable substitute (mm/hr)/(W/m^2)
        real :: hi, cos_hi                   ! intersection hour angle, degrees


        ! Set default values for optional arguments
        if (.not.present(elf)) elv = 0.0
        if (.not.present(year)) year = 0
        if (.not.present(sf)) sf = 1.0
        if (.not.present(tc)) tc = 23.0
        if (.not.present(sw)) sw = 1.0

        ! Error handle and assign required public variables:
        if (lon>180.0 .or. lon<(-180.0)) then
          print*,"Longitude outside range of validity (-180 to 180)!"
          stop
        endif

        if (lat>90.0 .or. lat<(-90.0)) then
          print*,"Latitude outside range of validity (-90 to 90)!"
          stop
        endif

        if (n<1 .or. n>366) then
          print*,"Day of year outside range of validity (1 to 366)!"
          stop
        endif


        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 1. Calculate number of days in year (kN), days
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (y==0) then
          kN = 365
          user_year = 2001
        else
          kN = julian_day((user_year+1),1,1) - julian_day(user_year, 1, 1)
        endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Berger (1978)
        out_berger = berger_tls( n )
        my_nu      = out_berger(1)
        my_lambda  = out_berger(2)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 3. Calculate distance factor (dr), unitless
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Berger et al. (1993)
        my_rho = (1.0 - ke**2)/(1.0 + ke * dcos( my_nu ))
        dr = (1.0/my_rho)**2

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 4. Calculate declination angle (delta), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Woolf (1968)
        delta = arcsin( dsin( my_lambda ) * dsin( keps ) )
        delta = degrees( delta )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 5. Calculate variable substitutes (u and v), unitless
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = dsin(delta) * dsin(lat)
        rv = dcos(delta) * dcos(lat)

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
          hs = acos(-1.0*ru/rv)
          hs = degrees(hs)
        endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 1.10.3, Duffy & Beckman (1993)
        ra_d = (86400.0/pi)*kGsc*dr*(degrees(ru)*hs + rv * dsin(hs))

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 8. Calculate transmittivity (tau), unitless
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 11, Linacre (1968)
        tau_o = (kc + kd*sf)

        ! Eq. 2, Allen (1996)
        tau = tau_o*(1.0 + (2.67e-5)*elv)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 9. Calculate daily PPFD (ppfd_d), mol/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 57, STASH 2.0 Documentation
        ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra_d

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 10. Estimate net longwave radiation (rnl), W/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        rnl = (kb + (1.0 - kb)*sf)*(kA - tc)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 11. Calculate variable substitute (rw), W/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = (1.0-kalb_sw)*tau*kGsc*dr

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 12. Calculate net radiation cross-over hour angle (hn), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ((self.rnl - rw*ru)/(rw*rv) >= 1.0) then
          ! Net radiation negative all day
          hn = 0.0
        elseif ((rnl - rw*ru)/(rw*rv) <= -1.0) then
          ! Net radiation positive all day
          hn = 180.0
        else
          hn = degrees( acos((rnl - rw*ru)/(rw*rv)) )
        endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 13. Calculate daytime net radiation (rn_d), J/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 53, STASH 2.0 Documentation
        rn_d = (86400.0/pi) * (hn*(pi/180.0)*(rw*ru - rnl) + rw*rv*dsin(hn))

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 14. Calculate nighttime net radiation (rnn_d), J/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 56, STASH 2.0 Documentation
        rnn_d = (86400.0/pi)*(radians(rw*ru*(hs-hn)) + rw*rv*(dsin(hs)-dsin(hn)) + rnl*(pi - 2.0*radians(hs) + radians(hn)))

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 15. Calculate water-to-energy conversion (econ), m^3/J
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Slope of saturation vap press temp curve, Pa/K
        s = sat_slope(tc)
        ! Enthalpy of vaporization, J/kg
        lv = enthalpy_vap(tc)
        ! Density of water, kg/m^3
        pw = density_h2o(tc, elv2pres(elv))
        ! Psychrometric constant, Pa/K
        g = psychro(tc, elv2pres(elv))

        ! Eq. 58, STASH 2.0 Documentation
        econ = s/(lv*pw*(s + g))

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 16. Calculate daily condensation (wc), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 68, STASH 2.0 Documentation
        wc = 1000.0 * econ * abs(rnn_d)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 17. Estimate daily EET (eet_d), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 70, STASH 2.0 Documentation
        eet_d = 1000.0*econ*(rn_d)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 18. Estimate daily PET (pet_d), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 72, STASH 2.0 Documentation
        pet_d = (1.0+kw)*eet_d

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = 1000.0*3600.0*(1.0+kw)*econ

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
          hi = degrees(acos(cos_hi))
        endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 21. Estimate daily AET (aet_d), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 81, STASH 2.0 Documentation
        aet_d = (24.0/pi)*(radians(sw*hi) + rx*rw*rv*(dsin(hn) - dsin(hi)) + radians((rx*rw*ru - rx*rnl)*(hn - hi)))


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

      end subroutine evap


      subroutine stash()




      end subroutine stash






      end program
