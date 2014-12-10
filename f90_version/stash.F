      module modelparams
        !////////////////////////////////////////////////////////////////
        ! Module contains model parameters
        !----------------------------------------------------------------      
        real, parameter :: kA = 107             ! constant for Rnl (Monteith & Unsworth, 1990)
        real, parameter :: ka = 1.49598e8       ! semi-major axis, km (Allen, 1973)
        real, parameter :: kalb_sw = 0.17       ! shortwave albedo (Federer, 1968)
        real, parameter :: kalb_vis = 0.03      ! visible light albedo (Sellers, 1985)
        real, parameter :: kb = 0.20            ! constant for Rnl (Linacre, 1968)
        real, parameter :: kc = 0.25            ! cloudy transmittivity (Linacre, 1968)
        real, parameter :: kCw = 1.05           ! supply constant, mm/hr (Federer, 1982)
        real, parameter :: kd = 0.50            ! angular coefficient of transmittivity (Linacre, 1968)
        real, parameter :: ke = 0.0167          ! eccentricity for 2000 CE (Berger, 1978)
        real, parameter :: keps = 23.44         ! obliquity for 2000 CE, degrees (Berger, 1978)
        real, parameter :: kerror = -9999       ! error value
        real, parameter :: kfFEC = 2.04         ! from flux to energy conversion, umol/J (Meek et al., 1984)
        real, parameter :: kG = 9.80665         ! gravitational acceleration, m/s^2 (Allen, 1973)
        real, parameter :: kGM = 1.32712e11     ! standard gravity of the sun, km^3/s^2
        real, parameter :: kGsc = 1360.8        ! solar constant, W/m^2 (Kopp & Lean, 2011)
        real, parameter :: kL = 0.0065          ! temperature lapse rate, K/m (Cavcar, 2000)
        real, parameter :: kMa = 0.028963       ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
        real, parameter :: kMv = 0.01802        ! molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
        real, parameter :: kPo = 101325         ! standard atmosphere, Pa (Allen, 1973)
        real, parameter :: kR = 8.3143          ! universal gas constant, J/mol/K (Allen, 1973)
        real, parameter :: ksb = 5.670373e-8    ! Stefan-Boltzman constant, W/m^2/K^4
        real, parameter :: kTo = 298.15         ! base temperature, K (Prentice, unpublished)
        real, parameter :: kWm = 150            ! soil moisture capacity, mm (Cramer & Prentice, 1988)
        real, parameter :: kw = 0.26            ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
        real, parameter :: komega = 283.0       ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)

      end module modelparams


      module auxfunc
        !////////////////////////////////////////////////////////////////
        ! Module contains auxiliary functions required by stash
        !----------------------------------------------------------------   

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

          ! FUNCTION RETURN VALUES
          real :: jde                   ! Julian Ephemeris Day

          if (m <= 2) then
            y = y - 1
            m = m + 12
          endif

          a = int(real(y)/real(100))
          b = 2 - a + int(real(a)/real(4))
          jde = real(int( 365.25 * real(y+4716) )) + real(int( 30.6001 * real(m+1) )) + real(i) + real(b) - 1524.5

        end function julian_day


        function map_days( n, y )
          !----------------------------------------------------------------  
          ! Computes earth's longitude relative to the vernal equinox 
          ! for a given day of the year.
          ! Depends on functions:
          !        - earth_period
          !        - earth_velocity
          !        - equinox
          !        - full_kepler
          !----------------------------------------------------------------
          use modelparams, only: komega   
          implicit none

          ! ARGUMENTS
          integer, intent(in) :: n     ! day of the year
          integer, intent(in) :: y     ! year

          ! LOCAL VARIABLES
          real,dimension(366) :: lon_nu  ! longitude vector
          real :: wp                   ! angle of the vernal equinox w.r.t. the perihelion
          integer :: i

          ! FUNCTION RETURN VALUES
          real :: lamda_doy            ! longitude relative to the vernal equinox  d

          ! Create longitude field and compute the days for orbit:
          lon_nu = (/(i, i=0,365,1)/)*(360.0/365.0)
          day_nu = full_kepler(lon_nu)
          
          ! Compute the angle of the vernal equinox w.r.t. the perihelion
          ! i.e., the explementary angle of komega:
          wp = (360.0 - komega)
          
xxx

          ! Calculate the length of time it takes earth to travel from the
          ! perihelion (t=0) to the vernal equinox:
          if (wp==180.0) then
            wp = wp + 1.e-6
          endif
          days_to_ve = full_kepler(wp)[0]  ! xxx first value of bituple
          
          ! Get day of year of vernal equinox
          if (y==0) then
            day_ve = 80.0
          else
            day_ve = equinox(y)
          endif

          ! Calculate the offset between days to and day of vernal equinox:
          offset = (day_ve - days_to_ve)
          
          ! Calculate the calendar days and set between 0 and kN:
          calendar_days = (day_nu + offset)
          calendar_days[numpy.where(calendar_days >= self.kN)] -= self.kN
          
          ! Check to see if n is listed in calendar:
          if n in calendar_days:
              icalendar = numpy.where(calendar_days==n)[0][0]
              nu_doy = lon_nu[icalendar]
          else:
              ! Find the calendar day the precedes doy:
              calendar_temp = numpy.sort(numpy.append(calendar_days, [n,]))
              dbefore = calendar_temp[numpy.where(calendar_temp==n)[0][0]-1]
              
              ! Get the index of the preceding day:
              ibefore = numpy.where(calendar_days == dbefore)[0][0]
              
              ! Get the angular velocity for the longitude of the preceding day:
              vbefore = self.earth_velocity(lon_nu[ibefore])
              
              ! Calculate the delta time
              dt = (n - dbefore)
              
              ! Calculate the new longitude, degrees:
              nu_doy = lon_nu[ibefore] + (vbefore*dt)*180.0/numpy.pi
          
          if nu_doy >= 360:
              nu_doy -= 360
          
          ! Convert nu to lambda:
          lamda_doy = nu_doy + komega
          if lamda_doy >= 360:
              lamda_doy -= 360
          
          return (nu_doy, lamda_doy)


    def full_kepler(self, lon):
        """
        Name:     EVAP.full_kepler
        Input:    float, longitude w.r.t. the perihelion, deg (lon)
        Output:   float, days traveled from perihelion (t)
        Features: Returns the days traveled since the earth past the perihelion
        Depends:  -ke (ellipticity)
        Ref:      Kepler's Second Law
        """
        # Make lon into numpy array, if it is not already:
        if isinstance(lon, numpy.float) or isinstance(lon, numpy.int):
            lon = numpy.array([lon,])
        elif not isinstance(lon, numpy.ndarray):
            lon = numpy.array(lon)
        #
        xee = ke**2
        xte = ke**3
        xse = numpy.sqrt(1.0 - xee)
        xco = self.dcos(lon) + 1
        #
        A = (0.5*self.kN/numpy.pi)*(1 - xee)**1.5
        B = 2*numpy.arctan((ke - 1.)*self.dsin(lon)/(xse*xco))
        C = xse*(xee - 1)
        D = 2.*ke*self.dsin(lon)
        E = xco*((xte - xee - ke + 1)*self.dsin(lon)**2/xco**2 - 
                 xte - xee + ke + 1)
        kepler_t = A*(B/C - D/E)
        #
        # Correct negative points due to inverse tan function:
        neg_points = numpy.where(lon > 180)[0]
        kepler_t[neg_points] -= 2*numpy.pi*A/C
        #
        return (kepler_t)


    def equinox(self, year, opt=0):
        """
        Name:     EVAP.equinox
        Input:    - int, year (year)
                  - int, option (opt)
                    0: vernal equinox     1: summer solstice
                    2: autumnal equinox   3: winter solstice
        Output:   float, day of the year
        Features: Calculates the day of the year on which seasonal dates fall
        Depends:  julian_day
        Ref:      J. Meeus (1991), Ch.26 "Equinoxes and solstices," 
                  Astronomical Algorithms
        """
        # Table 26.C (Meeus, 1991)
        periodic_terms = numpy.array([
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
            ])
        #
        # Table 26.A (Meeus, 1991)
        jde_table_a = numpy.array([
            numpy.array([1721139.29189, 365242.13740,  0.06134,  
                         0.00111, -0.00071]), # March equinox
            numpy.array([1721233.25401, 365241.72562, -0.05323,  
                         0.00907,  0.00025]), # June solstice
            numpy.array([1721325.70455, 365242.49558, -0.11677, 
                         -0.00297,  0.00074]), # September equinox
            numpy.array([1721414.39987, 365242.88257, -0.00769, 
                         -0.00933, -0.00006])  # December solstice
            ])
        #
        # Table 26.B (Meeus, 1991)
        jde_table_b = numpy.array([
            numpy.array([2451623.80984, 365242.37404,  0.05169, 
                         -0.00411, -0.00057]), # March equinox
            numpy.array([2451716.56767, 365241.62603,  0.00325,  
                         0.00888, -0.00030]), # June solstice
            numpy.array([2451810.21715, 365242.01767, -0.11575,  
                         0.00337,  0.00078]), # September equinox
            numpy.array([2451900.05952, 365242.74049, -0.06223, 
                         -0.00823,  0.00032])  # December solstice
            ])
        #
        if year < 1000:
            # Use Table 26.A for years -1000 to 1000
            jde_table = jde_table_a
            y = year/1000.0
        else:
            # Use Table 26.B for years 1000 to 3000
            jde_table = jde_table_b
            y = (year-2000.0)/1000.0
        #
        # Calculate the mean equinox based on Table 26.A or 26.B:
        jde0 = (
        jde_table[opt,0] +
        jde_table[opt,1]*y +
        jde_table[opt,2]*y*y +
        jde_table[opt,3]*y*y*y +
        jde_table[opt,4]*y*y*y*y
        )
        #
        # Calculate the other three terms:
        t = (jde0 - 2451545.0)/36525.0
        w = (35999.373*t) - 2.47
        dl = (
            1.0 + 
            (0.0334*self.dcos(w)) + 
            (0.0007*self.dcos(2*w))
        )
        #
        # Compute the sum of the 24 periodic terms:
        s = 0
        j = 0
        for i in xrange(24):
            s += periodic_terms[j]*self.dcos(
                (periodic_terms[j+1] + (periodic_terms[j+2]*t))
            )
            j += 3
        #
        # Calculate the JDE (Julian Ephemeris Day) of the equinox/solstice:
        jde = jde0 + ((s*0.00001)/dl)
        #
        # Calcuate the JDE for the first day of this year:
        jde_year = self.julian_day(year, 1, 1)
        #
        # Return the day of year:
        return (jde - jde_year + 1)
   
xxx

      end module auxfunc


      program stash
      !////////////////////////////////////////////////////////////////
      ! main program 
      !----------------------------------------------------------------
      use modelparams

      implicit none

      lon = 2.0
      lat = 46.5
      elv = 300

      call stash_thiscell( lon, lat, n, elv=0.0, y=0, sf=1.0, tc=23.0, sw=1.0, drm='loutre', lamm = 'kepler', delm = 'loutre')

      contains

      subroutine stash_thiscell( lon, lat, n, elv, year, sf, tc, sw, drm, lamm, delm )

      implicit none

      ! ARGUMENTS
      real, intent(in) :: lon                           ! longitude, degrees
      real, intent(in) :: lat                           ! latitude, degrees
      integer, intent(in) :: n                          ! day of the year
      real, intent(in) :: elv                           ! elevation, metres
      integer, intent(in) :: year                       ! year
      real, intent(in), optional :: sf                  ! fraction of sunshine hours 
      real, intent(in), optional :: tc                  ! mean daily air temperature, C
      real, intent(in), optional :: sw                  ! evaporative supply rate, mm/hr
      character(len=*), intent(in), optional :: drm     ! distance method 
      character(len=*), intent(in), optional :: lamm    ! lambda method
      character(len=*), intent(in), optional :: delm    ! delta method

      ! LOCAL VARIABLES
      integer :: kN            ! number of days in year

      ! Set default values for optional arguments
      if (.not.present(elf)) elv = 0.0
      if (.not.present(year)) year = 0
      if (.not.present(sf)) sf = 1.0
      if (.not.present(tc)) tc = 23.0
      if (.not.present(sw)) sw = 1.0
      if (.not.present(drm)) drm = 'loutre'
      if (.not.present(lamm)) lamm = 'kepler'
      if (.not.present(delm)) delm = 'loutre'

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 1. Calculate number of days in year (kN), days
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (y == 0) then
        kN = 365
      else
        kN = ( julian_day( (year+1), 1, 1 ) - julian_day( year, 1, 1 ) )
      endif

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (lamm == 'kepler') then
        ! Kepler Method:

xxx started working here, but got distracted following up function map_days ...

        self.my_nu, self.my_lambda = self.map_days(self.user_day, 
                                                   self.user_year)

      else if (lamm == 'woolf') then
        ! Woolf Method:
        woolf_b = (self.user_day - 1.0)*(360.0/self.kN)
        self.my_lambda = (
            279.9348 + woolf_b + 1.914827*self.dsin(woolf_b) - 
            0.079525*self.dcos(woolf_b) + 0.019938*self.dsin(2*woolf_b) - 
            0.00162*self.dcos(2*woolf_b)
        )
        if self.my_lambda < 0:
            self.my_lambda += 360.0
        elif self.my_lambda > 360:
            self.my_lambda -= 360.0
        self.my_nu = (self.my_lambda - komega)
        if self.my_nu < 0:
            self.my_nu += 360.0

      else if (lamm == 'berger') then
        ! Berger'78 Method:
        xee = ke**2 
        xec = ke**3
        xse = numpy.sqrt(1.0 - xee)
        xlam = ((ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(komega) - 
            xee/4.0*(0.5 + xse)*self.dsin(2.0*komega) + 
            xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*komega))
        xlam = numpy.degrees(2.0*xlam)
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        anm = dlamm - komega
        ranm = numpy.radians(anm)
        ranv = (ranm + (2.0*ke - xec/4.0)*numpy.sin(ranm) + 
            5.0/4.0*xee*numpy.sin(2.0*ranm) + 
            13.0/12.0*xec*numpy.sin(3.0*ranm))
        anv = numpy.degrees(ranv)
        self.my_lambda = anv + komega
        if self.my_lambda < 0:
            self.my_lambda += 360.0
        elif self.my_lambda > 360:
            self.my_lambda -= 360.0
        self.my_nu = (self.my_lambda - komega)
        if self.my_nu < 0:
            self.my_nu += 360.0
      else
        print*,"Heliocentric longitude method not recognized!"
        stop
      endif

      end subroutine stash_thiscell

      end program
