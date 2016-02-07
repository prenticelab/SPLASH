module modelparams
  !////////////////////////////////////////////////////////////////
  ! Module contains model parameters
  !----------------------------------------------------------------   
  implicit none

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

  ! Chose whether to use monthly or daily input data
  logical, parameter :: use_daily_input = .true.

  ! Chose whether to do the standard consistency check (see Wiki). This overrides 'use_daily_input'.
  logical, parameter :: do_consistency_check = .true.

end module modelparams


module outputvars
  !////////////////////////////////////////////////////////////////
  ! Module contains output variables
  !----------------------------------------------------------------   
  implicit none

  ! daily totals
  real, dimension(366) :: dho             ! daily solar irradiation, J/m2
  real, dimension(366) :: dhn             ! daily net radiation, J/m2
  real, dimension(366) :: dqn             ! daily PPFD, mol/m2
  real, dimension(366) :: dcn             ! daily condensation water, mm
  real, dimension(366) :: dwn             ! daily soil moisture, mm
  real, dimension(366) :: dpn             ! daily precipitation, mm
  real, dimension(366) :: dro             ! daily runoff, mm
  real, dimension(366) :: deq_n           ! daily equilibrium ET, mm
  real, dimension(366) :: dep_n           ! daily potential ET, mm
  real, dimension(366) :: dea_n           ! daily actual ET, mm

  ! monthly totals
  real, dimension(12) :: meq_m
  real, dimension(12) :: mep_m
  real, dimension(12) :: mea_m
  real, dimension(12) :: mcpa
  real, dimension(12) :: mcwd
  real, dimension(12) :: mqm

contains

  subroutine initdaily
    !----------------------------------------------------------------   
    ! Initialises daily output variables
    !----------------------------------------------------------------   
    dho(:) = 0.0
    dhn(:) = 0.0
    dqn(:) = 0.0
    dcn(:) = 0.0
    dwn(:) = 0.0
    dpn(:) = 0.0
    dro(:) = 0.0
    deq_n(:) = 0.0
    dep_n(:) = 0.0
    dea_n(:) = 0.0

  end subroutine initdaily


  subroutine initmonthly
    !----------------------------------------------------------------   
    ! Initialises monthly output variables
    !----------------------------------------------------------------   
    meq_m(:) = 0.0
    mep_m(:) = 0.0
    mea_m(:) = 0.0
    mcpa(:) = 0.0
    mcwd(:) = 0.0
    mqm(:) = 0.0

  end subroutine initmonthly 

end module outputvars


module daily_input
  !////////////////////////////////////////////////////////////////
  ! Module contains function to read daily values of one year 
  ! (365/366) from a text file.
  !----------------------------------------------------------------   
  implicit none

contains

  function read1year_daily( filename, ndayyear ) result ( dval )
    !////////////////////////////////////////////////////////////////
    ! Function reads a file that contains 365 lines, each line for
    ! a daily value. 
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: ndayyear   ! 366 in leap years, 365 otherwise

    ! function return value
    real, allocatable, dimension(:) :: dval    ! Daily value to be read in

    ! allocate lengt of vector
    allocate( dval(ndayyear) )

    open( 20, file='../data/'//filename, status='old', form='formatted', action='read', err=888 )
    read( 20, *) dval
    close( 20 )

    return

    600 format (F9.7)
    888 write(0,*) 'READ1YEAR_DAILY: error opening file '//trim(filename)//'. Abort. '
    stop

  end function read1year_daily

end module daily_input


program splash
  !////////////////////////////////////////////////////////////////
  ! SPLASH MAIN PROGRAM 
  ! This subroutine maintains daily, monthly and annual quantities of 
  ! radiation, evapotranspiration, and soil moisture based on the 
  ! SPLASH methods.
  !----------------------------------------------------------------
  use modelparams
  use outputvars
  use daily_input

  implicit none

  ! local variables
  real :: my_lon
  real :: my_lat
  real :: my_elv

  integer :: yr        ! year AD (CE)
  integer :: inlen     ! =12 if monthly input data, =ndayyear if daily input data
  integer :: ndayyear  ! 366 in leap years, 365 otherwise

  ! state variables of SPLASH (daily updated, used by different SRs)
  real :: ra_d         ! daily extraterrestrial solar radiation (J/m^2)
  real :: rn_d         ! daytime net radiation, J/m^2

  ! Climate input variables
  real, allocatable, dimension(:) :: insf
  real, allocatable, dimension(:) :: intc
  real, allocatable, dimension(:) :: inppt

  ! Set current year (at which spinup is executed), and calculate number of days in year
  yr = 2000 ! year 2000 data is read in => leap year!
  ndayyear = julian_day(yr+1, 1, 1) - julian_day(yr, 1, 1)

  if (do_consistency_check) then

    inlen = ndayyear

    ! allocate size of array
    allocate( insf(ndayyear) )
    allocate( intc(ndayyear) )
    allocate( inppt(ndayyear) )

    ! Do consistency check with fixed daily temperature, sunshine fraction and precipitation
    print*, 'consistency check with fixed daily temperature, sunshine fraction and precipitation ...'
    insf(:) =  1.0
    intc(:) = 23.0
    inppt(:) = 0.9
  
  else if (use_daily_input) then

    inlen = ndayyear

    ! allocate size of array
    allocate( insf(ndayyear) )
    allocate( intc(ndayyear) )
    allocate( inppt(ndayyear) )

    ! Reading daily input data from file
    print*, 'reading daily climate input from files ...'

    print*, '   ../data/daily_sf_2000_cruts.txt'
    insf(:)  = read1year_daily( "../data/daily_sf_2000_cruts.txt", ndayyear )
    
    print*, '   ../data/daily_tair_2000_wfdei.txt'
    intc(:)  = read1year_daily( "../data/daily_tair_2000_wfdei.txt", ndayyear )
    
    print*, '   ../data/daily_pn_2000_wfdei.txt'
    inppt(:) = read1year_daily( "../data/daily_pn_2000_wfdei.txt", ndayyear )

    print*, '... done.'

  else

    inlen = nmonth

    ! allocate size of array
    allocate( insf(nmonth) )
    allocate( intc(nmonth) )
    allocate( inppt(nmonth) )
  
    ! Example data (Met Office average climate for this gridcell)
    insf  = (/0.21, 0.27, 0.30, 0.40, 0.39, 0.39, 0.40, 0.43, 0.36, 0.32, 0.23, 0.19/)
    intc  = (/4.80, 4.85, 7.10, 9.10, 12.4, 15.3, 17.6, 17.3, 14.6, 11.2, 7.55, 5.05/)
    inppt = (/61.0, 41.2, 44.5, 48.0, 46.4, 44.6, 46.0, 52.3, 50.3, 71.8, 66.3, 62.9/)

  end if

  ! Define input variables
  my_lon = -0.00 ! longitude, degrees
  my_lat = 37.7   ! latitude, degrees
  my_elv = 142.0  ! elevation, m

  print*,'point scale simulation for ...'
  print*,'    longitude =', my_lon
  print*,'    latitude  =', my_lat
  print*,'    elevation =', my_elv

  print*,'Running SPLASH ...'

  ! sanity check
  if (my_lon>180.0 .or. my_lon<(-180.0)) then
    print*,"Longitude outside range of validity (-180 to 180)!"
    stop
  endif

  if (my_lat>90.0 .or. my_lat<(-90.0)) then
    print*,"Latitude outside range of validity (-90 to 90)!"
    stop
  endif

  ! Initialize daily totals (output variables)
  print*,'initilasing daily variables ...'
  call initdaily
  print*,'... done'

  ! Initialize monthly totals (output variables)
  print*,'initialising monthly variables ...'
  call initmonthly
  print*,'... done'

  ! soil moisture spinup
  ! output variables are over-written in each year of spinup
  ! hence, values at the end of the spinup represent equilibrium
  print*,'spinning up soil moisture ...'
  call spin_up( yr, my_lon, my_lat, my_elv, inppt, intc, insf, inlen, ndayyear )
  print*,'... done'

  print*,'SPLASH sucessfully completed.'

contains

  subroutine spin_up( yr, lon, lat, elv, ppt, tc, sf, inlen, ndayyear )
    !----------------------------------------------------------------   
    ! Spins up the daily soil moisture
    !----------------------------------------------------------------   
    ! arguments
    integer, intent(in)                :: yr    ! year AD  
    real, intent(in)                   :: lon   ! longitude (degrees)
    real, intent(in)                   :: lat   ! latitude (degrees)
    real, intent(in)                   :: elv   ! altitude (m)
    real, dimension(inlen), intent(in) :: ppt   ! monthly precip (mm) 
    real, dimension(inlen), intent(in) :: tc    ! mean monthly temperature (deg C)
    real, dimension(inlen), intent(in) :: sf    ! mean monthly sunshine fraction (unitless)
    integer, intent(in)                :: inlen ! =12 if monthly input data, =ndayyear if daily input data
    integer, intent(in)                :: ndayyear ! number of days in this year

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
      start_sm = dwn(365)
      call run_one_year( lon, lat, elv, yr, ppt, tc, sf, inlen, ndayyear )
      end_sm   = dwn(365)

      ! Check to see if 1 Jan soil moisture matches 31 Dec:
      diff_sm  = abs(end_sm - start_sm)

      print*,'soil moisture in equil. when diff < 0.0001'
      print*,'diff ',diff_sm

      ! Increase counter
      spin_count = spin_count + 1

    enddo

    ! Once in equilibrium, write daily and monthly values of one year to file
    call write_to_file

  end subroutine spin_up


  subroutine write_to_file( )
    !----------------------------------------------------------------   
    ! Writes daily and monthly values to files
    !----------------------------------------------------------------   
    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR DAILY OUTPUT
    !----------------------------------------------------------------
    prefix = "./output/"

    !----------------------------------------------------------------
    ! WRTIING DAILY TOTALS
    !----------------------------------------------------------------
    print*,"writing daily totals to files"

    ! HO: daily solar irradiation, J/m2
    filnam=trim(prefix)//'ho.d.out'
    open(101,file=filnam,err=888,status='unknown')
    write(101,999) dho

    ! HN: daily net radiation, J/m2
    filnam=trim(prefix)//'hn.d.out'
    open(102,file=filnam,err=888,status='unknown')
    write(102,999) dhn

    ! QN: daily PPFD, mol/m2
    filnam=trim(prefix)//'qn.d.out'
    open(103,file=filnam,err=888,status='unknown')
    write(103,999) dqn

    ! CN: daily condensation water, mm
    filnam=trim(prefix)//'cn.d.out'
    open(104,file=filnam,err=888,status='unknown')
    write(104,999) dcn

    ! WN: daily soil moisture, mm
    filnam=trim(prefix)//'wn.d.out'
    open(105,file=filnam,err=888,status='unknown')
    write(105,999) dwn

    ! PN: daily precipitation, mm
    filnam=trim(prefix)//'pn.d.out'
    open(106,file=filnam,err=888,status='unknown')
    write(106,999) dpn

    ! RO: daily runoff, mm
    filnam=trim(prefix)//'ro.d.out'
    open(107,file=filnam,err=888,status='unknown')
    write(107,999) dro

    ! EQ_N: daily equilibrium ET, mm
    filnam=trim(prefix)//'eq_n.d.out'
    open(108,file=filnam,err=888,status='unknown')
    write(108,999) deq_n

    ! EP_N: daily potential ET, mm
    filnam=trim(prefix)//'ep_n.d.out'
    open(109,file=filnam,err=888,status='unknown')
    write(109,999) dep_n

    ! EA_N: daily actual ET, mm
    filnam=trim(prefix)//'ea_n.d.out'
    open(110,file=filnam,err=888,status='unknown')
    write(110,999) dea_n


    !----------------------------------------------------------------
    ! WRTIING MONTHLY TOTALS
    !----------------------------------------------------------------
    print*,"writing monthly totals to files"

    ! eq_m
    filnam=trim(prefix)//'eq_m.m.out'
    open(111,file=filnam,err=888,status='unknown')
    write(111,999) meq_m

    ! ep_m
    filnam=trim(prefix)//'ep_m.m.out'
    open(112,file=filnam,err=888,status='unknown')
    write(112,999) mep_m

    ! ea_m
    filnam=trim(prefix)//'ea_m.m.out'
    open(113,file=filnam,err=888,status='unknown')
    write(113,999) mea_m

    ! cpa
    filnam=trim(prefix)//'cpa.m.out'
    open(114,file=filnam,err=888,status='unknown')
    write(114,999) mcpa

    ! cwd
    filnam=trim(prefix)//'cwd.m.out'
    open(115,file=filnam,err=888,status='unknown')
    write(115,999) mcwd

    ! qm
    filnam=trim(prefix)//'qm.m.out'
    open(116,file=filnam,err=888,status='unknown')
    write(116,999) mqm

    return

  888    stop 'WRITE_TO_FILE: error opening output file'
  999    format (F20.8)

  end subroutine write_to_file


  subroutine run_one_year( lon, lat, elv, yr, ppt, tc, sf, inlen, ndayyear )
    !----------------------------------------------------------------   
    ! Calculates daily and monthly quantities for one year
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in)                   :: lon   ! longitude (degrees)
    real, intent(in)                   :: lat   ! latitude (degrees)
    real, intent(in)                   :: elv   ! altitude (m)
    integer, intent(in)                :: yr    ! year AD
    real, dimension(inlen), intent(in) :: ppt   ! monthly precip (mm) 
    real, dimension(inlen), intent(in) :: tc    ! mean monthly temperature (deg C)
    real, dimension(inlen), intent(in) :: sf    ! mean monthly sunshine fraction (unitless)

    integer, intent(in) :: inlen                         ! =12 if monthly input data, =ndayyear if daily input data
    integer, intent(in) :: ndayyear                      ! number of days in this year

    ! local variables
    real :: use_tc                       ! mean monthly or daily temperature (deg C)
    real :: use_sf                       ! mean monthly or daily sunshine fraction (unitless)

    real :: ppfd_d                       ! daily PPFD (mol/m^2)
    real :: eet_d                        ! daily EET (mm)
    real :: pet_d                        ! daily PET (mm)
    real :: aet_d                        ! daily AET (mm)
    real :: wc                           ! daily condensation (mm)
    real :: sw                           ! evaporative supply rate (mm/h)

    integer :: doy                       ! day of year
    integer :: moy                       ! month of year
    integer :: ndaymonth                 ! number of days in month (formerly 'nm')
    integer :: idx                       ! day of year corresponding to yesterday
    integer :: dm                        ! day of month

    ! Reset monthly totals
    call initmonthly

    ! Iterate through months
    do moy=1,nmonth

      ndaymonth = julian_day(yr, moy+1, 1) - julian_day(yr, moy, 1)

      ! Iterate through days in this month
      do dm=1,ndaymonth

        ! Calculate the day of the year
        ! xxx is this '+1' really necessary here?
        doy = julian_day(yr, moy, dm) - julian_day(yr, 1, 1) + 1

        ! Select daily or monthly input data
        if (use_daily_input) then
          use_tc   = tc(doy)
          use_sf   = sf(doy)
          dpn(doy) = ppt(doy)
        else
          use_tc   = tc(moy)
          use_sf   = sf(moy)
          dpn(doy) = ppt(moy)/ndaymonth
        end if

        ! Get index for yesterday (Note: zero indexing xxx)
        idx = int(doy-1)
        if (idx==0) idx = int(ndayyear)

        ! Calculate evaporative supply rate, mm/h
        sw = kCw * dwn(idx) / kWm

        ! Calculate radiation and evaporation quantities
        !print*,'calling evap for doy', doy, '...'
        call evap( ndayyear, lon, lat, doy, ppfd_d, eet_d, pet_d, aet_d, wc, elv=elv, yr=yr, sf=use_sf, tc=use_tc, sw=sw )
        !print*,'... done'

        ! Update soil moisture
        dwn(doy) = dwn(idx) + dpn(doy) + wc - aet_d

        if (dwn(doy)>kWm) then
          ! Bucket is full 
          ! * set soil moisture to capacity
          ! * add remaining water to monthly runoff total
          dro(doy) = dwn(doy)
          dro(doy) = dro(doy) - kWm
          dwn(doy) = kWm

        elseif (dwn(doy)<0) then
          ! Bucket is empty
          ! * set soil moisture to zero
          dwn(doy) = 0.0
          dro(doy) = 0.0

        else
          ! No runoff occurrs
          dro(doy) = 0.0

        endif

        ! Save the daily totals:
        dho(doy) = ra_d
        dhn(doy) = rn_d
        dqn(doy) = ppfd_d
        dcn(doy) = wc
        deq_n(doy) = eet_d
        dep_n(doy) = pet_d
        dea_n(doy) = aet_d

        ! Update monthly totals:
        meq_m(moy) = meq_m(moy) + eet_d
        mep_m(moy) = mep_m(moy) + pet_d
        mea_m(moy) = mea_m(moy) + aet_d
        mqm(moy)   =   mqm(moy) + ppfd_d

      enddo

      ! Calculate other monthly totals:
      mcpa(moy) = mea_m(moy)
      mcpa(moy) = mcpa(moy) / mep_m(moy)
      mcwd(moy) = mep_m(moy)
      mcwd(moy) = mcwd(moy) - mea_m(moy)

    enddo

  end subroutine run_one_year


  subroutine evap( ndayyear, lon, lat, doy, ppfd_d, eet_d, pet_d, aet_d, wc, elv, yr, sf, tc, sw )
    !----------------------------------------------------------------   
    ! This subroutine calculates daily radiation and evapotranspiration
    ! quantities
    ! - daily PPFD (ppfd_d), mol/m^2
    ! - daily EET (eet_d), mm
    ! - daily PET (pet_d), mm
    ! - daily AET (aet_d), mm
    ! - daily condensation (wc), mm
    !-------------------------------------------------------------  
    ! arguments

    ! inputs
    integer, intent(in) :: ndayyear      ! number of days in this year
    real, intent(in) :: lon              ! longitude, degrees
    real, intent(in) :: lat              ! latitude, degrees
    integer, intent(in) :: doy           ! day of the year (formerly 'n')

    ! outputs
    real, intent(out) :: ppfd_d          ! daily PPFD (ppfd_d), mol/m^2
    real, intent(out) :: eet_d           ! daily EET (eet_d), mm
    real, intent(out) :: pet_d           ! daily PET (pet_d), mm
    real, intent(out) :: aet_d           ! daily AET (aet_d), mm
    real, intent(out) :: wc              ! daily condensation (wc), mm

    ! optional inputs
    real,    intent(in), optional :: elv ! elevation, metres
    integer, intent(in), optional :: yr  ! year
    real,    intent(in), optional :: sf  ! fraction of sunshine hours 
    real,    intent(in), optional :: tc  ! mean daily air temperature, C
    real,    intent(in), optional :: sw  ! evaporative supply rate, mm/hr


    ! local variables
    real     :: my_elv                   ! elevation, metres; aux. var. for local use        
    integer  :: my_yr                    ! year; aux. var. for local use
    real     :: my_sf                    ! fraction of sunshine hours ; for local use
    real     :: my_tc                    ! mean daily air temperature, C; for local use
    real     :: my_sw                    ! evaporative supply rate, mm/hr; for local use

    real, dimension(2) :: out_berger     ! temporary var to store output of berger_tls
    real :: my_nu                        ! heliocentric longitudes
    real :: my_lambda                    ! heliocentric longitudes
    real :: my_rho
    real :: dr                           ! distance factor
    real :: delta                        ! declination angle 
    real :: ru                           ! variable substitute for u
    real :: rv                           ! variable substitute for v
    real :: hs                           ! sunset hour angle
    real :: tau_o, tau                   ! transmittivity (unitless)
    real :: rnl                          ! net longwave radiation (W/m^2)
    real :: rnn_d                        ! nighttime net radiation (J/m^2)
    real :: rw                           ! variable substitute (W/m^2)
    real :: hn                           ! net radiation cross-over hour angle
    real :: s                            ! slope of saturation vap press temp curve, Pa/K
    real :: pw                           ! density of water, kg/m^3
    real :: lv                           ! enthalpy of vaporization, J/kg
    real :: g                            ! psychrometric constant, Pa/K
    real :: econ                         ! Eq. 58, Documentation
    real :: rx                           ! variable substitute (mm/hr)/(W/m^2)
    real :: hi, cos_hi                   ! intersection hour angle, degrees


    ! Set default values for optional arguments
    if (.not.present(elv)) then
      my_elv = 0.0
    else
      my_elv = elv
    endif

    if (.not.present(yr)) then
      my_yr = 0
    else
      my_yr = yr
    endif

    if (.not.present(sf)) then
      my_sf = 1.0
    else
      my_sf = sf
    endif

    if (.not.present(tc)) then 
      my_tc = 23.0
    else
      my_tc = tc
    endif

    if (.not.present(sw)) then 
      my_sw = 1.0
    else
      my_sw = sw
    endif

    ! Error handle and assign required public variables:
    if (lon>180.0 .or. lon<(-180.0)) then
      print*,"Longitude outside range of validity (-180 to 180)!"
      stop
    endif

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
    ! Nothing done here anymore, because this is calculated 
    ! at higher level and passed onto here as argument.

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger (1978)
    out_berger = berger_tls( doy )
    
    ! consistency check
    if (doy==172) print*,'out_berger', out_berger

    ! xxx my_nu is not perfectly identical to Python version
    my_nu      = out_berger(1)
    my_lambda  = out_berger(2)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3. Calculate distance factor (dr), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger et al. (1993)
    my_rho = (1.0 - ke**2)/(1.0 + ke * dgcos( my_nu ))

    ! consistency check
    if (doy==172) print*,'my_rho', my_rho

    dr = (1.0/my_rho)**2

    ! consistency check
    if (doy==172) print*,'dr', dr

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 4. Calculate declination angle (delta), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Woolf (1968)
    delta = dasin( dgsin( my_lambda ) * dgsin( keps ) )   ! xxx arcsin is asin in Fortran?
    delta = degrees( delta )

    ! consistency check
    if (doy==172) print*,'delta', delta

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 5. Calculate variable substitutes (u and v), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ru = dgsin(delta) * dgsin(lat)
    rv = dgcos(delta) * dgcos(lat)

    ! consistency check
    if (doy==172) print*,'ru, rv ',ru,rv

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
      hs = dacos(-1.0*ru/rv)
      hs = degrees(hs)
    endif

    ! consistency check
    if (doy==172) print*,'hs',hs

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 7. Calculate daily extraterrestrial solar radiation / irradiation (ra_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 1.10.3, Duffy & Beckman (1993)
    ra_d = (86400.0/pi)*kGsc*dr*(ru*radians(hs) + rv * dgsin(hs))

    ! consistency check
    if (doy==172) print*,'daily extraterrestrial solar radiation / irradiation (ra_d), J/m^2', ra_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 8. Calculate transmittivity (tau), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Linacre (1968)
    tau_o = (kc + kd*my_sf)

    ! Eq. 2, Allen (1996)
    tau = tau_o*(1.0 + (2.67e-5)*my_elv)

    ! consistency check
    if (doy==172) print*,'tau_o ',tau_o
    if (doy==172) print*,'tau   ',tau

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 9. Calculate daily PPFD (ppfd_d), mol/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 57, Documentation
    ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra_d

    ! consistency check
    if (doy==172) print*,'ppfd_d',ppfd_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 10. Estimate net longwave radiation (rnl), W/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
    rnl = (kb + (1.0 - kb)*my_sf)*(kA - my_tc)

    ! consistency check
    if (doy==172) print*,'rnl',rnl

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 11. Calculate variable substitute (rw), W/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw = (1.0-kalb_sw)*tau*kGsc*dr

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
      hn = degrees( dacos((rnl - rw*ru)/(rw*rv)) )
    endif

    ! consistency check
    if (doy==172) print*,'cross-over hour angle, hn ',hn

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 13. Calculate daytime net radiation (rn_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 53, Documentation
    rn_d = (86400.0/pi) * (hn*(pi/180.0)*(rw*ru - rnl) + rw*rv*dgsin(hn))

    ! consistency check
    if (doy==172) print*,'daytime net radiation, rn_d',rn_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 14. Calculate nighttime net radiation (rnn_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 56, Documentation
    rnn_d = (86400.0/pi)*(radians(rw*ru*(hs-hn)) + rw*rv*(dgsin(hs)-dgsin(hn)) + rnl*(pi - 2.0*radians(hs) + radians(hn)))

    ! consistency check
    if (doy==172) print*,'nighttime net radiation, rnn_d ',rnn_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 15. Calculate water-to-energy conversion (econ), m^3/J
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Slope of saturation vap press temp curve, Pa/K
    s = sat_slope(my_tc)

    ! consistency check
    if (doy==172) print*,'Slope of saturation vap press temp curve, Pa/K, hn ', s

    ! Enthalpy of vaporization, J/kg
    lv = enthalpy_vap(my_tc)

    ! consistency check
    if (doy==172) print*,'Enthalpy of vaporization, J/kg, lv ', lv

    ! Density of water, kg/m^3
    pw = density_h2o(my_tc, elv2pres(my_elv))

    ! consistency check
    if (doy==172) print*,'Density of water, kg/m^3, pw ', pw

    ! Psychrometric constant, Pa/K
    g = psychro(my_tc, elv2pres(my_elv))

    ! consistency check
    if (doy==172) print*,'Psychrometric constant, Pa/K ', g

    ! Eq. 58, Documentation
    econ = s/(lv*pw*(s + g))

    ! consistency check
    if (doy==172) print*,'Econ ',econ

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 16. Calculate daily condensation (wc), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 68, Documentation
    wc = 1000.0 * econ * abs(rnn_d)

    ! consistency check
    if (doy==172) print*,'daily condensation (mm) ',wc

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 17. Estimate daily EET (eet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 70, Documentation
    eet_d = 1000.0*econ*(rn_d)

    ! xxx debug
    !print*,'eet_d',eet_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 18. Estimate daily PET (pet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 72, Documentation
    pet_d = (1.0+kw)*eet_d

    ! consistency check
    if (doy==172) print*,'daily PET (mm) ',pet_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx = 1000.0*3600.0*(1.0+kw)*econ

    ! consistency check
    if (doy==172) print*,'variable substitute rx ',rx

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 20. Calculate the intersection hour angle (hi), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cos_hi = my_sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
    if (cos_hi >= 1.0) then
      ! Supply exceeds demand:
      hi = 0.0
    elseif (cos_hi <= -1.0) then
      ! Supply limits demand everywhere:
      hi = 180.0
    else
      hi = degrees(acos(cos_hi))
    endif

    ! consistency check
    if (doy==172) print*,'intersection hour angle ',hi


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 21. Estimate daily AET (aet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 81, Documentation
    aet_d = (24.0/pi)*(radians(my_sw*hi) + rx*rw*rv*(dgsin(hn) - dgsin(hi)) + radians((rx*rw*ru - rx*rnl)*(hn - hi)))

    ! consistency check
    if (doy==172) print*,'daily AET (mm) ', aet_d


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


  function dgcos( x )
    !----------------------------------------------------------------   
    ! Calculates the cosine of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    use modelparams, only: pi

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
    use modelparams, only: pi

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
    use modelparams, only: pi

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

    ! arguments
    integer, intent(in) :: day   ! day of the year

    ! local variables
    real :: anm, ranm, anv, ranv
    real :: dlamm                ! Mean longitude for day of year
    real :: my_nu
    real :: my_tls
    real :: xee, xec, xse        ! variable substitutes
    real :: xlam                 ! Mean longitude for vernal equinox
    real :: tmp1, tmp2, tmp3     ! variable substitutes

    ! function return value
    real, dimension(2) :: berger_tls ! return "tuple": (my_nu = true anomaly, degrees; my_tls = true longitude, degrees)

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

    ! xxx debug
    !print*,'xlam', xlam

    ! Mean longitude for day of year:
    dlamm = xlam + (day - 80.0)*(360.0/ndayyear)

    ! Mean anomaly:
    anm = dlamm - komega
    ranm = radians(anm)

    ! True anomaly:
    ranv = (ranm + (2.0*ke - xec/4.0)*dsin(ranm) + 5.0/4.0*xee*dsin(2.0*ranm) + 13.0/12.0*xec*dsin(3.0*ranm))
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


  function julian_day( yr, moy, dom )
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
    integer :: julian_day                   ! Julian Ephemeris Day


    ! xxx debug
    !print*,"IN: yr, mo, dm :", yr, moy, dom

    if (moy <= 2) then
      my_yr = yr - 1
      my_moy = moy + 12
    else
      my_yr = yr
      my_moy = moy
    endif

    a = int(real(my_yr)/real(100))
    b = 2 - a + int(real(a)/real(4))
    julian_day=int(real(int(365.25*real(my_yr+4716)))+real(int(30.6001*real(my_moy+1)))+real(dom)+real(b)-1524.5)  ! xxx use moy instead of (moy+1) without zero-indexing as in Python?

    ! xxx debug
    !print*,"OUT: julian_day:", julian_day

  end function julian_day


  function sat_slope( tc )
    !----------------------------------------------------------------   
    ! Calculates the slope of the sat pressure temp curve, Pa/K
    ! Ref:      Eq. 13, Allen et al. (1998)
    !----------------------------------------------------------------   

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

    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real ::  enthalpy_vap ! enthalpy of vaporization, J/kg

    enthalpy_vap = 1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2

  end function enthalpy_vap


  function elv2pres( alt )
    !----------------------------------------------------------------   
    ! Calculates atm. pressure for a given elevation
    ! Ref:      Allen et al. (1998)
    !----------------------------------------------------------------   
    use modelparams, only: kPo, kTo, kL, kMa, kG, kR

    ! arguments
    real, intent(in) :: alt ! elevation above sea level, m

    ! function return value
    real ::  elv2pres ! atm. pressure for a given elevation

    elv2pres = kPo*(1.0 - kL*alt/kTo)**(kG*kMa/(kR*kL))

  end function elv2pres


  function density_h2o( tc, press )
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
    ! if (doy==172) print*,'atmospheric pressure, bar', pbar

    density_h2o = 1000.0*po*(ko + ca*pbar + cb*pbar**2.0)/(ko + ca*pbar + cb*pbar**2.0 - pbar)

  end function density_h2o


  function psychro( tc, press )
    !----------------------------------------------------------------   
    ! Calculates the psychrometric constant for a given temperature and pressure
    ! Ref: Allen et al. (1998); Tsilingiris (2008) 
    !----------------------------------------------------------------   
    use modelparams, only: kMa, kMv

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
    lv = enthalpy_vap(tc)

    ! Calculate psychrometric constant, Pa/K
    ! Eq. 8, Allen et al. (1998)
    psychro = cp*kMa*press/(kMv*lv)

  end function psychro

end program
