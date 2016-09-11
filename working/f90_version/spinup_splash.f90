program test_spinup_splash
  !////////////////////////////////////////////////////////////////
  ! EXECUTE TEST 4 FOR SPLASH (spinup with daily input)
  ! see https://bitbucket.org/labprentice/splash/wiki/Home
  !
  ! LAST UPDATED: 2016-09-11
  !
  ! Copyright (C) 2016 Prentice Lab
  !
  ! This file is part of the SPLASH model.
  !
  ! SPLASH is free software: you can redistribute it and/or modify it under
  ! the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 2.1 of the License, or
  ! (at your option) any later version.
  !
  ! SPLASH is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU Lesser General Public License for more details.
  !
  ! You should have received a copy of the GNU Lesser General Public License
  ! along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
  !----------------------------------------------------------------
  use splash

  implicit none

  ! Climate input variables
  real, allocatable, dimension(:) :: insf
  real, allocatable, dimension(:) :: intc
  real, allocatable, dimension(:) :: inppt

  integer :: use_yr   ! year for which soil is to be spun up
  integer :: inlen    ! length of input vectors

  ! Chose whether to use monthly or daily input data
  use_daily_input = .true.

  ! initialise water balance
  waterbal%sm = 0.0
  waterbal%ro = 0.0

  use_yr = 2000

  ! number of days in this year
  inlen = get_julian_day( use_yr+1, 1, 1 ) - get_julian_day( use_yr, 1, 1 )

  ! Test 4
  print*,'Executing TEST 4 ...'

  ! get daily or monthly climate inputs (temp, prec, sunshine fraction)
  if (use_daily_input) then
    print*,'getting daily input variables ...'
    call get_input_daily1year()
  else
    print*,'getting monthly input variables ...'
    call get_input_monthly1year()
  end if

  ! Initialize daily totals (output variables)
  print*,'initilasing daily output variables ...'
  call initdaily()
  print*,'... done'

  ! Initialize monthly totals (output variables)
  print*,'initialising monthly output variables ...'
  call initmonthly()
  print*,'... done'

  ! don't print various splash variables from whithin functions
  verbose = .false.

  print*,'spinning up soil moisture ...'
  call spin_up_sm( yr=use_yr, lat=37.7, elv=142.0, ppt=inppt, tc=intc, sf=insf, inlen=inlen )
  print*,'... done'

  ! Once in equilibrium, write daily and monthly values of one year to file
  print*,'writing output ...'
  call write_to_file()

  print*, '--------------- '
  print*, 'TEST 4 results (daily soil moisture) are in output/*.d.out '
  print*, '--------------- '

  print*,'SPLASH sucessfully completed.'


contains


  subroutine get_input_monthly1year()
    !----------------------------------------------------------------
    ! Populates monthly climate input vectors with values (obsolete)
    !----------------------------------------------------------------
    inlen = nmonth

    ! allocate size of array
    allocate( insf(nmonth) )
    allocate( intc(nmonth) )
    allocate( inppt(nmonth) )

    ! Example data (Met Office average climate for this gridcell)
    insf  = (/0.21, 0.27, 0.30, 0.40, 0.39, 0.39, 0.40, 0.43, 0.36, 0.32, 0.23, 0.19/)
    intc  = (/4.80, 4.85, 7.10, 9.10, 12.4, 15.3, 17.6, 17.3, 14.6, 11.2, 7.55, 5.05/)
    inppt = (/61.0, 41.2, 44.5, 48.0, 46.4, 44.6, 46.0, 52.3, 50.3, 71.8, 66.3, 62.9/)

  end subroutine get_input_monthly1year


  subroutine get_input_daily1year()
    !----------------------------------------------------------------
    ! Reads daily climate input from files
    !----------------------------------------------------------------
    ! allocate size of array

    inlen = 366

    allocate( insf(inlen) )
    allocate( intc(inlen) )
    allocate( inppt(inlen) )

    ! Reading daily input data from file
    print*, 'reading daily climate input from files ...'

    print*, '   daily_sf_2000_cruts.txt'
    insf(:)  = read1year_daily( "daily_sf_2000_cruts.txt", inlen )

    print*, '   daily_tair_2000_wfdei.txt'
    intc(:)  = read1year_daily( "daily_tair_2000_wfdei.txt", inlen )

    print*, '   daily_pn_2000_wfdei.txt'
    inppt(:) = read1year_daily( "daily_pn_2000_wfdei.txt", inlen )

    print*, '... done.'

  end subroutine get_input_daily1year


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

    ! allocate length of vector
    allocate( dval(ndayyear) )

    open( 20, file='../../data/'//filename, status='old', form='formatted', action='read', err=888 )
    read( 20, *) dval
    close( 20 )

    return

    600 format (F9.7)
    888 write(0,*) 'READ1YEAR_DAILY: error opening file '//trim(filename)//'. Abort. '
    stop

  end function read1year_daily


end program test_spinup_splash
