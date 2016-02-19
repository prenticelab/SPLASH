program test_splash
  !////////////////////////////////////////////////////////////////
  ! EXECUTE TESTS 1, 2, AND 3 FOR SPLASH (without spinup)
  ! see https://bitbucket.org/labprentice/splash/wiki/Home
  ! Copyright (C) 2015, see LICENSE, Benjamin D. Stocker
  ! contact: benjamin.stocker@gmail.com
  !----------------------------------------------------------------
  use _splash

  implicit none

  ! initialise water balance
  waterbal%sm = 75.0
  waterbal%ro = 0.0

  ! Test 3 (see https://bitbucket.org/labprentice/splash/wiki/Home)
  print*,'Running evap(), TEST 2 ...'

  ! print various splash variables from whithin functions
  verbose = .true.

  print*, 'TEST 1 and 2 results: '
  print*, '--------------- '
  out_evap = evap( lat=37.7, doy=172, elv=142.0, yr=2000, sf=1.0, tc=23.0, sw=0.9 )

  print*, '--------------- '
  print*, 'daily TOA solar irradiation (J/m2) ', out_evap%ra
  print*, 'daily net radiation (J/m2)         ', out_evap%rn
  print*, 'daily PPFD (mol/m^2)               ', out_evap%ppfd
  print*, 'daily EET (mm)                     ', out_evap%eet
  print*, 'daily PET (mm)                     ', out_evap%pet
  print*, 'daily AET (mm)                     ', out_evap%aet
  print*, 'daily condensation (mm)            ', out_evap%cn
  print*, '--------------- '


  ! Test 3 (see https://bitbucket.org/labprentice/splash/wiki/Home)
  print*,'Running run_one_day(), TEST 3 ...'

  ! don't print various splash variables from whithin functions
  verbose = .false.

  call run_one_day( lat=37.7, elv=142.0, doy=172, yr=2000, pr=5.0, tc=23.0, sf=1.0 )

  print*, 'TEST 3 results: '
  print*, '--------------- '
  print*, 'soil moisture (mm)                 ', waterbal%sm
  print*, 'runoff (mm)                        ', waterbal%ro
  print*, '--------------- '


end program test_splash
