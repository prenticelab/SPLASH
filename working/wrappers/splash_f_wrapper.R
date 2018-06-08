###############################################
##Interface for the Fortran version of SPLASH##
###############################################
dyn.load("../f90_version/libsplash.so")

#Initialise arrays that the Fortran code uses to save computations 
initialise_work_arrays <- function() {
    null <- .C('__splash_MOD_initdaily', PACKAGE='libsplash')
    null <- .C('__splash_MOD_initmonthly', PACKAGE='libsplash')
}

#Spin up SPLASH with the given forcings and parameters 
spin_up <- function(year, lat, elv, ppt, tc, sf, len) {
    null <- .C('__splash_MOD_spin_up_sm',
               as.integer(year),
               as.double(lat),
               as.double(elv),
               as.double(ppt),
               as.double(tc),
               as.double(sf),
               as.integer(len),
               as.double(0.0167),
               as.double(23.44),
               as.double(283.0),
               PACKAGE='libsplash')
}

#Run SPLASH for a particular day
run_one_day <- function( lat, elv, doy, yr, wn, pr, tc, sf ) {
    #TODO:Figure out how to set wn
    null <- .C('__splash_MOD_run_one_day',
               as.double(lat),
               as.double(elv),
               as.integer(doy),
               as.integer(yr),
               as.double(pr),
               as.double(tc),
               as.double(sf),
               PACKAGE='libsplash')
}

#Run SPLASH for a particular year
run_one_year <- function(year, lat, elv, ppt, tc, sf, len) {
    null <- .C('__splash_MOD_run_one_year',
               as.double(lat),
               as.double(elv),
               as.integer(year),
               as.double(ppt),
               as.double(tc),
               as.double(sf),
               as.integer(len),
               PACKAGE='libsplash')
}

#Write the held variables for SPLASH to predetermined files
write_to_file <- function() {
    null <- .C('__splash_MOD_write_to_file',
               PACKAGE='libsplash')
}

#Only run if not in interactive mode so that we can import the interface without running this test
if(!interactive()) {
    source("../r_version/data.R") #So we can import test data

    # Location constants:
    my_lat <- 37.7
    my_elv <- 142

    # Calculate days in the year
    y <- 2000
    #We can't provide an interface for the Fortran julian day interface as R only supports subroutines in Fortran, not functions
    ny <- 365 #julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)

    # Example data (San Francisco, 2000 CE)
    #   $file_name, file name
    #   $sf, fractional sunshine hours, unitless
    #   $tair, air temperature, deg C
    #   $pn, daily precipitation, mm
    #   $num_lines, number of lines of data
    #   $year, year (Gregorian calendar)
    my_file <- "../../data/example_data.csv"
    my_data <- read_csv(my_file, y)
    my_data$lat_deg <- my_lat
    my_data$elv_m <- my_elv

    # Spin up the soil moisture content
    null <- spin_up(my_data$year,
                    my_data$lat_deg,
                    my_data$elv_m,
                    my_data$pn,
                    my_data$tair,
                    my_data$sf,
                    ny)

    # Run SPLASH for a year with the test forcings
    null <- run_one_year(my_data$year,
                         my_data$lat_deg,
                         my_data$elv_m,
                         my_data$pn,
                         my_data$tair,
                         my_data$sf,
                         ny)

    #Save the outputs for the Fortran code 
    null <- write_to_file()
}
