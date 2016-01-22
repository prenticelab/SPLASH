# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# data.R
#
# last updated: 2016-01-22
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Robust indices of radiation
# evapo-transpiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions to handle the file IO for reading and writing
# data, i.e.:
#   read_csv(character fname, double y=-1)
#   read_txt(list my_data, character fname, character var, double y=-1)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - created read_csv and read_txt functions [15.02.23]
#
#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     read_csv
# Inputs:   - character, file name (fname)
#           - double, year (y)
# Returns:  list object (data)
#             $file_name ............ file name
#             $sf ................... sunshine fraction
#             $tair ................. air temperature
#             $pn ................... precipitation
#             $num_lines ............ number of data points
#             $year ................. year of data
# Features: Reads all three daily input variables (sf, tair, and pn)
#           for a single year from a CSV file that includes a header
# ************************************************************************
read_csv <- function(fname, y=-1) {
    my_data <- list()
    my_data$file_name <- fname

    # Read data from CSV file and save to return list:
    DATA <- read.csv(fname)
    my_data$sf <- DATA$sf
    my_data$tair <- DATA$tair
    my_data$pn <- DATA$pn
    my_data$num_lines <- dim(DATA)[1]

    my_year <- y
    if (y == -1) {
        if (dim(DATA)[1] == 366) {
            my_year <- 2000
        } else if (dim(DATA)[1] == 365) {
            my_year <- 2001
        }
    }
    my_data$year <- my_year
    return (my_data)
}


# ************************************************************************
# Name:     read_txt
# Inputs:   - list object (my_data)
#           - character, file name (fname)
#           - character, variable name (var)
#           - double, year (y)
# Returns:  list object (my_data)
# Features: Reads plain text file (no header) of one of the input
#           arrays
# ************************************************************************
read_txt <- function(my_data, fname, var, y=-1) {
    my_data$file_name <- c(my_data$file_name, fname)
    DATA <- scan(fname)
    if (var == "sf") {
        my_data$sf <- DATA
    } else if (var == "tair") {
        my_data$tair <- DATA
    } else if (var == "pn") {
        my_data$pn <- DATA
    }
    my_data$num_lines <- length(DATA)

    my_year <- y
    if (y == -1) {
        if (length(DATA) == 366) {
            my_year <- 2000
        } else if (length(DATA) == 365) {
            my_year <- 2001
        }
    }
    my_data$year <- my_year
    return (my_data)
}
