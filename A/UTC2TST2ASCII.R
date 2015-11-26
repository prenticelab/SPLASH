# From formula in Excel available from 
# http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
# 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: UTC_TST1
# *
# * Input: d1900 = Julian Days from 01-01-1900 (use Excel) 
# *        hour = Time in hr in UTC
# *        longitude = Longitude 
# *        
# *
# * Return: matrix of Local Solar Time
# *
# ************************************************************************
library(circular)
library(SDMTools)
library(compiler)

UTC_TST1 <- function(d1900, hour, longitude, start_date)
{
  J_D <- d1900 + 2415018.5 + hour/24  #Julian day, hr in UTC
  G2 <- (J_D - 2451545)/36525
  I2 <- (280.46646 + G2*(36000.76983 + G2*0.0003032))%%360
  J2 <- 357.52911 + G2*(35999.05029 - 0.0001537*G2)
  K2 <- 0.016708634 - G2*(0.000042037+0.0000001267*G2)
  Q2 <- 23 + (26 + ((21.448 - G2*(46.815 + G2*(0.00059 - G2*0.001813))))/60)/60
  R2 <- Q2 + 0.00256*cos(rad(125.04 - 1934.136*G2))
  U2 <- tan(rad(R2/2))*tan(rad(R2/2))
  Eq_time <- 4*deg(U2*sin(2*rad(I2)) - 2*K2*sin(rad(J2)) + 
                     4*K2*U2*sin(rad(J2))*cos(2*rad(I2)) - 
                     0.5*U2*U2*sin(4*rad(I2)) - 1.25*K2*K2*sin(2*rad(J2)))
  day_of_year = d1900 - start_date         
  time_offset = Eq_time + 4*longitude 
  tst_hr = hour + time_offset/60               #tst in hours as in LST MODIS day_view_time
  tst_norm = (day_of_year-1)*24 + tst_hr       #calculate number of hours from first day of year
  return(tst_norm)
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: UTC_ASCI1
# *
# * Input: d1900 = list of numbers of Julian Days from 01-01-1900 (use Excel) 
# *        time_hour = list of time in hr in UTC
# *        lon_map = matrix of longitude 
# *                 (transformed to suit nc_files, and somehow opens properly in ArcGIS)
# *        lon = longitude of centre of lower left cell 
# *        lat = latitude of centre of lower left cell 
# *             
# *        
# *        
# *
# * Return: matrix of Local Solar Time
# *
# *
# * Depend: UTC_TST1        
# *
# *
# ************************************************************************
UTC2ASCI1 <- function(d1900, time_hour, lon_map, lon, lat)
{
  ctr = 0                  #set counter
  start_date <- d1900[1]
  for (i in d1900)
  {
    for (j in time_hour)
    {
      tst_norm<-UTC_TST1(d1900 = i, hour = j, longitude = lon_map, start_date = start_date)  
      ctr = ctr + 1
      tst <- as.asc(tst_norm, xll = lon, yll = lat, cellsize = 0.5)   #due to the extracted nc files having swapped dimensions, but opens fine in ArcGIS
      write.asc(tst, file = paste("Time_steps/","t", formatC(ctr, width = 4,flag = 0),".asc", sep = ""))
    }
  }
}


# Remarks: Compiled UTC2ASCI1 for increased speed 
UTC2ASCI_COMP1 <- cmpfun(UTC2ASCI1)

# Execution
# Set working directory
setwd("~/EEC/Winter_Research_project/WFDEI/2013")

# Julian days since 1900, 31-12-2012 to 01-01-2014
d1900 <- seq(from = 41274, to = 41640, by = 1)

# Time each time step represents
time_hour <- list(-1.5, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5)

# Define coordinate of midpoint of cell at bottom left corner
lat = min(seq(-89.75, 89.75, 0.5))
lon = min(seq(-179.75, 179.75, 0.5))

#generate matrix with longitude, but transformed as the matrix generated from nc-files were transformed, but opens fine in ArcGIS
lon_map <- matrix(data = rep(x = seq(from = -179.75, to = 179.75, by = 0.5), 
               each = 360), nrow = 720, ncol = 360, byrow = T) 

# Execute
system.time(UTC2ASCI_COMP1(d1900 = d1900, time_hour = time_hour, lon_map = lon_map, lon = lon, lat = lat))






