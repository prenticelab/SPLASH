#Calculate mean temp
library(SDMTools)
library(ncdf)
library(compiler)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: nc2asc
# *
# * Input: folder directory containing 14 months data (e.g.Air_temp/)
# *   
# *
# * Return: dataframe containing all variables of each timestep
# *
# *
# ************************************************************************

nc2ascii <- function(file_name, var_name){
  setwd(file_name)
  NODATA = 1e+19
  ctr = 0
  ctr1 = 0
  lat = min(seq(-89.75, 89.75, 0.5))
  lon = min(seq(-179.75, 179.75, 0.5))
  cellsize = 0.5
  l <- list.files(pattern = "\\.nc$")  
  for (i in l){
    A<-open.ncdf(i)
    ctr = ctr + 1
    x<-length(A$dim$tstep$vals)
    if (ctr == 1)
      {
      n <- seq(from = x-7, to = x)
      for (j in n)
        {
        ctr1 = ctr1 + 1
        temp = get.var.ncdf( nc = A, varid = var_name, start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        temp[temp>NODATA] <- NA
        temp <- as.asc(temp, xll=lon, yll = lat, cellsize = cellsize)    #convert matrix to ASCII
        write.asc(temp, file = paste(var_name, formatC(ctr1,width = 4,flag = 0),".asc", sep = "")) 
      }
    }
    else
    {
      if (ctr == 14)
        {
        for (j in 1:8)
        {
        ctr1 = ctr1 + 1
        temp = get.var.ncdf( nc = A, varid = var_name, start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        temp[temp>NODATA] <- NA
        temp <- as.asc(temp, xll=lon, yll = lat, cellsize = cellsize)    #convert matrix to ASCII
        write.asc(temp, file = paste(var_name, formatC(ctr1,width = 4,flag = 0),".asc", sep = ""))
        }
          
      }
      else
      {
        for (j in 1:x){
        ctr1 = ctr1 + 1
        temp = get.var.ncdf( nc = A, varid = var_name, start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        temp[temp>NODATA] <- NA
        temp <- as.asc(temp, xll=lon, yll = lat, cellsize = cellsize)    #convert matrix to ASCII
        write.asc(temp, file = paste(var_name, formatC(ctr1, width = 4,flag = 0),".asc", sep = ""))
        }
      }
    }
  }
}

# nc2ascii_comp <- cmpfun(nc2ascii) compile incase it is faster. Uncompiled is faster in this case

#execution
setwd("~/EEC/Winter_Research_project/WFDEI/2013")
system.time(nc2ascii(file_name = "~/EEC/Winter_Research_project/WFDEI/2013/Air_temp/", var_name = "Tair"))
system.time(nc2ascii(file_name = "~/EEC/Winter_Research_project/WFDEI/2013/Wind_speed", var_name = "Wind"))



## test for speed
setwd("~/EEC/Winter_Research_project/WFDEI/2013")
system.time(nc2ascii(file_name = "Air_temp/"))
system.time(nc2ascii_uncomp(file_name = "Air_temp/"))







