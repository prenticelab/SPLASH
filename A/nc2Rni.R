# LW and SW to NetRad
library(SDMTools)
library(ncdf)
library(compiler)
setwd("~/EEC/Winter_Research_project/WFDEI/2013")

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: nc2NetRad
# *
# * Input: folder directory containing 14 months data (e.g.Air_temp/)
# *   
# *
# * Return: dataframe containing all variables of each timestep
# *
# *
# ************************************************************************

nc2NetRad <- function(file_name){
  setwd(file_name)
  ctr = 0
  ctr1 = 0
  lat = min(seq(-89.75, 89.75, 0.5))
  lon = min(seq(-179.75, 179.75, 0.5))
  cellsize = 0.5
  l <- list.files(pattern = "\\.nc$")
  LW <- grep(pattern = "LWdown", l, value = T)
  SW <- grep(pattern = "SWdown", l, value = T)
  Tair <- grep(pattern = "Tair", l, value = T)
  albedo <- 0.17 #standard SPLASH for broadband albedo
  emissivity <- 0.99 # assumed for leaves with Leaf Area Index > 2
  sigma <- 5.670367 #Stefan Boltzman constant, Standard uncertainty 0.000013
  for (i in LW){
    A <- open.ncdf(i)
    ctr = ctr + 1
    B <- open.ncdf(SW[ctr])
    C <- open.ncdf(Tair[ctr])
    x<-length(A$dim$tstep$vals)
    if (ctr == 1)
    {
      n <- seq(from = x-7, to = x)
      for (j in n)
      {
        ctr1 = ctr1 + 1
        RadA = get.var.ncdf( nc = A, varid = "LWdown", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        RadB = get.var.ncdf( nc = B, varid = "SWdown", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        T_air = get.var.ncdf( nc = C, varid = "Tair", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        Rni <- (RadA - albedo*RadA) + (RadB - emissivity*sigma*(T_air^4))
        Rni <- as.asc(Rni, xll=lon, yll = lat, cellsize = cellsize)    #convert matrix to ASCII
        write.asc(Rni, file = paste("Rni", formatC(ctr1,width = 4,flag = 0),".asc", sep = "")) 
      }
    }
    else
    {
      if (ctr == 14)
      {
        for (j in 1:8)
        {
          ctr1 = ctr1 + 1
          RadA = get.var.ncdf( nc = A, varid = "LWdown", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
          RadB = get.var.ncdf( nc = B, varid = "SWdown", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
          T_air = get.var.ncdf( nc = C, varid = "Tair", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
          Rni <- (RadA - albedo*RadA) + (RadB - emissivity*sigma*(T_air^4))
          Rni <- as.asc(Rni, xll=lon, yll = lat, cellsize = cellsize)    #convert matrix to ASCII
          write.asc(Rni, file = paste("Rni", formatC(ctr1,width = 4,flag = 0),".asc", sep = ""))
        }
        
      }
      else
      {
        for (j in 1:x){
          ctr1 = ctr1 + 1
          RadA = get.var.ncdf( nc = A, varid = "LWdown", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
          RadB = get.var.ncdf( nc = B, varid = "SWdown", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
          T_air = get.var.ncdf( nc = C, varid = "Tair", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
          Rni <- (RadA - albedo*RadA) + (RadB - emissivity*sigma*(T_air^4))
          Rni <- as.asc(Rni, xll=lon, yll = lat, cellsize = cellsize)    #convert matrix to ASCII
          write.asc(Rni, file = paste("Rni", formatC(ctr1,width = 4,flag = 0),".asc", sep = "")) 
        }
      }
    }
  }
}

# nc2NetRad_comp <- cmpfun(nc2NetRad) compile incase it is faster. Uncompiled is faster

#execution
system.time(nc2NetRad(file_name = "~/EEC/Winter_Research_project/WFDEI/2013/Rni"))







