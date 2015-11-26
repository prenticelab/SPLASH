#Calculate mean temp
library(SDMTools)
library(ncdf)
library(compiler)
setwd("~/EEC/Winter_Research_project/WFDEI/2013")

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

nc2ascii_uncomp <- function(file_name){
  setwd(file_name)
  ctr = 0
  ctr1 = 0
  l <- list.files(pattern = "\\.nc$")  
  for (i in l){
    A<-open.ncdf(i)
    ctr = ctr + 1
    x<-length(A$dim$tstep$vals)
    if (ctr == 1)
      {
      n <- seq(from = x-8, to = x)
      for (j in n)
        {
        ctr1 = ctr1 + 1
        temp = get.var.ncdf( nc = A, varid = "Tair", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        temp <- as.asc(temp, xll=-179.75, yll = -89.75, cellsize = 0.5)    #convert matrix to ASCII
        write.asc(temp, file = paste("temp", formatC(ctr1,width = 4,flag = 0),".asc", sep = "")) 
      }
    }
    else
    {
      if (ctr == 14)
        {
        for (j in 1:8)
        {
        ctr1 = ctr1 + 1
        temp = get.var.ncdf( nc = A, varid = "Tair", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        temp <- as.asc(temp, xll=-179.75, yll = -89.75, cellsize = 0.5)    #convert matrix to ASCII
        write.asc(temp, file = paste("temp", formatC(ctr1,width = 4,flag = 0),".asc", sep = ""))
        }
          
      }
      else
      {
        for (j in 1:x){
        ctr1 = ctr1 + 1
        temp = get.var.ncdf( nc = A, varid = "Tair", start=c(1, 1, j), count = c(-1, -1, 1))  #read timestep j
        temp <- as.asc(temp, xll=-179.75, yll = -89.75, cellsize = 0.5)    #convert matrix to ASCII
        write.asc(temp, file = paste("temp", formatC(ctr1,width = 4,flag = 0),".asc", sep = ""))
        }
      }
    }
  }
}

nc2ascii <- cmpfun(nc2ascii_uncomp)

#eg
nc2ascii(file_name = "Air_temp/")            
nc2ascii(file_name = "Wind_speed/")

setwd("~/EEC/Winter_Research_project/WFDEI/2013")
system.time(nc2ascii(file_name = "Air_temp/"))

lat = min(seq(-89.75, 89.75, 0.5))
lon = min(seq(-179.75, 179.75, 0.5))
cellsize = 0.5

## we can then convert the matrices into ASCII and then into dataframe. The dataframe can then be named accordingly

l_file<-list.files( pattern = "\\.asc$")
tdata=asc2dataframe(l_file)


write.asc(temp, file = paste("temp", formatC(ctr1,width = 4,flag = 0),".asc", sep = ""))