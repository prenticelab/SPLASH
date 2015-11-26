setwd("~/EEC/Winter_Research_project/WFDEI")

s <- 2007           #start year
e <- 2014           #end year
sm <- 1             #start month
em <- 12             #end month

for (i in s:e){
  for (x in sm:em)
    download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/Tair_WFDEI/Tair_WFDEI_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                destfile = paste("Air_temp/temp",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}
for (i in s:e){
  for (x in sm:em)
  download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/Tair_daily_WFDEI/Tair_daily_WFDEI_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                destfile = paste("Daily/Air_temp_daily/temp_daily",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}
for (i in s:e){
  for (x in sm:em)
  download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/Wind_WFDEI/Wind_WFDEI_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                destfile = paste("Wind_speed/wind",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}  
for (i in s:e){
  for (x in sm:em)
  download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/LWdown_WFDEI/LWdown_WFDEI_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                destfile = paste("Net_rad/LW",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}  
for (i in s:e){
  for (x in sm:em)
  download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/Rainf_daily_WFDEI_CRU/Rainf_daily_WFDEI_CRU_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                destfile = paste("Daily/Precipitation/RF",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}
for (i in s:e){
  for (x in sm:em)
    download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/Snowf_daily_WFDEI_CRU/Snowf_daily_WFDEI_CRU_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                  destfile = paste("Daily/Precipitation/Snowf",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}
for (i in s:e){
  for (x in sm:em)
  download.file(paste("ftp://rfdata:forceDATA@ftp.iiasa.ac.at/WFDEI/SWdown_daily_WFDEI/SWdown_daily_WFDEI_",i,formatC(x,width = 2,flag = 0),".nc.gz",sep = ""), 
                destfile = paste("Daily/Shortwave_rad/SW",i,formatC(x,width = 2,flag = 0),".nc.gz", sep = ""), method = "libcurl")
}

## blah blah
## blah blah
## blah blah 
# blah
# blah
