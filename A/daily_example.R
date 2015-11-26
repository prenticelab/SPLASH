# EXAMPLE
# assume calculating day 1 of feb
# using "daily" function to linearly interopolate monthly climatologies to daily 
# mly (monthly value) as inputs, and return dly (daily values)
daily<-function(mly,grid)
{
  tmly <- matrix(nrow=grid,ncol=14)
  dly <- matrix(nrow=grid,ncol=365)
  
  tmly[,1]<-mly[,12]        #A# 14 months to interpolate for the 12 months, 12 months plus months before and after the period (or assume value for Dec and Jan to be the same in different years)
  tmly[,14]<-mly[,1]
  tmly[,2]<-mly[,1]
  tmly[,3]<-mly[,2]
  tmly[,4]<-mly[,3]
  tmly[,5]<-mly[,4]
  tmly[,6]<-mly[,5]
  tmly[,7]<-mly[,6]
  tmly[,8]<-mly[,7]
  tmly[,9]<-mly[,8]
  tmly[,10]<-mly[,9]
  tmly[,11]<-mly[,10]
  tmly[,12]<-mly[,11]
  tmly[,13]<-mly[,12]
  
  mtd<-function(w,i,k)  # parameter w: 1- in the first half of month
    #              0- in the second half of month
    # parameter i: in which month
    # parameter k: in which day of that month
  {
    i<-i+1 # month loop controlor from 2 to 13
    
    days<-c(31,31,28,31,30,31,30,31,31,30,31,30,31,31)
    
    if(w==1)                      #first half of month
    {
      mbef<-tmly[,i-1]         #values of Jan
      maft<-tmly[,i]           #values of Feb
      dbef<-days[i-1]          #number of days in Jan
      daft<-days[i]            #number of days in Feb
      d15<-k+floor(dbef/2)     #the day of the month + number of days in Jan/2 = 1 + 15 = 16     
    }
    else if(w==0)                 #second half of month
    {
      mbef<-tmly[,i]           #values of Feb
      maft<-tmly[,i+1]         #values of March
      dbef<-days[i]            #values of days in Feb
      daft<-days[i+1]          #values of days in March
      d15<-k-floor(daft/2)     #the day of the month + number of days in March/2 = 1-15 = -14
    }
    
    inc<-(maft-mbef)/((dbef+daft)/2)      #calculate increase per day
    
    dly<-mbef+(inc%*%t(d15))              #daily values equals the values of month before + (increase*number of days from middle of previous month)
  }
  
  dly[,1:15]<-mtd(1,1,c(1:15))        # January
  dly[,16:31]<-mtd(0,1,c(16:31))
  dly[,32:44]<-mtd(1,2,c(1:13))       # Feberary
  dly[,45:59]<-mtd(0,2,c(14:28))
  dly[,60:74]<-mtd(1,3,c(1:15))       # March
  dly[,75:90]<-mtd(0,3,c(16:31))
  dly[,91:104]<-mtd(1,4,c(1:14))      # April
  dly[,105:120]<-mtd(0,4,c(15:30))
  dly[,121:135]<-mtd(1,5,c(1:15))     # May
  dly[,136:151]<-mtd(0,5,c(16:31))
  dly[,152:165]<-mtd(1,6,c(1:14))     # June
  dly[,166:181]<-mtd(0,6,c(15:30))    
  dly[,182:196]<-mtd(1,7,c(1:15))     # July
  dly[,197:212]<-mtd(0,7,c(16:31))
  dly[,213:227]<-mtd(1,8,c(1:15))     # August
  dly[,228:243]<-mtd(0,8,c(16:31))   
  dly[,244:257]<-mtd(1,9,c(1:14))     # September
  dly[,258:273]<-mtd(0,9,c(15:30))
  dly[,274:288]<-mtd(1,10,c(1:15))    # October
  dly[,289:304]<-mtd(0,10,c(16:31))
  dly[,305:318]<-mtd(1,11,c(1:14))    # November
  dly[,319:334]<-mtd(0,11,c(15:30))
  dly[,335:349]<-mtd(1,12,c(1:15))    # December  
  dly[,350:365]<-mtd(0,12,c(16:31))
  
  return(dly)
}