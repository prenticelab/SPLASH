
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: interpolate
# *
# * Input: dataframe of columns of time, temp, and modis view time (A)
# *        A is generated generated from UTC_TST and nc2dataframe
# *
# * Return: list, meteorological data
# *
# * Features: Interpolates and calculate monthly means
# *
# * Depends: UTC_TST, nc2dataframe
# ************************************************************************

interpolate <- function(B)
  {
  A <- read.csv("Air_temp/tempx.csv")  #dataframe generated from nc2dataframe
  A1 <- read.csv("Air_temp/tx.csv")     #dataframe of solar time generated from UTC_TST from New_methods
  B #merged dataframe A and A1
  month_day <- data.frame(start = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 336), end = c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365))
  month_day$start[2]
  for (i in 1:12)
    {
    sequence<- seq(from = month_day$start[i], to = month_day$end[i])*24    
    assign(x = paste("seq",i, sep = ""), value = sequence)
    time2interpolate <- B[,2] + sequence #What ever column Modis time is in)
    ctrx = 0      #reset counter

    for (j in time2interpolate)
      {
      ctrx = ctrx + 1
      Wt<-B[1:368]    #1:368 represents number of columns (check again)
      t_b <- max(Wt[Wt<j] )  # time before
      t_a <- min(Wt[Wt>=j])  # time after
      Ta_b <- which(B == t_b, arr.ind = T)[2]  #columns before start of Tair
      Ta_a <- which(B == t_a, arr.ind = T)[2]  #columns in which air temp value should be taken from
      Tair <- B[368:700]    #whatever column Tair represents
      Tair_b <- Tair[Ta_b]  #Air temp value before j
      Tair_a <- Tair[Ta_b]  #Air temp value after j
      Tair_int <- Tair_b + (Tair_a - Tair_b)*((j - t_b)/(t_a - t_b))   #calculate Tair for particular j
      Tot_Tair <- M_Tair + Tair_int    #calculate total Tair in that month
    } 
    M_tair <- Tot_Tair/ctrx
    assign(x = paste("Mean_tair", i, sep = ""), value = M_Tair)

  }
  return(data.frame(Mean_tair1, Mean_tair2, Mean_tair3, Mean_tair4, Mean_tair5, Mean_tair6, Mean_tair7, Mean_tair8, Mean_tair9, Mean_tair10, Mean_tair11, Mean_tair12))  
 }

x <- adply(.data = A, .margins = 1, .fun = interpolate)    #adply faster, but generates entire files
x1 <- ddply(.data = A, .variables = .(lon, lat), .fun = interpolate)  # slower, but only get our results 
write.csv(x = x1, file = "Air_temp/Tair_database.csv", col.names = T)



for (i in 1:720){        #will need modification as lat, lon represents the actual degree not the row/column number
  for (j in 1:360){
    M_t[lon==i & lat==j]                   #get MODIS solar time
    t_b <- max( WFDEI_t[WFDEI_t<M_t] )  #get WFDEI solar time before MODIS solar time
    t_a <- min( WFDEI_t[WFDEI_t>=M_t])   #get WFDEI solar time after MODIS solar time
    Tair_b <- (paste("Tair_",(which(A == t_b & lon == i & lat == j, arr.ind = TRUE))[2], sep = ""))[lon == i & lat ==j]##need mod
    Tair_a <- (paste("Tair_",(which(A == t_a & lon == i & lat == j, arr.ind = TRUE))[2], sep = ""))[lon == i & lat ==j]##need mod
    Tair_M <- Tair_b + (Tair_a - Tair_b)*((M_t - t_b)/(t_a - t_b))
    A$M_Tair <- Tair_M[lon==i & lat==j]
  }
}


A<-matrix(data = rep(1:10, each = 4), nrow = 10, ncol = 4, dimnames = "X","y")

