install.packages("ncdf")
require(ncdf)


a<-list.files("C:\\Users\\Data", "*.nc", full.names = TRUE)
for(i in 1:length(a)){
  f <- open.ncdf(a[i])
  A = get.var.ncdf(nc=f,varid="So",verbose=TRUE)
  B <- get.var.ncdf(nc=f,varid="il")
  C <- get.var.ncdf(nc=f,varid="fg")
  write.table(t(rbind(A,B,C)),file="output-all.txt")}


A<-open.ncdf("Air_temp/Tair_WFDEI_201401.nc")
temp1 = get.var.ncdf( nc = A, varid = "Tair", start=c(1, 1, 1), count = c(-1, -1, 1))  #read timestep 1
temp2 = get.var.ncdf( nc = A, varid = "Tair", start=c(1, 1, 2), count = c(-1, -1, 1)) #read timestep 2
y = seq(-89.75, 89.75, 0.5)
x = seq(-179.75, 179.75, 0.5)
cellsize=0.5
out1.asc=as.asc(temp1, xll=min(x), yll=min(y), cellsize=cellsize)
out2.asc=as.asc(temp2, xll=min(x), yll=min(y), cellsize=cellsize)
#write the ascii files to the work directory
write.asc(out1.asc, 'out1.asc')
write.asc(out2.asc, 'out2.asc')
#list the ascii files
ascfiles=c('out1.asc', 'out2.asc')

#generate a dataframe from the ascii files
tdata=asc2dataframe(ascfiles)
tdata






temp1 = get.var.ncdf( nc = A, varid = "time")  #read time


create.ncdf("Ex", vars = )
n <- dim.def.ncdf("Timezone", "hours", c(1:720), unlim=FALSE, create_dimvar=TRUE )  #create new dimension for timezone but unable to input values yet
b<- var.def.ncdf("tz", "hours", n, missval = 0, longname = "hours from UTC", prec = "single")
create.ncdf("Ex", vars = b)

image(temp1)           #trying to get them into raster format
m <- t(temp1)
m <- m[nrow(m):1,]
r <- raster(m)
plot(r)
r <-raster(temp1)
r
plot(r)
