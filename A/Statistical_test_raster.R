#extracting raster values and performing correlation rasters
#as we can extract raster values, we should be able to perform paired t-test, two sample t-test, and other similar statistical approaches on the raster
#Stats practical 4 - 7

install.packages("raster")
library(raster)

# read rasters
r1 = raster("/dir/dir/file1.tif")
r2 = raster("/dir/dir/file2.tif")
cor(getValues(r1), getValues(r2))        


# Resample r2 to r1
r2.samp = round(resample(r2, r1, "bilinear"))        #resample if the resolution is different

# plot results
# Points
plot(getValues(r2.samp) ~ getValues(r1))
# Linear regression
abline(lm(getValues(r2.samp) ~ getValues(r1)))
# (Pearson) Correlation
legend("topleft", legend=paste("Correlation =", round(cor(getValues(r1), getValues(r2.samp), use="complete.obs"), 2)))