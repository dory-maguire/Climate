##Basic Information####

##' @title FORMATTING THE CLIMATE ENVIRONMENTAL DATA FOR USE IN BIOMOD
##' @date 18/10/2016
##' @author maguire.dory@gmail.com
##' @description formatting all of the explainatory variables so that they work in Biomod2.
##########################################################################@

##Housekeeping####
##' Load the library
rm(list = ls())
setwd("/Users/dorothy_maguire/Documents/The R Folder/Preliminary Analyses Vedu")

##' Load the packages
library(biomod2)
library(dplyr)
require(raster)
require(rgdal)
library(raster)
############################################################################################################################################
##' Data Formatting 1. "Historical.R" 
##' This environmental data corresponds to the historical analyses of Ventenata dubia in its native range. We used historical and recent 30 year mean climate data from the CRU database (http://www.ipcc-data.org/observ/clim/get_30yr_means.html). Data had to be changed from .dat files to .grd data files.

# the path to the .dat file (need to do this for each file)
f <- '/Users/dorothy_maguire/Documents/The R Folder/Data/Climate data/cwet0130.dat'

# read in the file
d <- readLines(f)

# get the raster attributes
d.head.colnames <- as.character(strsplit(d[1], "[[:space:]]{1,10}")[[1]])
d.head.values <- as.numeric(strsplit(d[2], "[[:space:]]{1,10}")[[1]])
names(d.head.values) <- d.head.colnames
d.head.values

# extract the data you need
d <- d[-c(1:2)]
d <- lapply(d, function(r) lapply(0:(d.head.values['n_cols']-1), function(i) substr(r, i*5+1,(i+1)*5)))
d <- as.numeric(unlist(d))
d[d==-9999] <- NA
d <- matrix(d, ncol=d.head.values['n_months'])

# put the values into a raster
b <- brick(xmn=d.head.values['xmin'] - d.head.values['grd_sz'] / 2, xmx=d.head.values['xmax'] + d.head.values['grd_sz'] / 2,
           ymn=d.head.values['ymin'] - d.head.values['grd_sz'] / 2, ymx=d.head.values['ymax'] + d.head.values['grd_sz'] / 2,
           nl=d.head.values['n_months'])
res(b) <- d.head.values['grd_sz']
values(b) <- d
b <- rotate(b)

# check if the data are correctly read
plot(b)

##' Combine the monthly values into annual mean
b2=stackApply(b, indices=c(1,1,1,1,1,1,1,1,1,1,1,1), fun=mean, filename='b2', overwrite=T)
plot(b2mean)

##' Save the raster stack as a .grd file
writeRaster(b2, filename = sub(".dat$", ".grd", f), overwrite = TRUE)

############################################################################################################################################
##' Data Formatting 2: "Invasive.R"
##' This environmental data corresponds to the analyses of of Ventenata dubia in the invasive vs native range. These data are from the National Geographic: Human Footprint database (http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic/data-download). This is how to ready the files for use in biomod. 

##' Converting the "HFP" or Human Footprint data (.adf files originally) to use in Biomod.
library(rgdal)
library(RColorBrewer)
dpath<-"//Users/dorothy_maguire/Documents/The R Folder/Data/Climate data/HFP/hfp_global_geo_grid/hf_v2geo/hdr.adf"
x <- new("GDALReadOnlyDataset", dpath)
getDriver(x)
getDriverLongName(getDriver(x))
xx<-asSGDF_GROD(x)
r <- raster(xx)
r
plot(r)

##' Reprojects the HFP data to be the same extent and resolution as the CRU data (MUST do this to use in Biomod), and saves the data as .grd files. 

# "w" is the CRU data.
raster_data <- r
w <- raster("//Users/dorothy_maguire/Documents/The R Folder/Data/Climate data/ccld0130.grd", crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")

rp <- projectRaster(from = r, to = w,
            method = "bilinear",
            format = "raster",
            overwrite = TRUE)
rp

#then write "rp" to new raster:
writeRaster(rp, filename="HFP.grd")

############################################################################################################################################
##################################################################################' Data Formatting 3: "Future.R"
##'This environmental data corresponds to the analyses of of Ventenata dubia in given future climate scenarios. The following code can be used to make sure the future climate data from the Worldclim database (http://www.worldclim.org/CMIP5v1) are in the proper format and projected in the same format as .grd files.

##' Reprojects the future .tif bioclim data to same format, and saves to a .grd file.
r=raster("/Users/dorothy_maguire/Documents/The R Folder/Preliminary Analyses VeDu/gs60bi7012.tif")
w <- raster("//Users/dorothy_maguire/Documents/The R Folder/Preliminary Analyses VeDu/bio_1.grd")

rp <- projectRaster(from = r, to = w,
                    method = "bilinear",
                    format = "raster",
                    overwrite = TRUE)
writeRaster(rp, format="raster", filename="gs60bi7012.grd")




