##Basic Information####
#Trying with Terminal
#Trying with Colin
##' @title SDM analyses and projection of VeDu in future climate scenarios.
##' @date 09/06/2016
##' @author maguire.dory@gmail.com
##' @description Based on native and invasive XY data. Projecting based on climate variables alone (Bio 1, 3, 4, 5, 6, 12).
##########################################################################

# Section One Housekeeping ---------------------------------
##' Load the library
rm(list = ls())

##' Load the packages
install.packages("biomod2",dependencies=TRUE)
install.packages("dplyr",dependencies=TRUE)
install.packages("rgeos",dependencies=TRUE)
install.packages("sp",dependencies=TRUE)
library(biomod2)
library(dplyr)

#Load the data####

# CURRENT species occurrences in both invasive and native range. Invasive refers to USA data, Current referes to current native range data. Historicl data are not included in these analyses.
DataSpeciesInvasive <- read.table("Invasive.txt", header=T)
DataSpeciesCurrent <- read.table("Current2.txt", header=T)
DataSpeciesBoth=read.table("All XY Data.txt", header=T)

attach(DataSpeciesInvasive)
head(DataSpeciesInvasive)
str(DataSpeciesInvasive)

attach(DataSpeciesCurrent)
head(DataSpeciesCurrent)
str(DataSpeciesCurrent)

attach(DataSpeciesBoth)
head(DataSpeciesBoth)
str(DataSpeciesBoth)

##GET RID OF SPATIAL SAMPLING BIAS####
##' To reduce spatial sampling bias in the occurence data, do the following:
require(rgeos)
require(sp)
# The Function:
filterByProximity <- function(xy, dist, mapUnits = F) {
  #xy can be either a SpatialPoints or SPDF object, or a matrix
  #dist is in km if mapUnits=F, in mapUnits otherwise
  if (!mapUnits) {
    d <- spDists(xy,longlat=T)
  }
  if (mapUnits) {
    d <- spDists(xy,longlat=F)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    return(xy[-discard[,1],])
  }
  if (nrow(closePts) == 0) {
    return(xy)
  }
}

##' Test the function (need to make the data.frame into matrix first)
#Need to make DataSpeciesInvasive into a simple XY matrix

######FILTERED CURRENT XY Coords in BOTH NATIVE AND INVASIVE RANGE
matrixXYBoth=as.matrix(DataSpeciesBoth[,c("Y4","X4")])
FilteredBoth <- filterByProximity(matrixXYBoth, dist=5, mapUnits=F)
FilteredBoth
plot(FilteredBoth)
#Need to add a column called "Ventenata" full of 1's.
DataSpeciesBothFiltered=cbind(FilteredBoth, rep(1,663))
DataSpeciesBothFiltered
colnames(DataSpeciesBothFiltered)[3:3]=c("Ventenata4")
DataSpeciesBothFiltered
#put into dataframe
DataSpeciesBothFiltered=as.data.frame(DataSpeciesBothFiltered)
DataSpeciesBothFiltered
plot(DataSpeciesBothFiltered)

##' GET DATA READY FOR BIOMOD####
##' Put occurence data into proper format for Biomod
# the name of studied species
myRespName <- "Ventenata"

# the presence/absences data for our species
myRespBoth <- as.numeric(DataSpeciesBothFiltered[,"Ventenata4"])

# the XY coordinates of species data
myRespXYBoth <- DataSpeciesBothFiltered[,c("Y4","X4")]


# transform the presence/absence dataset into a presence only dataset to match with your own input data

pres.idBoth <- which(myRespBoth == 1)
myRespBoth <- myRespBoth[pres.idBoth]
myRespXYBoth <- myRespXYBoth[pres.idBoth,]
myRespXYBoth

################# SELECT THE ENVIRONMENTAL VARIABLES
##' Stack all the data (same climate data for both this time)
##' Data from Worldclim dataset (http://www.worldclim.org/CMIP5v1). For information on formatting these data for use in Biomod see "Data Preparation.R" Data Formatting 3.
## Select the climate variables you wish to run analyses with. Be cautious not to select too many correlated variables to avoid overfitting the model. Select variables based on biological significance, and variables that are unlikely to be heavily correlated. Data must be a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to build your models.

## CURRENT CLIMATE
# select uncorrelated climate data (SAME FOR PAST AND FUTURE!)

library(rgdal)
require(raster)
myExplCurrent=stack(c("bio_1.grd", "bio_3.grd", "bio_4.grd", "bio_5.grd", "bio_6.grd", "bio_12.grd"))

##'FUTURE CLIMATE SCENARIO (have to do this for each scenario data)
##'Here we use as an example GCM: GISS-E2-R, for RCP 60, in year 2070. 
##'This should be done for 2-3 GCMs, and 2-4 RCPs. 

myExplFuture1GS6070=stack(c("./gs60bi701.grd", "./gs60bi703.grd", "./gs60bi704.grd", "./gs60bi705.grd", "./gs60bi706.grd", "./gs60bi7012.grd"))

##' need to ensure the layer names match the current data layer names
names(myExplFuture1GS6070)=c("bio_1", "bio_3", "bio_4", "bio_5", "bio_6", "bio_12")
print(myExplFuture1GS6070)

##' Here you will add all the future climate scenarios data...

###################################################################################################################################################################
# Section Two BIOMOD Modeling ---------------------------------

#### First need to do the modeling based on the current climate data and Vedu localities. Then we will reproject the model on the future climate data after.
##'CURRENT NATIVE AND INVASIVE RANGES TOGETHER
myBiomodDataBoth <- BIOMOD_FormatingData(resp.var = myRespBoth,
                                         expl.var = myExplCurrent,
                                         resp.xy = myRespXYBoth,
                                         resp.name = "FinalBoth.Disk",
                                         PA.nb.rep = 2,
                                         PA.strategy = 'disk', ## here you can choose the type of PA sampling you want to do
                                         PA.nb.absences = 663, ## the following args will only be considered if you are using 'disk' strategy
                                         PA.dist.min = 50000,
                                         PA.dist.max = 200000)


# check that the data are read correctly
myBiomodDataBoth
plot(myBiomodDataBoth)

##'Defining Model Options

##' to see all individual model options
#Print_Default_ModelingOptions()

### to play around with parameters in each individual model.
### e.g. myBiomodOption2 <- BIOMOD_ModelingOptions((GLM = list( type = 'quadratic', interaction.level = 0)), (RF=...) , (SRE=...))

# here we will consider a quadratic GLM with first order interaction
myBiomodOption1 <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',
                                                      interaction.level = 1))



##'Modelling
##' Create the myBiomodModelOut objects for current data

# both native and invasive range data (because from prior analyses these models best explained the data).
myBiomodModelOutBoth <- BIOMOD_Modeling( myBiomodDataBoth,
                                         models = c('SRE','RF', 'GLM', 'GAM', 'MARS', 'MAXENT.Tsuruoka', 'CTA'),
                                         models.options = myBiomodOption1,
                                         NbRunEval=2,
                                         DataSplit=80,
                                         VarImport=3,
                                         models.eval.meth = c('TSS','ROC'),
                                         do.full.models=FALSE,
                                         modeling.id="test")


# print modelling summaries

myBiomodModelOutBoth

##' get model evaluations
# past
eval.dfBoth <- get_evaluations(myBiomodModelOutBoth, as.data.frame = TRUE)
# print the mean evaluation by type of model
eval.dfBoth %>% mutate(model = sub("_.*$", "", Model.name)) %>%
  group_by(model, Eval.metric) %>%
  summarise(mean.score = mean(Testing.data))



##' get variable importance for each model
# past
viBoth <- get_variables_importance(myBiomodModelOutBoth, as.data.frame = TRUE)
# print the mean vi by model
apply(viBoth, c(1,2), mean)
# print the mean vi lumping all the models together
apply(viBoth, c(1), mean)


##'Ensemble modelling
### Here we will construct the full ensemble models (based on all built models). We can keep only models having a TSS > 0.8 (could also do by weights...)

##' creates the myBiomodEM object
# past
myBiomodEMBoth <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOutBoth,
                                           chosen.models = 'all',
                                           em.by = 'all',
                                           eval.metric = c('TSS'),
                                           eval.metric.quality.threshold = c(0.75),
                                           models.eval.meth = c('TSS','ROC'),
                                           prob.mean = TRUE,
                                           prob.cv = TRUE,
                                           prob.ci = FALSE,
                                           prob.ci.alpha = 0.05,
                                           prob.median = FALSE,
                                           committee.averaging = TRUE,
                                           prob.mean.weight = TRUE,
                                           prob.mean.weight.decay = 'proportional' )


# get model evaluations
# past
eval.emBoth <- get_evaluations(myBiomodEMBoth, as.data.frame = TRUE)
eval.emBoth

#################### PROJECTIONS ##########################
##'Create Individual Model Projections
##'
##'BASED ON THE FUTURE CLIMATE SCENARIO

# future
myBiomodProjectionBoth <- BIOMOD_Projection(modeling.output = myBiomodModelOutBoth,
                                            new.env = myExplFuture1GS6070,
                                            proj.name = 'future projections 4b',
                                            selected.models = 'all',
                                            binary.meth = 'TSS',
                                            compress = FALSE,
                                            output.format = ".img",
                                            do.stack = FALSE,
                                            build.clamping.mask = FALSE)

##'Do Ensemble Model Projections
##'
# PROJECTIONS ON FUTURE CLIMATE SCENARIO
bm.efBoth <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProjectionBoth,
                                         EM.output = myBiomodEMBoth,
                                         output.format = ".img",
                                         do.stack = FALSE, binary=T)


# get the rasterstack of ensemble model projections
pred.efBoth <- get_predictions(bm.efBoth)
pred.efBoth


##########################################################################################################################################
##########################
##### FINAL PLOTS #####

library(rasterVis)
library(ggplot2)

##'Future Scenario 1a
##'
predFuture1 <- raster("./FinalBoth.Disk/proj_future projections 4b/individual_projections/FinalBoth.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
plot(predFuture1/10)
#map.r <- rasterToPoints(predFuture1)
#df <- data.frame(map.r)
#str(df)

new.raster=predFuture1
res(new.raster)=0.25 
B=resample(predFuture1, new.raster, method="bilinear")
plot(B)
map.r <- rasterToPoints(B)
df <- data.frame(map.r)
str(df)




#Make appropriate column headings
colnames(df)=c("Longitude", "Latitude", "MAP")

#Make lims
#xlimsN = c(-10, 50)
#ylimsN = c(30,60)
xlimsI=c(-130, -60)
ylimsI=c(20, 60)


#the XY coordinates of points (not necessary)
#myRespXYBoth
#colnames(myRespXYBoth)=c("x", "y")

#### Plots the whole world predicted distribution in future climate scenarios. *Change the x and ylims to display native and invasive ranges only.
#1a. Future scenario 1a
setEPS()
postscript("Future4b.eps")
ggplot(data=df, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=MAP), interpolate = TRUE) +
  borders(fill=NA, col = "grey20", size=0.05) +
  theme_bw() +
  coord_equal() +
  scale_fill_gradient("CA", low="yellow", high="red") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16, angle=90),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = ("right"),
        legend.key = element_blank()) +
  coord_fixed(xlim=xlimsI, ylim=ylimsI)
dev.off()

