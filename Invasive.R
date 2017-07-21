##Basic Information####

##' @title SDM analyses of VeDu in NATIVE vs INVASIVE range. Data from BioClim, and HFP datasets. 
##' @date 02/11/2016
##' @author maguire.dory@gmail.com 
##' @description Modeling the distribution of the plant Ventenata dubia in its native VS invasive range. 1) Creating species distribution models based on native range alone, invasive range alone, and then both together, 2) Evaluations and plots of those models, 3) Comparing niches of native populations to invasive populations using PCA analyses. 
##########################################################################
# Section One Housekeeping ---------------------------------

##' Load the library
rm(list = ls())


##' Load the packages
library(biomod2)
library(dplyr)

#Load the data####

# species occurrences (data.frame with four columns: Subject ID, latitude, longitude, Presence(1)/Absence(0))
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
##' To reduce spatial sampling bias in the occurence data, do the following to remove points that are too close to each other (e.g. within 5km radius):
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

##' To use the function above, (need to make the data.frame into matrix first)
#e.g. need to make DataSpeciesInvasive into a simple XY matrix

#INVASIVE
matrixXYInvasive=as.matrix(DataSpeciesInvasive[,c("Y3","X3")])
FilteredInvasive <- filterByProximity(matrixXYInvasive, dist=5, mapUnits=F)
str(FilteredInvasive)  
#Need to add a column called "Ventenata" full of 1's. 
DataSpeciesInvasiveFiltered=cbind(FilteredInvasive, rep(1,304))
DataSpeciesInvasiveFiltered
colnames(DataSpeciesInvasiveFiltered)[3:3]=c("Ventenata")
DataSpeciesInvasiveFiltered
#put into dataframe
DataSpeciesInvasiveFiltered=as.data.frame(DataSpeciesInvasiveFiltered)
DataSpeciesInvasiveFiltered

#NATIVE
matrixXYCurrent=as.matrix(DataSpeciesCurrent[,c("Y2","X2")])
FilteredCurrent <- filterByProximity(matrixXYCurrent,dist=5, mapUnits=F)
FilteredCurrent
#Need to add a column called "Ventenata" full of 1's. 
DataSpeciesCurrentFiltered=cbind(FilteredCurrent, rep(1,361))
DataSpeciesCurrentFiltered
colnames(DataSpeciesCurrentFiltered)[3:3]=c("Ventenata2")
DataSpeciesCurrentFiltered
#put into dataframe
DataSpeciesCurrentFiltered=as.data.frame(DataSpeciesCurrentFiltered)
DataSpeciesCurrentFiltered

######BOTH
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

##' GET OCCURENCE DATA READY FOR USE IN BIOMOD####
##' Put occurence data into proper format for Biomod: 2 columns matrix containing the X and Y coordinates of resp.var
# the name of studied species

myRespName <- "Ventenata"

# just the presence/absences column for our species (just the 1's)
myRespInvasive <- as.numeric(DataSpeciesInvasiveFiltered[,"Ventenata"])
myRespCurrent <- as.numeric(DataSpeciesCurrentFiltered[,"Ventenata2"])
myRespBoth <- as.numeric(DataSpeciesBothFiltered[,"Ventenata4"])

# just the XY coordinates of species data
myRespXYInvasive <- DataSpeciesInvasiveFiltered[,c("Y3","X3")]
myRespXYCurrent <- DataSpeciesCurrentFiltered[,c("Y2","X2")]
myRespXYBoth <- DataSpeciesBothFiltered[,c("Y4","X4")]

# transform the presence/absence formatted dataset (if you have absences) into a presence only dataset (which is what we have anyway)
pres.idInvasive <- which(myRespInvasive == 1)
myRespInvasive <- myRespInvasive[pres.idInvasive]
myRespXYInvasive <- myRespXYInvasive[pres.idInvasive,]
myRespXYInvasive

pres.idCurrent <- which(myRespCurrent == 1)
myRespCurrent <- myRespCurrent[pres.idCurrent]
myRespXYCurrent <- myRespXYCurrent[pres.idCurrent,]
myRespXYCurrent

pres.idBoth <- which(myRespBoth == 1)
myRespBoth <- myRespBoth[pres.idBoth]
myRespXYBoth <- myRespXYBoth[pres.idBoth,]
myRespXYBoth

########CLIMATE DATA ##########################################################
##' Data from Worldclim dataset (http://www.worldclim.org/current), and Global Human Footprint (Geographic) database (http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic/data-download). For information on formatting these data for use in Biomod see "Data Preparation.R" Data Formatting 2.  

## Select the climate variables you wish to run analyses with. Be cautious not to select too many correlated variables to avoid overfitting the model. Select variables based on biological significance, and variables that are unlikely to be heavily correlated. Data must be a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to build your models.

##' Stack all the data (same climate data for both this time)


myExplCurrent=stack(c("./bio_1.grd", "./bio_3.grd", "./bio_4.grd", "./bio_5.grd", "./bio_6.grd", "./bio_12.grd", "./HFP.grd"))


###################################################################################################################################################################
# Section Two Doing the BIOMOD Modeling ---------------------------------
##'A) INVASIVE RANGE USA

# Format the data for biomod 
myBiomodDataInvasive <- BIOMOD_FormatingData(resp.var = myRespInvasive,
                                             expl.var = myExplCurrent,
                                             resp.xy = myRespXYInvasive,
                                             resp.name = "FinalInvasive.Disk",
                                             PA.nb.rep = 2,
                                             PA.strategy = 'disk', ## here you can choose the type of PA sampling you want to do
                                             PA.nb.absences = 304, ## the following args will only be considered if you are using 'disk' strategy
                                             PA.dist.min = 50000,
                                             PA.dist.max = 200000)


# check that the data are read correctly
myBiomodDataInvasive
plot(myBiomodDataInvasive)

##' Defining Model Options

##' to see all individual model options
#Print_Default_ModelingOptions()

### to play around with parameters in each individual model.
### e.g. myBiomodOption2 <- BIOMOD_ModelingOptions((GLM = list( type = 'quadratic', interaction.level = 0)), (RF=...) , (SRE=...))

# here we will consider a quadratic GLM with first order interaction
myBiomodOption1 <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',
                                                      interaction.level = 1))



##'Modelling
##' Create the myBiomodModelOut objects for past and current data

# invasive
myBiomodModelOutInvasive <- BIOMOD_Modeling( myBiomodDataInvasive,
                                             models = c('SRE','RF', 'GLM', 'GAM', 'MARS', 'MAXENT.Tsuruoka', 'CTA'),
                                             models.options = myBiomodOption1,
                                             NbRunEval=2,
                                             DataSplit=80,
                                             VarImport=3,
                                             models.eval.meth = c('TSS','ROC'),
                                             do.full.models=FALSE,
                                             modeling.id="test")


# print modelling summaries
myBiomodModelOutInvasive

##' get model evaluations
# invasive
eval.dfInvasive <- get_evaluations(myBiomodModelOutInvasive, as.data.frame = TRUE)
# print the mean evaluation by type of model
eval.dfInvasive %>% mutate(model = sub("_.*$", "", Model.name)) %>%
  group_by(model, Eval.metric) %>%
  summarise(mean.score = mean(Testing.data))

##' get variable importance for each model
# invasive
viInvasive <- get_variables_importance(myBiomodModelOutInvasive, as.data.frame = TRUE)
# print the mean vi by model
apply(viInvasive, c(1,2), mean)
# print the mean vi lumping all the models together
apply(viInvasive, c(1), mean)


##'Ensemble modelling
### Here we will construct the full ensemble models (based on all built models). We can keep only models having a TSS > 0.75 (can also do by weights...).

##' First, need to create the myBiomodEM object

#invasive
myBiomodEMInvasive <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOutInvasive,
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
# invasive
eval.emInvasive <- get_evaluations(myBiomodEMInvasive, as.data.frame = TRUE)
eval.emInvasive

##'Create Individual Model Projections

##' do the projections for each individual model
myBiomodProjectionInvasive <- BIOMOD_Projection(modeling.output = myBiomodModelOutInvasive,
                                                new.env = myExplCurrent,
                                                proj.name = 'invasive projections',
                                                selected.models = 'all',
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                output.format = ".img",
                                                do.stack = FALSE,
                                                build.clamping.mask = FALSE)

##'Do Ensemble Model Projections
# invasive
bm.efInvasive <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProjectionInvasive,
                                             EM.output = myBiomodEMInvasive,
                                             output.format = ".img",
                                            do.stack = FALSE, binary=T)

 
# get the rasterstack of ensemble model projections
pred.efInvasive <- get_predictions(bm.efInvasive)
pred.efInvasive

######################################################################################################################################################
##' B) NATIVE RANGE EU

##'Formatting Data for Biomod
myBiomodDataNative <- BIOMOD_FormatingData(resp.var = myRespCurrent,
                                             expl.var = myExplCurrent,
                                             resp.xy = myRespXYCurrent,
                                             resp.name = "FinalNative.Disk",
                                             PA.nb.rep = 2,
                                             PA.strategy = 'disk', ## here you can choose the type of PA sampling you want to do
                                             PA.nb.absences = 361, ## the following args will only be considered if you are using 'disk' strategy
                                             PA.dist.min = 50000,
                                             PA.dist.max = 200000)


# check that the data are read correctly
myBiomodDataNative
plot(myBiomodDataNative)

##'Defining Model Options

##' to see all individual model options
#Print_Default_ModelingOptions()

### to play around with parameters in each individual model.
### e.g. myBiomodOption2 <- BIOMOD_ModelingOptions((GLM = list( type = 'quadratic', interaction.level = 0)), (RF=...) , (SRE=...))

# here we will consider a quadratic GLM with first order interaction
myBiomodOption1 <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',
                                                      interaction.level = 1))



##'Modelling
##' Create the myBiomodModelOut objects for native and invasive data

myBiomodModelOutNative <- BIOMOD_Modeling( myBiomodDataNative,
                                             models = c('SRE','RF', 'GLM', 'GAM', 'MARS', 'MAXENT.Tsuruoka', 'CTA'),
                                             models.options = myBiomodOption1,
                                             NbRunEval=2,
                                             DataSplit=80,
                                             VarImport=3,
                                             models.eval.meth = c('TSS','ROC'),
                                             do.full.models=FALSE,
                                             modeling.id="test")


# print modelling summaries
myBiomodModelOutNative

##' get model evaluations
# native
eval.dfNative <- get_evaluations(myBiomodModelOutNative, as.data.frame = TRUE)
# print the mean evaluation by type of model
eval.dfNative %>% mutate(model = sub("_.*$", "", Model.name)) %>%
  group_by(model, Eval.metric) %>%
  summarise(mean.score = mean(Testing.data))



##' get variable importance for each model
# native
viNative <- get_variables_importance(myBiomodModelOutNative, as.data.frame = TRUE)
# print the mean vi by model
apply(viNative, c(1,2), mean)
# print the mean vi lumping all the models together
apply(viNative, c(1), mean)


##'Ensemble modelling
### Here we will construct the full ensemble models (based on all built models). We can keep only models having a TSS > 0.8 (could also do by weights...)

##' creates the myBiomodEM object
# native
myBiomodEMNative <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOutNative,
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
# native
eval.emNative <- get_evaluations(myBiomodEMNative, as.data.frame = TRUE)
eval.emNative

##'Create Individual Model Projections

##' do the projections for each individual model
# native
myBiomodProjectionNative <- BIOMOD_Projection(modeling.output = myBiomodModelOutNative,
                                                new.env = myExplCurrent,
                                                proj.name = 'native projections',
                                                selected.models = 'all',
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                output.format = ".img",
                                                do.stack = FALSE,
                                                build.clamping.mask = FALSE)

##'Do Ensemble Model Projections
# native
bm.efNative <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProjectionNative,
                                             EM.output = myBiomodEMNative,
                                             output.format = ".img",
                                             do.stack = FALSE, binary=T)


# get the rasterstack of ensemble model projections
pred.efNative <- get_predictions(bm.efNative)
pred.efNative

######################################################################################################################################################
##'NATIVE AND INVASIVE RANGES TOGETHER
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
##' Create the myBiomodModelOut objects for past and current data

# both
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
# both
eval.dfBoth <- get_evaluations(myBiomodModelOutBoth, as.data.frame = TRUE)
# print the mean evaluation by type of model
eval.dfBoth %>% mutate(model = sub("_.*$", "", Model.name)) %>%
  group_by(model, Eval.metric) %>%
  summarise(mean.score = mean(Testing.data))



##' get variable importance for each model
# both
viBoth <- get_variables_importance(myBiomodModelOutBoth, as.data.frame = TRUE)
# print the mean vi by model
apply(viBoth, c(1,2), mean)
# print the mean vi lumping all the models together
apply(viBoth, c(1), mean)


##'Ensemble modelling
### Here we will construct the full ensemble models (based on all built models). We can keep only models having a TSS > 0.8 (could also do by weights...)

##' creates the myBiomodEM object
# both
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
# both
eval.emBoth <- get_evaluations(myBiomodEMBoth, as.data.frame = TRUE)
eval.emBoth

##'Create Individual Model Projections

##' do the projections for each individual model
# past
myBiomodProjectionBoth <- BIOMOD_Projection(modeling.output = myBiomodModelOutBoth,
                                                new.env = myExplCurrent,
                                                proj.name = 'both projections',
                                                selected.models = 'all',
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                output.format = ".img",
                                                do.stack = FALSE,
                                                build.clamping.mask = FALSE)

##'Do Ensemble Model Projections
# both
bm.efBoth <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProjectionBoth,
                                             EM.output = myBiomodEMBoth,
                                             output.format = ".img",
                                             do.stack = FALSE, binary=T)


# get the rasterstack of ensemble model projections
pred.efBoth <- get_predictions(bm.efBoth)
pred.efBoth

######################################################################################################################################################
###############
##'BOYCE INDEX####
##' To see how well our predictions match the actual XY coordinate data of our species in the invasive range, we can use the "Boyce Index". Does the evaluations of the boyce index for the invasive XY data, on the raster predictions (calibrated in invasive, native, and both ranges) provided by the ensemble model output. 
install.packages("~/Downloads/ecospat_2.1.1.tar", repos = NULL, type="binary", dependencies=TRUE)
library(ecospat)
require(ecospat)

##' 1.US XY data on Native ensemble model predictions. 

##'fit: ./Native.Disk/proj_native projections/individual_projections/Native.Disk_EMmeanByTSS_mergedAlgo_mergedRun_mergedData.img
##' obs: myRespInvasiveXY
##' 
XY=myRespXYInvasive
R=raster("./FinalNative.Disk/proj_native projections/individual_projections/FinalNative.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
R2=raster("./FinalInvasive.Disk/proj_invasive projections/individual_projections/FinalInvasive.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
R3=raster("./FinalBoth.Disk/proj_both projections/individual_projections/FinalBoth.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

#Quickly plot to check the data
plot(R)
points(XY, cex=0.3)

BoyceNative=ecospat.boyce (fit = R, obs= XY, nclass=0,
                             window.w="default", res=100, PEplot = TRUE)
summary(BoyceNative)
BoyceNative

BoyceInvasive=ecospat.boyce (fit = R2, obs= XY, nclass=0,
                           window.w="default", res=100, PEplot = TRUE)
summary(BoyceInvasive)
BoyceInvasive

BoyceBoth=ecospat.boyce (fit = R3, obs= XY, nclass=0,
                           window.w="default", res=100, PEplot = TRUE)
summary(BoyceBoth)
BoyceBoth

######################################################################################################################################################
##'PLOTS ###### 
##'Shows the raster maps of the ensemble model outputs with XY coordinate data. Plots are cropped in the native and invasive ranges. Plots are created for ensemble models calibrated in the invasive range, native range, and both. 

##' ## CALIBRATED WITH NATIVE RANGE
# To bring up the EM model CV for current based on native range
# plot a projection (just need to select a file to plot)
predNat <- raster("./FinalNative.Disk/proj_native projections/individual_projections/FinalNative.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
plot(predNat/10)


# To crop the map in EU
# Need to CROP the rasters
coords=c(-20, 50, 30, 60)
predNatCropped=crop(predNat, coords)
plot(predNatCropped)

new.raster=predNatCropped
res(new.raster)=0.25 
B=resample(predNatCropped, new.raster, method="bilinear")
plot(B)

# Make the final plot of the EU
setEPS()
postscript("NativeEM_EU.eps")
plot(B)
points(myRespXYCurrent, pch="+", cex = .5)
dev.off()

#To crop the map in the USA
coords2=c(-150, -60, 20, 60)
predInvCropped=crop(predNat, coords2)
plot(predInvCropped)
new.raster2=predInvCropped
res(new.raster2)=0.25 
C=resample(predInvCropped, new.raster2, method="bilinear")
plot(C) 

# Make the final plot of the USA
setEPS()
postscript("NativeEM_USA.eps")
plot(C)
points(myRespXYInvasive, pch="+", cex = .5)
dev.off()

##' ## CALIBRATED WITH INVASIVE RANGE 
##' 
predInv <- raster("./FinalInvasive.Disk/proj_invasive projections/individual_projections/FinalInvasive.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
plot(predInv/10)

# To crop the map in EU
# Need to CROP the rasters
coords=c(-20, 50, 30, 60)
predInvCropped=crop(predInv, coords)
plot(predInvCropped)

new.raster3=predInvCropped
res(new.raster3)=0.25 
B2=resample(predInvCropped, new.raster3, method="bilinear")
plot(B2)

# Make the final plot of the EU
setEPS()
postscript("InvasiveEM_EU.eps")
plot(B2)
points(myRespXYCurrent, pch="+", cex = .5)
dev.off()

#To crop the map in the USA
coords2=c(-150, -60, 20, 60)
predUSACropped=crop(predInv, coords2)
plot(predUSACropped)
new.raster4=predUSACropped
res(new.raster4)=0.25 
C2=resample(predUSACropped, new.raster4, method="bilinear")
plot(C2) 

# Make the final plot of the USA
setEPS()
postscript("InvasiveEM_USA.eps")
plot(C2)
points(myRespXYInvasive, pch="+", cex = .5)
dev.off()

##################################################
##'CALIBRATED IN BOTH RANGES
##' 
predBoth <- raster("./FinalBoth.Disk/proj_both projections/individual_projections/FinalBoth.Disk_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
plot(predBoth/10)

# To crop the map in EU
# Need to CROP the rasters
coords=c(-20, 50, 30, 60)
predBothCropped=crop(predBoth, coords)
plot(predBothCropped)

new.raster3=predBothCropped
res(new.raster3)=0.25 
B2=resample(predBothCropped, new.raster3, method="bilinear")
plot(B2)

# Make the final plot of the EU
setEPS()
postscript("BothEM_EU.eps")
plot(B2)
points(myRespXYCurrent, pch="+", cex = .5)
dev.off()

#To crop the map in the USA
coords2=c(-150, -60, 20, 60)
predUSACropped=crop(predBoth, coords2)
plot(predUSACropped)
new.raster4=predUSACropped
res(new.raster4)=0.25 
C2=resample(predUSACropped, new.raster4, method="bilinear")
plot(C2) 

# Make the final plot of the USA
setEPS()
postscript("BothEM_USA.eps")
plot(C2)
points(myRespXYInvasive, pch="+", cex = .5)
dev.off()


###################################################################################################################################################################

##' ## ASSESSING NICHE SHIFT BETWEEN CURRENT invasive vs native DATA####
##' @Description Here we try to assess NICHE shift between past and current data using PCA. We then use PERMANOVA, and between class inertia ratio to assess the significance of the shift. 

##' ## Trying out PCA stuff
##' Extract the NEW environmental variables as a matrix
# e.g. myExplCurrent is the matrix of environmental variables
# e.g. myRespXYInvasive is the data frame of XY coordinate data
# need to merge these data for native AND invasive ranges together into one matrix

myExplMatrixInvasive=extract(myExplCurrent, myRespXYInvasive)
str(myExplMatrixInvasive) #302
str(na.omit(myExplMatrixInvasive)) #301
myExplMatrixCurrent=extract(myExplCurrent, myRespXYCurrent)
str(myExplMatrixCurrent) #361
str(na.omit(myExplMatrixCurrent))#360
Matrix=rbind(myExplMatrixInvasive, myExplMatrixCurrent)
Matrix
str(Matrix)
head(Matrix)
M=na.omit(Matrix) 
str(M)

# need to make the groups that distinguish invasive range from native range data
I = c(rep('I',301)) 
N = c(rep('N',360))
Groups = c(I, N)
Groups
str(Groups) 
Matrix2=cbind(M, Groups)
str(Matrix2)
colnames(Matrix2)[8:8]=c("Groups")

################# PERFORM THE PCA

##' Does the PCA computation
PCA=prcomp(M, scale=T, center=T)
print(PCA)
plot(PCA)
biplot(PCA)

#Plot the PCA 
library("devtools")
install_github("kassambara/factoextra")
library("factoextra")

newplot <- fviz_pca_ind(PCA, label="none", habillage=Groups,
                        addEllipses=TRUE, ellipse.level=0.95, ellipse.alpha=0)


setEPS()
postscript("NEWPCAInvasive.eps")
print(newplot)
newplot + scale_color_brewer(palette="Set1") +
  theme_minimal()
dev.off()


# Plot the variable biplot
setEPS()
postscript("NEWVariableBiplotInvasive.eps")
fviz_pca_var(PCA)
fviz_pca_var(PCA, col.var="steelblue")+
  theme_minimal()
dev.off()

########### Additional stats for the PCA

## A) Performs a Monte-Carlo Test on the between-groups inertia percentage (in R).
library(ade4)
PCA2=dudi.pca(M, scannf = FALSE, nf = 8)
BCA=bca(PCA2, fac=as.factor(Groups), nf = 2)
rtest(BCA, nrepet=99)

## B) Performs a PERMANOVA test using the ADONIS function in Vegan package####
## Recall, the groups has to be a numeric value
library(vegan)
Data=as.data.frame(Matrix2)
G=as.numeric(Data[,8:8])
Matrix3=cbind(M, G)
Data3=as.data.frame(Matrix3)
class(Data3)
G3=Data3[,8:8]
class(G3)
M3=Data3[,1:7]
a=adonis(M3~G3, data=Data3, permutations=99)
summary(a)
a

############ Additional information for PCA plot
# To find where first introduction of Ventenata is on the PCA plot 
# First find the row of the earliest recorded observation data point in the raw data (ID=row188, )
ind <- get_pca_ind(PCA)
ind
head(ind$coord[187:188, 1:2])
head(ind$coord[27:28, 1:2])
#row188, -0.3031707  0.2880843 1960 Worley in Kootenai County
#row27, -0.008172106  1.070150, 1957 Beauty Bay/ Coeur D'Helene Spokane county



