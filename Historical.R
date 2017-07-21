##Basic Information####

##' @title SDM analyses of VeDu in native range
##' @date 02/06/2016
##' @author maguire.dory@gmail.com
##' @description This is the analyses of the modeled distribution of Ventenata dubia in its native range of Europe. Here we compare distributions of old and new plant populations based on historical records and surveys of plant population localities, and on the environmental data available for corresponding past and present time periods. Environmental data provided by the CRU dataset. Part 1) Do the biomod ensemble models for historical and current populations, 2) Assess range shift/expansion/contraction, 3) Assess niche shift between populations using PCA. 
##########################################################################@

##Housekeeping####
##' Load the library
rm(list = ls())


##' Load the packages
library(biomod2)
library(dplyr)

#Load the data####

# species occurrences (data.frame with four columns: Subject ID, latitude, longitude, Presence(1)/Absence(0))
DataSpeciesPast <- read.table("Past2.txt", header=T)
DataSpeciesCurrent <- read.table("Current2.txt", header=T)

attach(DataSpeciesPast)
head(DataSpeciesPast)
str(DataSpeciesPast)

attach(DataSpeciesCurrent)
head(DataSpeciesCurrent)
str(DataSpeciesCurrent)

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

#PAST
matrixXYPast=as.matrix(DataSpeciesPast[,c("Y","X")])
FilteredPast <- filterByProximity(matrixXYPast,dist=5, mapUnits=F)
FilteredPast
#Need to add a column called "Ventenata" full of 1's. 
DataSpeciesPastFiltered=cbind(FilteredPast, rep(1,174))
DataSpeciesPastFiltered
colnames(DataSpeciesPastFiltered)[3:3]=c("Ventenata")
DataSpeciesPastFiltered
#put into dataframe
DataSpeciesPastFiltered=as.data.frame(DataSpeciesPastFiltered)
DataSpeciesPastFiltered

#CURRENT
matrixXYCurrent=as.matrix(DataSpeciesCurrent[,c("Y2","X2")])
FilteredCurrent <- filterByProximity(matrixXYCurrent,dist=5, mapUnits=F)
FilteredCurrent
#Need to add a column called "Ventenata" full of 1's. 
DataSpeciesCurrentFiltered=cbind(FilteredCurrent, rep(1,334))
DataSpeciesCurrentFiltered
colnames(DataSpeciesCurrentFiltered)[3:3]=c("Ventenata2")
DataSpeciesCurrentFiltered
#put into dataframe
DataSpeciesCurrentFiltered=as.data.frame(DataSpeciesCurrentFiltered)
DataSpeciesCurrentFiltered

##' GET OCCURENCE DATA READY FOR USE IN BIOMOD####
##' Put occurence data into proper format for Biomod: 2 columns matrix containing the X and Y coordinates of resp.var
# the name of studied species
myRespName <- "Ventenata"

# the presence/absences data for our species
myRespPast <- as.numeric(DataSpeciesPastFiltered[,"Ventenata"])
myRespCurrent <- as.numeric(DataSpeciesCurrentFiltered[,"Ventenata2"])

# the XY coordinates of species data
myRespXYPast <- DataSpeciesPastFiltered[,c("Y","X")]
myRespXYCurrent <- DataSpeciesCurrent[,c("Y2","X2")]

# transform the presence/absence formatted dataset (if you have absences) into a presence only dataset (which is what we have anyway)
pres.idPast <- which(myRespPast == 1)
myRespPast <- myRespPast[pres.idPast]
myRespXYPast <- myRespXYPast[pres.idPast,]
myRespXYPast

pres.idCurrent <- which(myRespCurrent == 1)
myRespCurrent <- myRespCurrent[pres.idCurrent]
myRespXYCurrent <- myRespXYCurrent[pres.idCurrent,]
myRespXYCurrent

########CLIMATE DATA #########################################################
##' ----------------------------------------------------------------------------
##' Data from CRU dataset (http://www.ipcc-data.org/observ/clim/get_30yr_means.html). For information on formatting these data for use in Biomod see "Data Preparation.R" Data Formatting 1.  
## Select the climate variables you wish to run analyses with. Be cautious not to select too many correlated variables to avoid overfitting the model. Select variables based on biological significance, and variables that are unlikely to be heavily correlated. Data must be a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to build your models.

##' Stack all the data

# past
myExplPast=stack(c("ccld0130.grd", "cdtr0130.grd", "cpre0130.grd","ctmn0130.grd","ctmp0130.grd","ctmx0130.grd","cvap0130.grd","cwet0130.grd"))

# current
myExplCurrent=stack(c("ccld0130.grd", "cdtr6190.grd", "cpre6190.grd","ctmn6190.grd","ctmp6190.grd","ctmx6190.grd","cvap6190.grd","cwet6190.grd"))

#####################################################################################
# Section Two Doing the BIOMOD Modeling ---------------------------------
##'A) HISTORICAL (PAST) DISTRIBUTION 
#Formatting Data for Biomod####

# past
myBiomodDataPast <- BIOMOD_FormatingData(resp.var = myRespPast,
                                     expl.var = myExplPast,
                                     resp.xy = myRespXYPast,
                                     resp.name = "Past.Disk",
                                     PA.nb.rep = 2,
                                     PA.strategy = 'disk', ## here you can choose the type of PA sampling you want to do
                                     PA.nb.absences = 174, ## the following args will only be considered if you are using 'disk' strategy
                                     PA.dist.min = 50000,
                                     PA.dist.max = 200000)

# current
myBiomodDataCurrent <- BIOMOD_FormatingData(resp.var = myRespCurrent,
                                         expl.var = myExplCurrent,
                                         resp.xy = myRespXYCurrent,
                                         resp.name = "Current.Disk",
                                         PA.nb.rep = 2,
                                         PA.strategy = 'disk', ## here you can choose the type of PA sampling you want to do
                                         PA.nb.absences = 334, ## the following args will only be considered if you are using 'disk' strategy
                                         PA.dist.min = 50000,
                                         PA.dist.max = 200000)


# check that the data are read correctly
myBiomodDataPast
plot(myBiomodDataPast)

myBiomodDataCurrent
plot(myBiomodDataCurrent)

#Defining Model Options####

##' to see all individual model options
Print_Default_ModelingOptions()

### to play around with parameters in each individual model.
### e.g. myBiomodOption2 <- BIOMOD_ModelingOptions((GLM = list( type = 'quadratic', interaction.level = 0)), (RF=...) , (SRE=...))

# here we will consider a quadratic GLM with first order interaction
myBiomodOption1 <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',
                                                     interaction.level = 1))



#Modelling####
##' Create the myBiomodModelOut objects for past and current data

# past
myBiomodModelOutPast <- BIOMOD_Modeling( myBiomodDataPast,
                                     models = c('SRE','RF', 'GLM', 'GAM', 'MARS', 'MAXENT.Tsuruoka', 'CTA'),
                                     models.options = myBiomodOption1,
                                     NbRunEval=2,
                                     DataSplit=80,
                                     VarImport=3,
                                     models.eval.meth = c('TSS','ROC'),
                                     do.full.models=FALSE,
                                     modeling.id="test")

# ##'B) CURRENT DISTRIBUTION 
myBiomodModelOutCurrent <- BIOMOD_Modeling( myBiomodDataCurrent,
                                         models = c('SRE','RF', 'GLM', 'GAM', 'MARS', 'MAXENT.Tsuruoka', 'CTA'),
                                         models.options = myBiomodOption1,
                                         NbRunEval=2,
                                         DataSplit=80,
                                         VarImport=3,
                                         models.eval.meth = c('TSS','ROC'),
                                         do.full.models=FALSE,
                                         modeling.id="test")

# print modelling summaries
myBiomodModelOutPast
myBiomodModelOutCurrent

##' get model evaluations
# past
eval.dfPast <- get_evaluations(myBiomodModelOutPast, as.data.frame = TRUE)
# print the mean evaluation by type of model
eval.dfPast %>% mutate(model = sub("_.*$", "", Model.name)) %>%
  group_by(model, Eval.metric) %>%
  summarise(mean.score = mean(Testing.data))

# current
eval.dfCurrent <- get_evaluations(myBiomodModelOutCurrent, as.data.frame = TRUE)
# print the mean evaluation by type of model
eval.dfCurrent %>% mutate(model = sub("_.*$", "", Model.name)) %>%
  group_by(model, Eval.metric) %>%
  summarise(mean.score = mean(Testing.data))

##' get variable importance for each model
# past
viPast <- get_variables_importance(myBiomodModelOutPast, as.data.frame = TRUE)
# print the mean vi by model
apply(viPast, c(1,2), mean)
# print the mean vi lumping all the models together
apply(viPast, c(1), mean)

# current
viCurrent <- get_variables_importance(myBiomodModelOutCurrent, as.data.frame = TRUE)
# print the mean vi by model
apply(viCurrent, c(1,2), mean)
# print the mean vi lumping all the models together
apply(viCurrent, c(1), mean)

#Ensemble modelling####
### Here we will construct the full ensemble models (based on all built models). We can keep only models having a TSS > 0.75 (could also do by weights...)

##' First, need to create the myBiomodEM object

# past
myBiomodEMPast <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOutPast,
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

# current
myBiomodEMCurrent <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOutCurrent,
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
eval.emPast <- get_evaluations(myBiomodEMPast, as.data.frame = TRUE)
eval.emPast

# current
eval.emCurrent <- get_evaluations(myBiomodEMCurrent, as.data.frame = TRUE)
eval.emCurrent

#Create Individual Model Projections####

##' do the projections for each individual model
# past
myBiomodProjectionPast <- BIOMOD_Projection(modeling.output = myBiomodModelOutPast,
                                        new.env = myExplPast,
                                        proj.name = 'past projections',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        output.format = ".img",
                                        do.stack = FALSE,
                                        build.clamping.mask = FALSE)
# current
myBiomodProjectionCurrent <- BIOMOD_Projection(modeling.output = myBiomodModelOutCurrent,
                                            new.env = myExplCurrent,
                                            proj.name = 'current projections',
                                            selected.models = 'all',
                                            binary.meth = 'TSS',
                                            compress = FALSE,
                                            output.format = ".img",
                                            do.stack = FALSE,
                                            build.clamping.mask = FALSE)


#Do Ensemble Model Projections####
# past
bm.efPast <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProjectionPast,
                                     EM.output = myBiomodEMPast,
                                     output.format = ".img",
                                     do.stack = FALSE, binary=TRUE)


# current
bm.efCurrent <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProjectionCurrent,
                                         EM.output = myBiomodEMCurrent,
                                         output.format = ".img",
                                         do.stack = FALSE, binary=TRUE)

# get the rasterstack of ensemble model projections
pred.efPast <- get_predictions(bm.efPast)
pred.efCurrent <- get_predictions(bm.efCurrent)

################################################################################################################################################################### PLOTS



###################################################################################################################################################################

##' ## ASSESSING RANGE SHIFT BETWEEN PAST AND CURRENT DATA####
##' Here we assess RANGE shift between past and current data using the BIOMOD_RangeSize function. 

#Load the binary projections first

#Raster Binary Past
pastPred=raster("./Past.Disk/proj_past projections/individual_projections/Past.Disk_EMmeanByTSS_mergedAlgo_mergedRun_mergedData_TRUEbin.img")
plot(pastPred)

#Raster Binary Current
currentPred=raster("./Current.Disk/proj_current projections/individual_projections/Current.Disk_EMmeanByTSS_mergedAlgo_mergedRun_mergedData_TRUEbin.img")
plot(currentPred)

#Call Rangesize Function
myBiomodRangeSize=BIOMOD_RangeSize(CurrentPred=pastPred, FutureProj=currentPred)
#cropped
myBiomodRangeSize=BIOMOD_RangeSize(CurrentPred=pastCropped, FutureProj=currentCropped)

#See results
myBiomodRangeSize$Compt.By.Models
setEPS()
postscript("TrialRangeSizeMap1.eps")
plot(myBiomodRangeSize$Diff.By.Pixel)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
       "lightblue")
plot(myBiomodRangeSize$Diff.By.Pixel, add=T)
dev.off()

##' We then assess the DIRECTION of the range shift between past and current data using the technique that finds the centroids of the polygon representing the past and current range of the species. 

# Need to CROP the rasters
coords=c(-20, 50, 30, 60)
pastCropped=crop(pastPred, coords)
plot(pastCropped)
currentCropped=crop(currentPred, coords)
plot(currentCropped)

# Need to convert the two EM binary prediction raster files into polygons. 
pastPoly=rasterToPolygons(pastCropped, fun=function(x){x==1}, dissolve=T)
plot(pastPoly)
currentPoly=rasterToPolygons(currentCropped, function(x){x==1}, dissolve=T)
plot(currentPoly)

# Now want to find the centroids of the polygons
library(rgeos)
pastCentroids = gCentroid(pastPoly,byid=TRUE)
plot(pastPoly)
pastCentroids
points(pastCentroids,pch=2)
#16.09718, 49.84941
currentCentroids = gCentroid(currentPoly,byid=TRUE)
plot(currentPoly)
currentCentroids
points(currentCentroids,pch=2)
#15.12831, 44.52778

# plot the two centroids
library(rworldmap)
map("world", fill=TRUE, col="white", bg="lightblue", xlim=c(-10,50), ylim=c(30, 70), mar=c(0,0,0,0))
points(pastCentroids,pch=17, col='red')
points(currentCentroids,pch=16, col='red')
###################################################################################################################################################################

##' ## ASSESSING NICHE SHIFT BETWEEN PAST AND CURRENT DATA####
##' Here we assess the NICHE shift between past and current data using PCA. We then use PERMANOVA, and between class inertia ratio to assess the significance of the shift. 

##' ## Trying out PCA stuff
##' Extract the NEW environmental variables as a matrix
# e.g. myExplCurrent is the matrix of environmental variables
# e.g. myRespXYInvasive is the data frame of XY coordinate data
# need to merge these data for past AND current data together into one matrix

myExplMatrixPast=extract(myExplPast, myRespXYPast)
myExplMatrixPast
myExplMatrixCurrent=extract(myExplCurrent, myRespXYCurrent)
myExplMatrixCurrent
Matrix=rbind(myExplMatrixPast, myExplMatrixCurrent)
Matrix
# P= past C= current
P = c(rep('P',174)) 
C = c(rep('C',334))
Groups = c(P, C)
Groups
str(Groups) 
Matrix2=cbind(Matrix, Groups)
colnames(Matrix2)[9:9]=c("Groups")
Matrix2=data.frame(Matrix2)
str(Matrix2)

################# PERFORM THE PCA

##' Does the PCA computation
# CURRENT rows 1-716, PAST rows 717-988
ClimateVariables=Matrix2[, 1:8] 
str(ClimateVariables)
PCA=prcomp(ClimateVariables, scale=T, center=T)
print(PCA)
plot(PCA)
biplot(PCA)

#PLOT the PCA
library("devtools")
install_github("kassambara/factoextra")
library("factoextra")

newplot <- fviz_pca_ind(PCA, label="none", habillage=Groups,
  addEllipses=TRUE, ellipse.level=0.95, ellipse.alpha=0)

setEPS()
postscript("PCAPastCurrent.eps")
print(newplot)
newplot + scale_color_brewer(palette="Set1") +
  theme_minimal()
dev.off()

# To get the variable biplot
fviz_pca_var(PCA)
fviz_pca_var(PCA, col.var="steelblue")+
  theme_minimal()

########### Additional stats for the PCA

## A) Performs a Monte-Carlo Test on the between-groups inertia percentage.
library(ade4)
PCA2=dudi.pca(Matrix2[, 1:8], scannf = FALSE, nf = 8)
BCA=bca(PCA2, fac=as.factor(Matrix2[, 9:9]), nf = 2)
rtest(BCA, nrepet=99)

## B) Performs a PERMANOVA test using the ADONIS function in Vegan package####
## Recall, the groups has to be a numeric value
library(vegan)
Groups=Matrix2[, 9:9]
str(Groups)
adonis=adonis(Matrix~Groups, data=Matrix2, permutations=99)
adonis
summary(adonis)
