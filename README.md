# Climate

## Description

Example data and code for 1) building predictive species distribution models using the Biomod2 package in R (Thuiller et al. 2016), and 2) assessing niche and range shifts between historical and current populations (Historical.R), between invasive and native populations (Invasive.R), and between current and future species distributions (Future.R). 

## The Data 

XY Coordinate data: 
Presence only XY coordinates of historical and current populations of the grass species Ventenata dubia. This grass is native to Europe and invasive throughout the NW United States. 

Climate data: 
Historical data- CRU 30 year means (http://www.ipcc-data.org/observ/clim/get_30yr_means.html)
Current data- Worldclim v.1.4 current (http://www.worldclim.org/current), and Global Human Footprint (Geographic) database (http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic/data-download).
Future Climate Scenarios- Worldclim v.1.4 future (http://www.worldclim.org/CMIP5v1) database. 

## Setup

1. If you have never installed gdal, you must install it using brew or from source. If you use brew, open terminal and type:



    `brew install gdal`




2. Open file (e.g. Future.R) in R Studio
3. Set current working directory to repo directory
4. Run the code all at once, or section by section to see results and plots (explained throughout code in annotations). 
5. Follow the same procedure for Historical.R and Invasive.R. 
6. Apply similar analyses to other species. 
7. For further information on formatting data for use in Biomod2 see file "DataPreparation.R". 