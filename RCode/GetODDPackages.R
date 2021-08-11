library(dplyr)
library(magrittr)

GetODDPackages<-function(){

  list.of.packages <- c("ggplot2","sf","tidyverse","openxlsx","pracma","latex2exp",
                        "rJava","devtools","OpenStreetMap","sf","osmdata",
                        "tidyRSS","geojsonR", "wordcloud", "tiff", "gstat",
                        "RColorBrewer", "geosphere","GGally","FactoMineR","factoextra",
                        "gsubfn","mapsapi","leaflet", "ssh","RPostgres",
                        "rnaturalearth","rnaturalearthdata","wbstats",
                        "countrycode","rworldmap","rworldxtra","chron","ncdf4","xtable",
                        "GADMTools","akima","adehabitatMA","flexsurv")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  if(length(list.of.packages[!("ggmap" %in% installed.packages()[,"Package"])])){devtools::install_github("dkahle/ggmap")}
  if(length(list.of.packages[!("wbstats" %in% installed.packages()[,"Package"])])){devtools::install_github('nset-ornl/wbstats')}
  if(length(list.of.packages[!("wid" %in% installed.packages()[,"Package"])])){devtools::install_github("WIDworld/wid-r-tool")}
  
}

GetODDPackages_red<-function(){
  
  list.of.packages <- c("sf","tidyverse","pracma","gstat","flexsurv","mvtnorm")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
}

SourceFiles<-function(){
  
  source('RCode/Functions.R')
  source('RCode/GetInitialValues.R')
  # S4 object classes required:
  source('RCode/ODDobj.R')
  source('RCode/HAZARDobj.R')
  # Disaster related:
  source('RCode/GetUSGS.R')
  source('RCode/GetDisaster.R')
  # IDP estimate related:
  source('RCode/GetDisplacements.R')
  # Demography & population related:
  source('RCode/GetPopDemo.R')
  source('RCode/GetSocioEconomic.R')
  source('RCode/GetINFORM.R')
  # Damage estimate related:
  source('RCode/GetOSM.R')
  # Sourcing the data:
  source('RCode/GetData.R')
  # Extract model functions and priors
  source('RCode/Model.R')
  # Extract Monte Carlo algorithm
  source('RCode/Method.R')
  
}

GetODDPackages()
SourceFiles()
