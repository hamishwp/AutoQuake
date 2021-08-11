library(pracma)
library(ggmap)
library(OpenStreetMap)
library(osmdata)
library(ggplot2)
library(geosphere)
library(reshape2)
library(tidyverse)
library(maps)
library(magrittr)

###############################################################################################
# SortDemoData: filters out unwanted population data outside of bounding box provided by user
###############################################################################################
SortDemoData<-function(filer,bbox=NULL){
  
  info<-readLines(filer,n = 6); info<-strsplit(info, "\\s+")
  for (i in 1:6){assign(info[[i]][1],as.numeric(info[[i]][2]))} ; rm(info);
  popdemo<-read.csv(filer,header = FALSE,skip = 6,sep = " ",na.strings = NODATA_value,colClasses = "numeric")
  # dimensions of popdemo: [decreasing(latitude),longitude]
  sizer<-dim(popdemo)
  sizer2<-c(nrows,ncols+1)
  if(!all(length(sizer)==length(sizer2)) || !all(sizer==sizer2)){stop(paste0("ERROR! Incorrect dimensions: Check the population demography file for bounding box ",bbox," in file GetPopDemo.R"))}
  
  lat<-seq(from=yllcorner+90,by=-cellsize,length.out = nrows) #@@@ LATITUDE @@@#
  long<-seq(from=xllcorner,by=cellsize,length.out = ncols)    #@@@ LONGITUDE @@@#
  colnames(popdemo)<-long
  row.names(popdemo)<-lat
  
  if(is.null(bbox)) return(popdemo %>% as.matrix() %>% pracma::rot90(-1))
  
  imnlo<-which.min(abs(bbox[1]-long))
  imxlo<-which.min(abs(bbox[3]-long))
  imnla<-which.min(abs(bbox[2]-lat))
  imxla<-which.min(abs(bbox[4]-lat))
  
  popdemo<-popdemo[imxla:imnla,imnlo:imxlo] %>% as.matrix() %>% pracma::rot90(-1)
  
  return(popdemo)
}

########################################################################################
# GetSEDACfnum: Get SEDAC file number that contains the demography data in the bounding box
########################################################################################
GetSEDACfnum<-function(long,lat){
  
  longCo<-c(-90L,0L,90L,180L)
  Llo<-longCo-long
  
  if((Llo< -270L)||(Llo>360L)||(lat>90L)||(lat< -90L)){stop("Error in the longitude and latitude values for SEDACS data: see GetPopDemo.R\n")}
  
  SEDAC<-which.min(abs(Llo))
  llg<-longCo[SEDAC]
  # mLlo values have to be less than the box.
  if(Llo[SEDAC]<0L){SEDAC<-SEDAC+1L}
  if(lat<0L){SEDAC<-SEDAC+4L}
  
  return(as.integer(c(SEDAC,llg)))
}

######################################################################################
# Provided the bounding box, this function finds, filters and returns 
# the SEDAC population/demography data for that region.
######################################################################################
ExtractSEDACS<-function(strings, bbox){
  
  directory<-strings[1]
  loc<-strings[2]
  nom<-strings[3]
  
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  if(abs(LL[1]-LR[1])>1L) {
    print("Bounding box of SEDAC population/demography data is too large, using lower resolution")
    stop("Haven't modified code to run lowres yet")
  }
  
  SEDAC  <-unique(c(UL[1],UR[1],LL[1],LR[1]))
  uniquer<- length(SEDAC)
  if(uniquer==1L){
    
    filer<-paste0(directory,loc,nom,SEDAC,".asc")
    population<-SortDemoData(filer,bbox)
    
  } else if (uniquer==2L){
    
    if(((SEDAC[1]<5L)&&(SEDAC[2]<5L))|((SEDAC[1]>4L)&&(SEDAC[2]>4L))){
      # LEFT AND RIGHT SEDACS POPULATION DATA PANELS
      # bbox - [mnlo,mnla,mxlo,mxla]
      bb1<-c(bbox[1],bbox[2],LL[2],bbox[4]) # LEFT PANEL
      bb2<-c(LL[2],bbox[2],bbox[3],bbox[4]) # RIGHT PANEL
      functionner<-"rbind"
    } else if (SEDAC[1]==(SEDAC[2]+4L)|SEDAC[1]==(SEDAC[2]-4L)) {
      # UPPER AND LOWER SEDACS POPULATION DATA PANELS
      # bbox - [mnlo,mnla,mxlo,mxla]
      bb1<-c(bbox[1],0,bbox[3],bbox[4]) # UPPER PANEL
      bb2<-c(bbox[1],bbox[2],bbox[3],0) # LOWER PANEL
      functionner<-"cbind"
    } else {stop("Something went wrong.... check SEDAC file ordering, see GetPopDemo.R\n")}
    
    filer<-paste0(directory,loc,nom,SEDAC[1],".asc")
    pop1<-SortDemoData(filer,bb1)
    filer<-paste0(directory,loc,nom,SEDAC[2],".asc")
    pop2<-SortDemoData(filer,bb2)    
    
    bindy<-match.fun(functionner)
    population<-bindy(pop1,pop2)
    
    rm(pop1,pop2)
    
  } else {
    
    print("Population/Demography bounding box is in 4 separate quadrants... please be patient while this loads. Bejinhos")
    # ALL BOUNDING BOX CORNERS ARE IN DIFFERENT QUADRANTS OF THE SEDAC FILE
    # bbox - [mnlo,mnla,mxlo,mxla]
    bb1 <- c(bbox[1],  0,        LL[2],    bbox[4])   # Upper Left
    bb2 <- c(LL[2],    0,        bbox[3],  bbox[4])   # Upper Right
    bb3 <- c(bbox[1],  bbox[2],  LL[2],    0      )   # Lower Left
    bb4 <- c(LL[2],    bbox[2],  bbox[3],  0      )   # Lower Right
    
    filer<-paste0(directory,loc,nom,UL[1],".asc")
    pop1<-SortDemoData(filer,bb1)
    filer<-paste0(directory,loc,nom,UR[1],".asc")
    pop2<-SortDemoData(filer,bb2)
    filer<-paste0(directory,loc,nom,LL[1],".asc")
    pop3<-SortDemoData(filer,bb3)
    filer<-paste0(directory,loc,nom,LR[1],".asc")
    pop4<-SortDemoData(filer,bb4)
    
    # LEFT AND RIGHT USE rbind, UPPER AND LOWER USE cbind
    population<-cbind(rbind(pop1,pop2),rbind(pop3,pop4))
    
    rm(pop1,pop2,pop3,pop4)
    
  }
  
  return(population)
  
}

######################################################################################
# Extracts the SEDAC population data for the bounding box and can also plot it
######################################################################################
GetPopulationBbox<-function(directory,bbox,density=F,lowres=FALSE,yr="2015",plotty=FALSE, namer="Population",ncity=1){
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  # DATA: NASA - SEDAC
  
  if(density){ctds<-"density"} else {ctds<-"count"}
  
  if(abs(bbox[2])>90 | abs(bbox[4])>90 | abs(bbox[1])>180 | abs(bbox[3])>180) {stop("Error: non-physical bounding box values in GetPopulationBbox")}
  
  yr%<>%as.numeric()
  if(is.na(yr)) ("Please provide year as either numeric (2015) or character value ('2015')")
  
  # if(!is.null(date)) {
  #   date%<>%try(as.Date,silent=T)
  #   if(class(date)=="try-error") {
  #     print("WARNING: date badly specified in GetPopulationBbox (should be in format '2015-01-25'), no interpolation will be performed")
  #     date<-NULL
  #   }
  # }
  
  #   GetSEDACfnum(LONG,   LAT)
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  # if(abs(LL[1]-LR[1])>1L) {
  #   print("Bounding box of SEDAC population data is too large, using lower resolution")
  #   lowres<-TRUE
  # }
  
  if (lowres){
    print("WARNING: using 2020 not 2015 lowres pop data ")
    filer<-paste0(directory,"Demography_Data/Population/",
                  "gpw_v4_population_",ctds,"_adjusted_to_2015_unwpp_country_totals_rev11_",yr,"_2pt5_min.asc")
    population<-SortDemoData(filer,bbox)
    
  } else {
    
    poploc<-paste0("Demography_Data/Population/gpw-v4-population-",ctds,"-",yr,"/")
    popnom<-paste0("gpw_v4_population_",ctds,"_adjusted_to_2015_unwpp_country_totals_rev11_",yr,"_30_sec_")
    strings<-c(directory,poploc,popnom)
    
    population<-ExtractSEDACS(strings, bbox) 
    
  }
  
  if(plotty){
    
    longData<-melt(population)
    longData<-longData[longData$value!=0,]
    
    cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
    if(ncity>1){wordcloud::wordcloud(words=cities$name,freq = cities$pop,max.words = 30,scale = c(2.5,0.2))}
    cities<-slice(cities,1:ncity)
    
    p<-ggplot(longData, aes(x = Var1, y = Var2)) + 
      geom_raster(aes(fill=log10(value))) + 
      scale_fill_gradient(low="yellow", high="red",name = "log_{10}(Population Count)") +
      labs(x="Longitude", y="Latitude" )+#, title="East Gippsland, Australia") +
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11)) +
      geom_label(data = cities, aes(long, lat, label = name), size = 4, fontface = "bold", nudge_x = 0.25,nudge_y = -0.25)
    print(p)
    ggsave(paste0(namer,"Population_",ctds,".eps"), plot=p,path = paste0(directory,'/'),width = 9,height = 7.)
    
    # mad_map <- get_map(bbox,maptype = "toner-background")
    # print(ggmap(mad_map))
    
  } 
  
  population%<>%convMat2SPDF(name="Population")
  
  return(population)
  
}

######################################################################################
# Extracts the SEDAC aging data for the bounding box and can also plot it
######################################################################################
GetAgingBbox<-function(directory,bbox,lowres=FALSE,sex="b",plotty=FALSE, namer="Aging",ncity=1){
  # bbox is bounding box in the form 'min lon, max lat, max lon, min lat'
  # DATA: NASA - SEDAC
  
  # PERFORM BASIC CHECKS
  sex<-substr(sex,1,1) ;  sex<-tolower(sex)
  if(!(sex %in% c("b","m","f"))){stop("Error: choice of aging population demography 'sex' is either 'b'=both, 'f'=female or m'=male")}
  if(abs(bbox[2])>90 | abs(bbox[4])>90 | abs(bbox[1])>180 | abs(bbox[3])>180) {stop("Error: non-physical bounding box values in GetAgingBbox")}
  
  subfiler<-paste0("gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_2pt5_min.asc")
  
  #   GetSEDACfnum(LONG,   LAT)
  LL<-GetSEDACfnum(bbox[1],bbox[2]) # Lower Left
  LR<-GetSEDACfnum(bbox[3],bbox[2]) # Lower right
  UL<-GetSEDACfnum(bbox[1],bbox[4]) # Upper Left
  UR<-GetSEDACfnum(bbox[3],bbox[4]) # Upper right
  
  if(abs(LL[1]-LR[1])>1L) {
    print("Bounding box of SEDAC aging population data is too large, using lower resolution")
    lowres<-TRUE
  }
  
  if (lowres){
    
    filer<-paste0(directory,"Demography_Data/Age/",
                  "gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_2pt5_min.asc")
    aging<-SortDemoData(filer,bbox)
    
  } else {
    
    poploc<-paste0("Demography_Data/Age/gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_30_sec")
    popnom<-paste0("gpw_v4_basic_demographic_characteristics_rev11_a065plus",sex,"t_2010_dens_30_sec_")
    strings<-c(directory,poploc,popnom)
    
    aging<-ExtractSEDACS(strings, bbox) 
    
  }
  
  if(plotty){
    
    longData<-melt(aging)
    longData<-longData[longData$value!=0,]
    
    cities<-maps::world.cities%>%filter(lat>bbox[2]&lat<bbox[4]&long>bbox[1]&long<bbox[3])%>%arrange(desc(pop))
    if(ncity>1){wordcloud::wordcloud(words=cities$name,freq = cities$pop,max.words = 30,scale = c(2.5,0.2))}
    cities<-slice(cities,1:ncity)
    
    p<-ggplot(longData, aes(x = Var1, y = Var2)) + 
      geom_raster(aes(fill=value)) + 
      scale_fill_gradient(low="grey90", high="red") +
      labs(x="Longitude", y="Latitude", title="Aging Demographic 65+") +
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11)) +
      geom_label(data = cities, aes(long, lat, label = name), size = 4, fontface = "bold", nudge_x = 1.5)#+
    print(p)
    ggsave(paste0(namer,"Aging60plusDemo.eps"), plot=p,path = paste0(directory,'/'),width = 9,height = 7.)
    
    # mad_map <- get_map(bbox,maptype = "toner-background")
    # print(ggmap(mad_map))
    
  }
  
  return(aging)
  
}

#############################################################################
# Check that place and country names can be used and don't generate rubbish
#############################################################################
PlaceCheck<-function(country,place,viewbox){
  # Place is village/town/city/region,etc but not country!
  if(strcmpi(place,country)){
    writeLines(paste0("Setting place=NULL for '",place,"' and using country '",country,"' in PlaceCheck function"))
    place=NULL
  }
  if(is.null(place)){bb <- getbb (country, viewbox = viewbox, silent = FALSE, format_out = 'polygon')
  }  else {bb <- getbb (paste0(place,", ",country), viewbox = viewbox, silent = FALSE, format_out = 'polygon')}
  if(length(bb)>1){
    for (i in 2:length(bb)){
      disty<-distm(c(mean(bb[[1]][,1]), mean(bb[[1]][,1])), c(mean(bb[[i]][,2]), mean(bb[[i]][,2])), fun = distHaversine) 
      if (disty>5000){    
        writeLines(paste0("
                       ######################################################################################### \n
                       Check the uniqueness of the place specified: '",place,", ",country,"' \n 
                       It appears that there are multiple places with that name, consider using 'viewbox' \n
                       Error location: function 'GetPopulationPlace' in file 'GetPopDemo.R' \n
                       ######################################################################################### \n"))
        stop()
      }
    }
  }
  
  return(bb[[1]])
  
}

######################################################################################
# Extracts the SEDAC population data for a given place and can also plots it
######################################################################################
GetPopulationPlace<-function(directory, country, place=NULL, lowres=FALSE,sex="b", viewbox=NULL,plotty=FALSE, namer=NULL,ncity=1){
  if(is.null(namer)){namer<-paste0(place,"_",country,"_population")}
  
  bb<-PlaceCheck(country,place,viewbox)
  
  # bbox - [mnlo,mnla,mxlo,mxla]
  bbox[1]<-min(bb[[1]][,1],na.rm = T)
  bbox[2]<-min(bb[[1]][,2],na.rm = T)
  bbox[3]<-max(bb[[1]][,1],na.rm = T)
  bbox[4]<-max(bb[[1]][,2],na.rm = T)
  # DATA: NASA - SEDAC
  population<-GetPopulationBbox(directory,bbox,lowres=FALSE,plotty=plotty,namer = namer,ncity=ncity)
  
  if(plotty){
    df<-data.frame(xx=bb[,1],yy=bb[,2])
    mad_map <- get_map(getbb(country),source="stamen",maptype = "toner-background")
    ggmap(mad_map) + coord_fixed() +
      geom_polygon(aes(x = xx, y = yy),data = df, colour = NA, fill = "red", alpha = .2)
  }
  
  stop("Error: code cannot yet integrate population in polygon")
  
  totalpop<-IntPoly(population,bb)
  
  return(totalpop)
  
}

######################################################################################
# Extracts the SEDAC population data for a given place and can also plots it
######################################################################################
GetAgingPlace<-function(directory, country, place=NULL, viewbox=NULL, plotty=FALSE, namer=NULL,ncity=1){
  if(is.null(namer)){namer<-paste0(place,"_",country,"_aging60plus")}
  
  bb<-PlaceCheck(country,place,viewbox)
  
  # bbox - [mnlo,mnla,mxlo,mxla]
  bbox[1]<-min(bb[[1]][,1],na.rm = T)
  bbox[2]<-min(bb[[1]][,2],na.rm = T)
  bbox[3]<-max(bb[[1]][,1],na.rm = T)
  bbox[4]<-max(bb[[1]][,2],na.rm = T)
  # DATA: NASA - SEDAC
  aging<-GetAgingBbox(directory,bbox,lowres=FALSE,plotty=plotty,namer = namer,ncity = ncity)
  
  if(plotty){
    df<-data.frame(xx=bb[[1]][,1],yy=bb[[1]][,2])
    mad_map <- get_map(getbb(country),source="stamen",maptype = "toner-background")
    ggmap(mad_map) + coord_fixed() +
      geom_polygon(aes(x = xx, y = yy),data = df, colour = NA, fill = "red", alpha = .2)
  }
  
  stop("Error: code cannot yet integrate aging population in polygon")
  
  avage<-Mode(aging)
  
  return(avage)
  
}
