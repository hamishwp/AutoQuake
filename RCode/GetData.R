library(dplyr)
library(magrittr)

# Extract the initial values
GetInitVals(directory, haz=Model$haz, Model, AlgoParams)

CheckBGData<-function(datadir){
  
  if (!file.exists(paste0(datadir,
                          "Demography_Data/SocioEconomic/KUMMU/",
                          "GDP_per_capita_PPP_1990_2015_v2.nc"))){
    cat(" Please ensure that the data are located in the correct folders! 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GDP:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          Must have the Kummu file under the folder 
          'Demography_Data/SocioEconomic/KUMMU/'
          with file name 'GDP_per_capita_PPP_1990_2015_v2.nc' 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ")
    
    stop("GDP data not found, please modify the datadir input in Main_Quake.R")
  }
  
  if (!file.exists(paste0(datadir,
                          "Demography_Data/Population/gpw-v4-population-count-2015/",
                          "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_1.asc"))){
    cat(" Please ensure that the data are located in the correct folders!
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEDACS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          Must have the SEDAC files under the folder
          'Demography_Data/Population/gpw-v4-population-count-2015/'
          with file names 'gpw_v4_population_count_[...]_2015_30_sec_*.asc' 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ")
    
    stop("Population data not found, please modify the datadir input in Main_Quake.R")
  }
  
  cat("Population & GDP Data found")
  
  folders<-c("Disaster_Data/USGS/","Plots")
  for (fol in folders){
    if(!dir.exists(paste0(datadir,fol))) dir.create(paste0(datadir,fol))
  }
  
}

ExtractData<-function(haz="EQ",dir="./",extractedData=T){
  
  if(extractedData) return(paste0(dir,"IIDIPUS_Input/ODDobjects/"))

  # Get the human displacement data from IDMC Helix & GIDD databases and other sources filtered by hazard
  DispData<-GetDisplacements(haz, saved=T)
  # Extract GDACS database on the event for further validation & alertscore benchmarking
  dfGDACS<-FilterGDACS(haz=haz,syear=min(AsYear(DispData$sdate)),fyear=max(AsYear(DispData$sdate)),red=T)
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (ev in unique(DispData$eventid)){
    # Subset displacement and disaster database objects
    miniDisp<-DispData%>%filter(eventid==ev)
    # Set some confining dates for the rest of the data to be assigned to this event
    maxdate<-miniDisp$sdate-5
    if(is.na(miniDisp$fdate)) mindate<-miniDisp$sdate+3 else mindate<-miniDisp$fdate+3
    # GDACS subset
    miniDACS<-dfGDACS%>%filter(iso3%in%unique(miniDisp$iso3) & 
                                 sdate<mindate & sdate>maxdate)
    # Match displacement and hazard data and extract hazard maps
    # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
    # (list of each, bbox-cropped to remove M < minmag)
    lhazSDF<-tryCatch(GetDisaster(miniDisp,miniDACS),error=function(e) NULL)
    if(is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDisp$iso3),
                   " ",unique(miniDisp$hazard), " ", min(miniDisp$sdate) ))
      next
    }
    
    # Create the ODD object:
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF,DispData=miniDisp),error=function(e) NULL)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDisp$iso3)[1]," ", unique(miniDisp$sdate)[1])) ;next}
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDisp$iso3)[1],
                  "_",ODDy@eventid)
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDisp$iso3) & 
                               sdate<mindate & sdate>maxdate)
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDisp$iso3)[1]," ", unique(miniDisp$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
    # Path to file and ID
    path%<>%rbind(data.frame(ODDpath=ODDpath,
                             BDpath=BDpath,
                             eventid=ODDy@eventid))
    # Save some RAM
    rm(ODDy,BDy,miniDam)
  }
  
  return(path)
  
}

FormODDyOmega<-function(dir,Model,proposed,AlgoParams){
  
  output<-data.frame()
  # Load ODD files
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects_count"),pattern=Model$haz,recursive = T,ignore.case = T)
  for(i in 1:length(ufiles)){
    ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects_count/",ufiles[i]))
    ODDy@fIndies<-Model$fIndies
    if(class(ODDy@gmax)=="list") ODDy@gmax%<>%as.data.frame.list()
    # Apply DispX
    ODDy<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center,LL = F,Method = AlgoParams),
                   error=function(e) NA)
    # if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));next}
    if(is.na(ODDy)) stop(paste0("Failed to calculate Disp LL of ",ufiles[i]))
    
    saveRDS(ODDy,paste0(dir,"IIDIPUS_Results/ODDobjects_final/",ufiles[i]))
    
    print(ODDy@predictDisp)
    output%<>%rbind(ODDy@predictDisp)
  }
  return(output)
  
}

# # Total LL MAP:
# vals<-c(0.02337845,0.03644427,0.06072838,1.61965282,4.70430026,0.59366736,0.11988817,3.81639425)
# 
# # LL_Disp MAP:
# vals<-c(0.01080671,0.08229228,0.12590510,2.28739344,3.53698688,0.52015751,0.13003755,6.02350025)
# 
# Omega$Lambda[1:3]<-vals[1:3]
# Omega$zeta[1:2]<-vals[4:5]
# Omega$theta[1]<-vals[6]
# Omega$eps[1:2]<-vals[7:8]
# 
# FormODDyOmega(dir,Model,Omega,AlgoParams)
  

ODDypreds<-function(dir,haz){
  
  output<-data.frame()
  # Load ODD files
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Results/ODDobjects_final"),pattern=haz,recursive = T,ignore.case = T)
  for(i in 1:length(ufiles)){
    # if(grepl(ufiles[i],pattern = "CHN")) next
    ODDy<-tryCatch(readRDS(paste0(dir,"IIDIPUS_Results/ODDobjects_final/",ufiles[i])),error=function(e) NA)
    if(is.na(ODDy)|nrow(ODDy@predictDisp)<1) {print(paste0("no information found for ",ufiles[i])) ; next}
    tmp<-cbind(ODDy@predictDisp,data.frame(eventid=extractnumbers(ufiles[i])[2],namer=ufiles[i]))
    
    ttt<-data.frame(iso3=unique(ODDy@cIndies$iso3))
    for (var in unique(ODDy@cIndies$variable)){
      vt<-ODDy@cIndies[ODDy@cIndies$variable==var,c("iso3","value")]
      names(vt)[2]<-var
      ttt%<>%merge(vt,by="iso3")
    }
    
    output%<>%rbind(merge(tmp,ttt,by="iso3"))
    print(ODDy@predictDisp)
  }
  
  return(output)
  
}  
  
