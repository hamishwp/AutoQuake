library(openxlsx)
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(dplyr)
# source('../../IIDIPUS/RCode/Functions.R')
library(magrittr)
library(wbstats)
library(wid)
library(reshape2)

FilterKummu<-function(GDP,bbox,melted=F){
  
  lat<-as.numeric(colnames(GDP))
  lon<-as.numeric(rownames(GDP))
  nlat<- length(lat)
  nlon<- length(lon)
  
  imnlon<-which.min(abs(lon-bbox[1]))
  if(imnlon>1) imnlon<-imnlon-1
  imxlon<-which.min(abs(lon-bbox[3]))
  if(imxlon<nlon) imxlon<-imxlon+1
  imnlat<-which.min(abs(lat-bbox[2]))
  if(imnlat>1) imnlat<-imnlat-1
  imxlat<-which.min(abs(lat-bbox[4]))
  if(imxlat<nlat) imxlat<-imxlat+1
  
  # GDP_PPP[longitude,latitude]
  if(!melted) return(GDP[imnlon:imxlon,imnlat:imxlat])
  GDP<-GDP[imnlon:imxlon,imnlat:imxlat]
  GDP<-melt(GDP);colnames(GDP)<-c("X","Y","data")
  return(GDP)
  
}

GetKummu<-function(dir,bbox=NULL,yr=2015L){
  
  iii<-yr-1989L
  
  # file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_PPP_30arcsec_v3.nc")
  file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc")
  GDP<-brick(file,varname="GDP_per_capita_PPP")
  GDP<-GDP[[iii]]
  
  if(!is.null(bbox)) {
    e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
    crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    GDP%<>%raster::crop(e)
  }
  
  GDP%<>%as('SpatialPixelsDataFrame')
  
  return(GDP)
  
}

WBcall<-function(syear,fyear=NULL,indicator){
  if(is.null(fyear)) fyear<-syear
  value<-wb_data(indicator = indicator, start_date = as.character(syear-5L), end_date = as.character(fyear),date_as_class_date = T)
  value%>%transmute(iso3=iso3c,date=date,value=get(indicator))
}

# "EN.POP.DNST" -> Population Density
GetWBPDens<-function(syear,fyear){
  return(WBcall(syear,fyear,"EN.POP.DNST"))
}

# "SP.POP.TOTL" -> Total Population
GetWBPop<-function(syear,fyear){
  return(WBcall(syear,fyear,"SP.POP.TOTL"))
}

# "NY.GDP.MKTP.PP.KD" -> GDP-PPP (2011 US$)
GetWBGDP<-function(syear,fyear){
  return(WBcall(syear,fyear,"NY.GDP.MKTP.PP.KD"))
}

# "GCI.2NDPILLAR.XQ" -> Physical Infrastructure
GetWBInfra<-function(syear,fyear){
  return(WBcall(syear,fyear,"GCI.2NDPILLAR.XQ"))
}

# lister<-c("SI.DST.04TH.20","SI.DST.10TH.10","SI.DST.05TH.20","SI.DST.FRST.10","SI.DST.FRST.20","SI.DST.02ND.20","SI.DST.03RD.20")
# for (ind in lister){
#   ttt<-WBcall(2008,2020,ind)
#   sumz<-ttt%>%group_by(iso3)%>%summarise(nans=sum(!is.na(value)),len=length(value))
#   print(ind)
#   print(mean(sumz$nans/sumz$len))
#   print("")
#   # ggplot(sumz,aes(nans/len))+geom_histogram()
# }
# ggplot(sumz,aes(100*nans/len))+geom_histogram() + xlab("Percentage of non-empty values") +
#   ylab("Frequency") + ggtitle("Income Dist. Entries in World Bank 2008-2020") +
#   theme(plot.title = element_text(hjust = 0.5))

normaliseWB<-function(ndata,iso,dater){
  
  mindate<-min(ndata$date)
  ndata%<>%filter(iso3%in%iso&!is.na(value))%>%mutate(day=as.numeric(date-min(date)))
  # if(length(ndata$value)<=1) return(data.frame(iso3=iso,factor=1.))
  val<-data.frame()
  for (iso3c in iso){
    nd<-filter(ndata,iso3==iso3c)
    if(all(is.na(nd$value))|length(nd$value)<=1) {
      print(paste0("Not enough data found for country ",iso3c," for normalisation spline for country indicators, factor set as 1"))
      # print(nd$value)
      val%<>%rbind(data.frame(iso3=iso3c,factor=1.))
      next
    }
    func = tryCatch(splinefun(x=nd$day,y=nd$value, method="natural"),error = function(e) NULL)
    if(is.null(func)) {
      stop(paste0("No spline function over WB data possible for country ",iso3c))
    }
    val%<>%rbind(data.frame(iso3=iso3c,factor=func(as.numeric(dater-min(nd$date)))/func(as.numeric(as.Date("2015-01-01")-min(nd$date)))))
  }
  return(val)
}

# if(!is.null(date)) {
#   year<-AsYear(date)
#   nPop<-GetWBPop(directory,(year-5),year)
#   factor<-normaliseWB(nPop,iso = iso,date = date)}
# else factor<-1
# population<-population*factor

InterpWB<-function(iso3c,date,funcy){
  year<-AsYear(date)
  dataz<-do.call(funcy,list(syear=(year-5),fyear=year))
  normaliseWB(dataz,iso = iso3c,dater = date)
}
InterpPopWB<-function(iso3c,date){
  return(InterpWB(iso3c,date,GetWBPop))
}
InterpGDPWB<-function(iso3c,date){
  return(InterpWB(iso3c,date,GetWBGDP))
}

GetWID_perc<-function(perc,iso3c,year){
  
  # Note that 'j' refers to the income divided equally between spouses 
  # (only chosen because it is the dataset with largest number of entries)
  WID<-download_wid(indicators = "sptinc",years=as.character(year),perc = as.character(perc), pop = "j")
  # Filter by most popular variable (usually 'sptinc992j')
  # WID%<>%filter(variable==names(which.max(table(WID$variable))))%>%
  WID%<>%filter(variable=='sptinc992j')%>%
    mutate(iso3=convIso2Iso3(country))%>%
    dplyr::select(-c(variable,country,year)) %>%
    filter(iso3%in%iso3c)
  
  WID$value<-100*(1-WID$value)
  names(WID)[names(WID)=="percentile"]<-"variable"
  
  mins<-WID%>%group_by(iso3)%>%summarise(mins=min(value),.groups = 'drop_last')
  for(iso3c in mins$iso3) WID$value[WID$iso3==iso3c & WID$variable!="p10p100"]%<>%subtract(mins$mins[mins$iso3==iso3c])
  
  return(WID)
  
}



# 
# GetWB<-function(param,iso,Idate){
#   
#   param<-str_to_upper(gsub(" ", "", param, fixed = TRUE))
#   listy<-list("INFRASTRUCTURE"="GCI.2NDPILLAR.XQ","GDP"="NY.GDP.MKTP.PP.KD","POPULATION"="SP.POP.TOTL",
#               "POPDENSITY"="EN.POP.DNST")
#   indicator<-listy[[param]]
#   
#   if(is.null(indicator)) {
#     
#     ndata<-tryCatch(WBcall(as.character(2005),as.character(AsYear(Sys.Date())),param),error = function(e) NULL)
#     # ndata<-tryCatch(wb_data(indicator = param, start_date =as.character(2005), end_date = as.character(AsYear(Sys.Date())),date_as_class_date = T)%>%
#     # transmute(iso3=iso3c,date=date_ct,value=value),error = function(e) NULL)
#     if(is.null(ndata)) stop(paste0("ERROR: WB indicator not found '",param,"' see GetWB in GetSocioEconomic.R for examples"))
#     
#   } else {
#     
#     ndata<-tryCatch(WBcall(as.character(2005),as.character(AsYear(Sys.Date())),indicator),error = function(e) NULL)
#     # ndata<-wb_data(indicator = indicator, start_date = as.character(2005), end_date = as.character(AsYear(Sys.Date())),date_as_class_date = T)%>%
#     # transmute(iso3=iso3c,date=date_ct,value=value)
#   }
#   # ndata%<>%filter(iso3==iso&!is.na(value))%>%mutate(day=as.numeric(date-min(ndata$date)))
#   # tmp<-data.frame(iso3=iso,date=Idate)
#   if(length(iso))
#     
#     if(length(ndata$value)==0L) return(NA)
#   if(length(ndata$value)==1L) {print("WARNING: one value for WB ",indicator,", country ",iso) ;return(ndata$value)}
#   
#   func = splinefun(x=ndata$day,y=ndata$value, method="natural",  ties = mean)
#   return(func(as.numeric(Idate-min(ndata$date))))
#   
# }
# 
# GetINFORMinfo<-function(dir,year=2010){
#   
#   file<-paste0(dir,"Demography_Data/SocioEconomic/INFORM2020_TREND_2010_2020_v039_ALL.xls")
#   INFORMiso<-read.xlsx(file,sheetName = "Sheet1",header = TRUE)
#   
#   listy<-c("Population","Corruption Perception Index","Net ODA received (% of GNI)","Access to electricity",
#            "Government Effectiveness","Income Gini coefficient","Lack of Coping Capacity Index",
#            "Physical Infrastructure", "Disaster Risk Reduction", "INFORM Risk Index", "Vulnerability Index",
#            "Socio-Economic Vulnerability","Aid Dependency",
#            "FTS Current year","Estimated HDI from GDP per capita","Human Develpment Index",
#            "Physical exposure to earthquake MMI VI (relative) - raw","Physical exposure to flood (relative) - raw",
#            "Physical exposure to tropical cyclone of Saffir-Simpson category 1 (relative) - raw",
#            "Physical exposure to tsunami (relative) - raw","People affected by droughts (relative) - raw")
#   
#   INFORMiso%<>%filter(IndicatorName %in% listy & INFORMYear>year) %>% droplevels()
#   
#   drops<-c("IndicatorId","IndicatorType")
#   colnames(INFORMiso)[colnames(INFORMiso)=="Iso3"]<-"iso3"
#   INFORMiso$IndicatorName<-plyr::revalue(INFORMiso$IndicatorName,c("U5M"="Under 5 Mortality"))
#   INFORMiso<-INFORMiso[ , !(names(INFORMiso) %in% drops)]
#   
#   listy<-c("Population (total)"="POP","Corruption Perception Index"="CPI","Net ODA received (% of GNI)"="ODA",
#            "Estimated HDI from GDP per capita"="HDIGDP","Human Develpment Index"="HDI",
#            "Government Effectiveness"="GovEff","Lack of Coping Capacity Index"="CC",
#            "Physical Infrastructure"="PhysInf", "Disaster Risk Reduction"="DRR", "INFORM Risk Index"="INFORM", "Vulnerability Index"="Vuln",
#            "Socio-Economic Vulnerability"="SEVuln","Aid Dependency"="AidDep",
#            "Access to electricity"="ElecAcc","Income Gini coefficient"="GINI", "FTS Current year"="FTS",
#            "Physical exposure to earthquake MMI VI (relative) - raw"="EQexp","Physical exposure to flood (relative) - raw"="FLexp",
#            "Physical exposure to tropical cyclone of Saffir-Simpson category 1 (relative) - raw"="TCexp",
#            "Physical exposure to tsunami (relative) - raw"="TSexp","People affected by droughts (relative) - raw"="DRexp")
#   
#   INFORMiso$IndicatorName<-plyr::revalue(INFORMiso$IndicatorName, listy)
#   pc<-c("ElecAcc","CPI")
#   INFORMiso$IndicatorScore[INFORMiso$IndicatorName%in%pc]<-INFORMiso$IndicatorScore[INFORMiso$IndicatorName%in%pc]/100
#   
#   mINFORM<-data.frame()
#   for (ind in unique(INFORMiso$IndicatorName)){
#     t1<-INFORMiso%>%filter(IndicatorName==ind)
#     for (iso in unique(INFORMiso$iso3)){
#       tmp<-t1%>%filter(iso3==iso)%>%arrange(INFORMYear)
#       if(length(tmp$IndicatorScore)>0){
#         mINFORM<-rbind(mINFORM,data.frame(iso3=iso,IndicatorName=ind,dScore=max(tmp$IndicatorScore,na.rm = T)-min(tmp$IndicatorScore,na.rm = T),
#                                           wmScore=weighted.mean(x = tmp$IndicatorScore,w = 1:length(tmp$IndicatorScore))))
#       }
#     }
#     
#   }
#   for (ind in unique(mINFORM$IndicatorName)){
#     p<-ggplot(filter(mINFORM,IndicatorName==ind),aes(wmScore))+geom_density(fill="black",alpha=0.3)+
#       ylab("Density")+xlab(names(listy[listy==ind]))
#     ggsave(paste0("Density_",ind,'.png'), plot=p,path = paste0(directory,'Plots/WorldBank'),width = 7,height = 5)
#   }
#   
#   return(mINFORM)
#   
#   # corrupt<-xls %>% filter(IndicatorName=="Corruption Perception Index")
#   
# }
# 
# iso2WB<-function(iso,WBiso,ind="GINI",Score="wmScore"){
#   val<-WBiso%>%filter(iso3==iso & IndicatorName==ind)%>%pull(Score)
#   if(length(val)==0) return(NA)
#   return(val)
# }
# Viso2WB<-unname(Vectorize(iso2WB,vectorize.args = c("iso")))

# GetAllSocioEconomic<-function(dir){
#   
#   # group all data by country
#   WB<-GetWBinfo(dir)
#   WIB<-GetWIB(dir)
#   
#   filer<-paste0(dir,"Demography_Data/SocioEconomic/")
#   save(dfSE,filer)
#   
#   return(dfSE)
#   
# }

# GetWIB<-function(dir,iso3=NULL){
#   
#   file<-paste0(dir,"Demography_Data/SocioEconomic/")
#   
#   if(!is.null(iso3)) {
#     
#     return(WIB)
#   }
#   
#   # OUTPUT ALL DATA
#   
# }



# FillWBGaps<-function(listy,CM){
#   
#   nGDP<-listy[[1]]
#   nPop<-listy[[2]]
#   nPDens<-listy[[3]]
#   
#   interpy<-data.frame()
#   for (yr in unique(AsYear(CM$sdate))){
#     
#     isos<-CM%>%filter(AsYear(sdate)==yr)%>%pull(iso3)%>%unique()
#     
#     iGDP<-nGDP%>%filter(iso3%in%isos)
#     iGDP<-iGDP$iso3[is.null(iGDP[[as.character(yr)]])]
#     iPop<-nPop%>%filter(iso3%in%isos)
#     iPop<-iPop$iso3[is.null(iPop[[as.character(yr)]])]
#     iPDens<-nPDens%>%filter(iso3%in%isos)
#     iPDens<-iPDens$iso3[is.null(iPDens[[as.character(yr)]])]
#     
#     tmp<-data.frame()
#     if(length(iGDP)>0) tmp<-rbind(tmp,data.frame(iso3=iGDP,year=rep(yr,length(iGDP)),data=rep("GDP",length(iGDP))))
#     if(length(iPop)>0) tmp<-rbind(tmp,data.frame(iso3=iPop,year=rep(yr,length(iPop)),data=rep("Pop",length(iPop))))
#     if(length(iPDens)>0) tmp<-rbind(tmp,data.frame(iso3=iPDens,year=rep(yr,length(iPDens)),data=rep("PDens",length(iPDens))))
#     interpy%<>%rbind(tmp)
#     
#   }
#   rm(tmp)
#   
#   if(length(interpy)==0) return(listy)
#   
#   stop("Interpolation of WB data not yet done")
#   
# }




# GetKummu<-function(dir,bbox=NULL,yr=2015L){
#   
#   if(yr==2015L){iii<-3} else if (yr==2000L){iii<-2} else if(yr==1990L){iii<-1} else {stop("ERROR: incorrect year for GDP-PPP data")}
#   
#   # file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_PPP_30arcsec_v3.nc")
#   file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc")
#   if(!file.exists(file)) stop("no file found Kummu GDP")
#   ncin <- nc_open(file)
#   lat <- ncvar_get(ncin,"latitude")
#   nlat<- dim(lat)
#   lon <- ncvar_get(ncin,"longitude")
#   nlon<-dim(lon)
#   time <- ncvar_get(ncin,"time")
#   ntime<-dim(time)
#   
#   # GDP_PPP[longitude,latitude,time]
#   # GDP<-ncvar_get(ncin,"GDP_PPP",start = c(1,1,iii),count = c(nlon,nlat,1),collapse_degen = T)
#   GDP<-ncvar_get(ncin,"GDP_per_capita_PPP",start = c(1,1,iii),count = c(nlon,nlat,1),collapse_degen = T)
#   
#   nc_close(ncin)
#   
#   # GDP_PPP[longitude,latitude,time]
#   # GDP<-GDP[,,iii]
#   
#   colnames(GDP)<-lat
#   rownames(GDP)<-lon
#   
#   if(is.null(bbox)) {return(GDP)} else {GDP<-FilterKummu(GDP,bbox)}
#   
#   GDP%<>%convMat2SPDF(name="GDP")
#   # GDP_PPP[longitude,latitude]
#   return(GDP)
#   
# }