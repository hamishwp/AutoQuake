#######################################################################
####################### ODDobj CLASS DEFINITION #######################
########### (Oxford-University Disaster Displacement Object) ##########
#######################################################################
# FIELDS:
#   - Gridded-Data - population, GDP, country ISO3C code, mean hazard 
#     intensity, stand.dev hazard intensity
#   - Hazard type
#   - Hazard mapping source
#   - Dates corresponding to hazard intensity (initialised with NA)
#   - Start date (could be evacuation initialisations)
#   - Per Iso3 extract and store indicators from World Bank, WID, ...
#   - Predicted displacement data
# METHODS:
#   - Sample hazard intensity (mean,sd) from truncated normal distribution
#   - plotGDP, plotPop, plotHaz
#   - ODDfill: AddGDP, FindISOs, ISO indicators, hazard(haztype)
#   - NOTE: hazard(haztype) does kriging for irregular grid or 
#     non-gridded data and cubic spline interpolation otherwise
#   - AffectedPop(magnitude)
#   - SampleBuildHeight(object@Population,object@buildheightpars)
#   - SampleBuildArea(object@Population,object@buildareapars)
#   - Extract hazard intensity values (I0)
#######################################################################
source('RCode/Functions.R')
source('RCode/Model.R')
source('RCode/GetPopDemo.R')
source('RCode/GetSocioEconomic.R')
source('RCode/GetINFORM.R')
library(parallel)
library(doParallel)
library(foreach)

checkODD<-function(object) {
  
  if(!is.null(object$Population) & any(object$Population<0,na.rm = T)) return(F) 
  if(!is.null(object$GDP) & any(object$GDP<0,na.rm = T)) return(F) 
  if(any(is.na(object@cIndies$value))) 
    print(paste0("WARNING: missing country indicator elements for ",object@cIndies$iso3[is.na(object@cIndies$value)]))
  
  TRUE
}

# Extract the 1D longitude & latitude elements from the grid-mesh
Genx0y0<-function(ODDobj){
  
  xo<-ODDobj@coords[1:ODDobj@grid@cells.dim[1],1]
  yo<-ODDobj@coords[1:ODDobj@grid@cells.dim[2]*ODDobj@grid@cells.dim[1]-ODDobj@grid@cells.dim[1]+1,2]
  
  return(list(xo=xo,yo=yo))
}

AddGDP<-function(ODDobj,inds=NULL,GDP=NULL){
  
  if(is.null(GDP)) GDP<-GetKummu(ODDobj@dir,c(ODDobj@bbox))
  # Minimise computation if only one GDP value is found
  if(length(unique(GDP@data$GDP))==1) {
    ODDobj$GDP<-rep(unique(GDP@data$GDP),length(ODDobj@data$Population))
    ODDobj$GDP[is.na(ODDobj$Population)]<-NA
    return(ODDobj)
  }
  ODDobj$GDP<-NA
  # interpolate data onto a regular grid
  if(!is.null(inds)) {
    ODDobj$GDP[inds]<-GDP%>%raster%>%raster::extract(ODDobj@coords[inds,])
  } else {
    ODDobj$GDP<-GDP%>%raster%>%raster::extract(ODDobj@coords)
  }
  
  # The Kummu dataset is not high resolution, 
  # Therefore, sometimes Pop data exists but GDP doesn't
  # We fix this by using the nearest-neighbour extrapolation
  GDP%<>%as.data.frame()
  for (i in which(!is.na(ODDobj$Population)&is.na(ODDobj$GDP))){
    # Find the index of the closest non-NA GDP value by longitude & latitude
    # NOTE: the end term is to ensure NA's are removed
    iminnie<-which.min((ODDobj@coords[i,1]-GDP[,2])^2*(ODDobj@coords[i,2]-GDP[,3])^2 + GDP[,1]/GDP[,1])
    ODDobj$GDP[i]<-GDP[iminnie,1]
  }
  
  return(ODDobj)
  
}

setClass("ODD", 
         slots = c(dir="character",
                   hazard="character",
                   cIndies="data.frame",
                   fIndies="list",
                   IDPs="data.frame", # includes measurement dates
                   gmax="list",
                   alerts="data.frame",
                   I0="numeric",
                   hazdates="Date",
                   eventid="numeric",
                   predictDisp="data.frame"),
         contains = "SpatialPixelsDataFrame")

# Remove all hazard intensity values that are lower than a minimum threshold
ExtractI0poly<-function(hsdf,ODD){
  
  # Extract contours
  pcontour<-adehabitatMA::getcontour(raster::subset(hsdf,mean>=ODD@I0))
  conts<-data.frame()
  id<-1
  # For each contour, extract points within only if it has a large enough area inside
  for(k in 1:length(pcontour@polygons)) {
    if(!pcontour@polygons[[k]]@Polygons[[1]]@area<1e-3) {
      conts%<>%rbind(data.frame(id=rep(id,length(pcontour@polygons[[k]]@Polygons[[1]]@coords[,1])),
                                Longitude=pcontour@polygons[[k]]@Polygons[[1]]@coords[,1],
                                Latitude=pcontour@polygons[[k]]@Polygons[[1]]@coords[,2]))
      id<-id+1
    }
    # Check for holes in the closed contours: WE DON'T LIKE DONUTS
    if(pcontour@polygons[[k]]@Polygons[[1]]@hole) 
      print(paste0("WARNING: hole in polygon of area ",
                   pcontour@polygons[[k]]@Polygons[[1]]@area,
                   " for ", ODD@hazard," event in countries: ",unique(ODD@data$ISO3C)))
  }
  
  return(conts)
}

# Add GDP data to the ODD object by interpolating onto the grid using cubic splines
setGeneric("AddHazSDF", function(ODD,lhazSDF) 
  standardGeneric("AddHazSDF") )
setMethod("AddHazSDF", "ODD", function(ODD,lhazSDF){
  
  ODD@I0<-lhazSDF$I0
  # interpolate data onto the grid
  coords<-Genx0y0(ODD)
  lenny<-length(lhazSDF) ; start<-lenny-lhazSDF$NumEvents+1
  alertscores<-alertlevels<-c() ; dates<-rep(lhazSDF$sdate,lenny-start+1)
  
  polysave<-array(F,dim=c(nrow(ODD),(lenny-start+1)))
  
  for (i in start:lenny){
    print(paste0("Hazard event number: ",i-start+1))
    hsdf<-lhazSDF[[i]]
    # Extract detail of specific hazard
    dates[i-start+1]<-hsdf@eventdate
    alertlevels%<>%c(hsdf@alertlevel)
    alertscores%<>%c(hsdf@alertscore)
    
    if(lhazSDF$hazard=="TC"){
      
      layer<-with(as.data.frame(hsdf),akima::interp(x=Longitude,y=Latitude,z=mean,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))

      layer<-c(layer$z)
      layer[layer<ODD@I0]<-NA
      
      ODD@data[paste0("hazMean",i-start+1)]<-layer
      
    } else {
      
      # extract polycontour of I<I0
      pcontour<-ExtractI0poly(hsdf=hsdf,ODD=ODD)
      # Find all ODD coordinates not inside polycontour
      insidepoly<-rep(F,nrow(ODD))
      for(p in 1:length(unique(pcontour$id))){
        tcont<-filter(pcontour,id==p)
        insidepoly<-insidepoly | sp::point.in.polygon(ODD@coords[,1],
                                                      ODD@coords[,2],
                                                      tcont$Longitude,
                                                      tcont$Latitude)>0
      }
      rm(tcont)
      hsdf%<>%as.data.frame
      # Interpolate BOTH MEAN & SD onto the ODD grid
      layer<-with(hsdf,akima::interp(x=Longitude,y=Latitude,z=mean,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))
      layer<-c(layer$z)
      layer[!insidepoly]<-NA
      if(all(is.na(layer))) next
      
      ODD@data[paste0("hazMean",i-start+1)]<-layer
      
      layer<-with(hsdf,akima::interp(x=Longitude,y=Latitude,z=sd,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))
      layer<-c(layer$z)
      layer[!insidepoly]<-NA
      ODD@data[paste0("hazSD",i-start+1)]<-layer
      
      polysave[,i-start+1]<-insidepoly
    }
  }
  
  if(lhazSDF$hazard=="TC") { 
    ind<-unname(apply(ODD@data,2, function(x) sum(!is.na(x))))
    ODD@data<-ODD@data[,ind>0]
  }else ODD$Population[rowSums(polysave)==0]<-NA
  
  # ODD@hazdates<-dates
  ODD@alerts<-data.frame(alertscores=alertscores,alertlevels=alertlevels)
  
  return(ODD)
  
})

# Initialisation of the ODD object
setMethod(f="initialize", signature="ODD",
          # definition=function(.Object,bbox,lhazSDF,dater=NULL,dir=directory,
          definition=function(.Object,lhazSDF=NULL,dir="./",Model=list(
            INFORM_vars=c("CC.INS.GOV.GE", # Government Effectiveness
                          "VU.SEV.AD", # Economic Dependency Vulnerability
                          "CC.INS.DRR", # Disaster Risk Reduction
                          "VU.SEV.PD", # Multi-dimensional Poverty
                          "CC.INF.PHY" # Physical Infrastructure
            ),
            fIndies=list(CC.INS.GOV.GE=returnX, # Government Effectiveness
                         VU.SEV.AD=returnX, # Economic Dependency Vulnerability
                         CC.INS.DRR=returnX, # Disaster Risk Reduction
                         VU.SEV.PD=returnX, # Multi-dimensional Poverty
                         CC.INF.PHY=returnX, # Physical Infrastructure
                         HA.NAT.EQ=function(x) x*x, # Exposure to the specific hazard
                         dollar=returnX, # IncomeDistribution*GDP
                         Pdens=returnX), # IncomeDistribution*GDP
            WID_perc=   c("p10p100", # top 90% share of Income Distribution
                          "p20p100", # top 80% share of Income Distribution
                          "p30p100", # top 70% share of Income Distribution
                          "p40p100", # top 60% share of Income Distribution
                          "p50p100", # top 50% share of Income Distribution
                          "p60p100", # top 40% share of Income Distribution
                          "p70p100", # top 30% share of Income Distribution
                          "p80p100", # top 20% share of Income Distribution
                          "p90p100" # top 10% share of Income Distribution
            ))) {
            if(is.null(lhazSDF)) return(.Object)
            if(!class(lhazSDF[[length(lhazSDF)]])[1]=="HAZARD") return(.Object)
            if(lhazSDF$hazard=="EQ") Model$INFORM_vars%<>%c("HA.NAT.EQ")
            else if(lhazSDF$hazard=="TC") Model$INFORM_vars%<>%c("HA.NAT.TC")
            else if(lhazSDF$hazard=="FL") Model$INFORM_vars%<>%c("HA.NAT.FL")
            else stop("Not currently prepared for hazards other than EQ, TC or FL")
            
            .Object@dir<-dir
            .Object@hazard<-lhazSDF$hazard
            
            # This bounding box is taken as the minimum region that encompasses all hazard events in HAZARD object:
            bbox<-lhazSDF$bbox
            dater<-min(lhazSDF$sdate)
            .Object@hazdates<-lhazSDF$eventdates
            
            year<-AsYear(dater)
            
            print("Fetching population data")
            obj<-GetPopulationBbox(.Object@dir,bbox=bbox)
            # I have to extract like this or else the S4 class initialisation fails!
            .Object@data <- obj@data
            .Object@coords.nrs <-obj@coords.nrs
            .Object@grid <-obj@grid
            .Object@grid.index <-obj@grid.index
            .Object@coords <-obj@coords
            .Object@bbox <-obj@bbox
            .Object@proj4string <-crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            print("Adding hazard events")
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-!is.na(.Object$Population)
            
            print("Fetching GDP-PPP data")
            .Object%<>%AddGDP(inds)
            
            print("Filter spatial data per country")
            .Object@data$ISO3C<-NA_character_
            .Object@data$ISO3C[inds]<-coords2country(.Object@coords[inds,])
            iso3c<-unique(.Object@data$ISO3C) ; iso3c<-iso3c[!is.na(iso3c)]
            
            print("Interpolate population & GDP values")
            # Note there are as many values returned as iso3c codes (returns as data.frame with columns 'iso3' and 'factor')
            Popfactors<-InterpPopWB(iso3c,dater)
            GDPfactors<-InterpGDPWB(iso3c,dater)
            for (iso in iso3c){
              indie<-.Object@data$ISO3C==iso & !is.na(.Object@data$ISO3C)
              .Object@data$Population[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
              .Object@data$GDP[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
            }
            
            # print("Extract country indicators - INFORM:")
            # INFORM (Joint Research Center - JRC) data:
            # INFORM<-InterpINFORMdata(Model$INFORM_vars,max(dater,as.Date("2014-10-22")),iso=iso3c)
            # World Income Database (WID) data:
            # if(year==AsYear(Sys.Date())) year<-AsYear(Sys.Date())-1
            # print("Extract country indicators - WID:")
            # WID<-GetWID_perc(Model$WID_perc,iso3c,year)
            # Ik<-CalcStochasticIk(WID)
            # Bind it all together!
            # .Object@cIndies<-rbind(INFORM,WID,Ik)
            # .Object@cIndies$variable%<>%as.character()
            # .Object@cIndies$iso3%<>%as.character()
            .Object@fIndies<-Model$fIndies
            
            print("Checking ODD values")
            checkODD(.Object)
            
            return(.Object)
          }
)

# setReplaceMethod( "$", signature("ODD"), function(ODD, name, value) {
#   if( name=="sides" ){
#     ODD@sides <- value
#   }
#   x
# })
# setReplaceMethod( "[", signature("ODD"), function(ODD, name, value) {
#   if( name=="sides" ){
#     ODD@data[] <- value
#   }
#   x
# })

# setGeneric("ExtractCIndy", function(ODD,iso,var)
#   standardGeneric("ExtractCIndy") )
# setMethod("ExtractCIndy", "ODD", function(ODD,iso = NULL,var=NULL){
#   cIndies<-ODD@cIndies
#   if(!is.null(iso)) cIndies%<>%filter(iso3%in%iso)
#   if(!is.null(var)) cIndies%<>%filter(variable%in%var)
#   cIndies
# })

ExtractCIndy<- function(ODD,iso = NULL,var=NULL){
  cIndies<-ODD@cIndies
  if(!is.null(iso)) cIndies%<>%filter(iso3%in%iso)
  if(!is.null(var)) cIndies%<>%filter(variable%in%var)
  cIndies
}

FormParams<-function(ODD,listy){
  listy%<>%c(list(I0=ODD@I0,fIndies=ODD@fIndies))
  Ivars<-unique(ODD@cIndies$variable)
  Params<-listy
  tParams<-list()
  for (iso3c in unique(ODD@cIndies$iso3)){
    # Extract the income distribution stochastic diffusion enhancement variable
    tParams$Ik<-ExtractCIndy(ODD,iso = iso3c,var = "Ik")$value
    # Extract all the constant country specific variables
    tParams$var<-ExtractCIndy(ODD,iso = iso3c,
                              var = Ivars[!(Ivars=="Ik" | endsWith(as.character(Ivars),"p100"))])%>%
      dplyr::select(-iso3)
    # tParams$var%<>%rbind(data.frame(value=NA,variable="dollar"))
    Params[[iso3c]]<-tParams
  }
  # names(Params)[(length(listy)+1):length(Params)]<-unique(ODD@cIndies$iso3)
  return(Params)
}

setGeneric("DispX", function(ODD,Omega,center,LL,Method)
  standardGeneric("DispX") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX", "ODD", function(ODD,Omega,center,LL=F,
                                   Method=list(Np=20,cores=8)
){
  # Extract 0D parameters & speed up loop
  Params<-FormParams(ODD,list(Np=Method$Np,center=center))
  # Income distribution percentiles & extract income percentile  
  SincN<-seq(20,90,by = 10); Sinc<-ExtractCIndy(ODD,var = paste0("p",SincN,"p100"))
  # Speed-up calculation (through accurate cpu-work distribution) to only values that are not NA
  notnans<-which(!(is.na(ODD$Population) | is.na(ODD$ISO3C) | is.na(ODD$GDP)))
  # Calculate non-local linear predictor values
  LP<-GetLP(ODD,Omega,Params,Sinc,notnans)
  # Speed things up a little
  hrange<-grep("hazMean",names(ODD),value = T) ; srange<-grep("hazSD",names(ODD),value = T)
  # Function to predict displacement per gridpoint
  CalcDisp<-function(ij){
    iso3c<-ODD@data$ISO3C[ij]
    # Calculate local linear predictor (NOTE: is a vector due to income distribution)
    locallinp<-LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]]
    # locallinp<-rep(1,10)
    # Sample population per income distribution (Assumes 9 percentiles)
    lPopS<-SplitSamplePop(Pop=ODD@data$Population[ij],Method$Np) 
    tDisp<-array(0,Method$Np)
    # for(h in hrange){
    for(ih in 1:length(hrange)){
      h<-hrange[ih]
      # for(h in c(1)){
      if(is.na(ODD@data[ij,h])) next
      # Resample population based on who is already displaced
      ind<-(colSums(lPopS)-tDisp)>0
      if(h!=hrange[1]) {
        if(sum(ind)==0) return(tDisp)
        if(length(lPopS[,ind])==0) return(tDisp)
        if(sum(ind)>1) sumz<-colSums(lPopS[,ind])
        else sumz<-sum(lPopS[,ind])
        lPopS[,!ind]<-0
        lPopS[,ind]<-SplitSamplePop(Pop=(sumz-tDisp[ind]))
      }
      # Sample hazard Intensity 
      # the uncertainty is too high... so I scale it to get some interpretable results (I know, I'm not really a statistician, I don't even have a degree, I was actually just having a look around the department when they confused me for the interviewee. I didn't have the heart to say anything. You don't hate me as much as I do)
      I_ij<-rnorm(n = Method$Np,
                  mean = ODD@data[ij,h],
                  sd = ODD@data[ij,srange[ih]]/10)
      # I_ij<-ODD@data[ij,h]
      
      # Separate into income distributions (as each have 10% of population, order doesn't matter)
      for (s in 1:length(SincN)){
        if(all(lPopS[s,]==0)) next
        # Predict damage at coordinate {i,j} (vector with MC particles)
        Damage<-tryCatch(fDamUnscaled(I_ij,Params[c("I0","Np")],Omega,Params[[iso3c]]$Ik)*locallinp[s], error=function(e) NA)
        if(any(is.na(Damage))) print(ij)
        # Scaled damage
        Dprime<-BinR(Omega$Lambda$kappa*Damage*Damage + 
                       Omega$Lambda$nu*Damage +
                       Omega$Lambda$omega,Omega$zeta)
        # Local displacement additions
        tDisp[ind]<-tDisp[ind]+Fbdisp(lPopS[s,ind],Dprime[ind])
      }  
    }
    tDisp[tDisp>ODD@data$Population[ij]]<-ODD@data$Population[ij]
    return(tDisp)
  }
  
  Disp<-array(0,c(nrow(ODD),Method$Np))
  Disp[notnans,]<-t(matrix(unlist(mclapply(X = notnans,FUN = CalcDisp,mc.cores = Method$cores)),ncol=length(notnans)))
  # Disp[notnans,]<-t(matrix(unlist(lapply(X = notnans,FUN = CalcDisp)),ncol=length(notnans)))
  
  funcy<-function(i,LLout=T) {
    tmp<-data.frame(iso3=ODD$ISO3C,IDPs=Disp[,i])%>%
      group_by(iso3)%>%summarise(predictor=floor(sum(IDPs,na.rm = T)),.groups = 'drop_last')
    tmp<-tmp[!is.na(tmp$iso3) & tmp$iso3%in%ODD@gmax$iso3,]

    if(LLout) return(sum(LL_IDP(merge(ODD@gmax,tmp,by="iso3")),na.rm = T))
    return(tmp)
  }
  
  if(LL){
    return(vapply(1:Method$Np,funcy,numeric(1)))
  }
  
  if(length(ODD@predictDisp)==0) {
    ODD@data$Disp<-rowMeans(Disp)
    tmp<-data.frame(iso3=ODD@data$ISO3C,IDPs=ODD@data$Disp)%>%
      group_by(iso3)%>%summarise(predictor=floor(sum(IDPs,na.rm = T)),.groups = 'drop_last')
    ODD@predictDisp<-tmp[!is.na(tmp$iso3)]
    return(ODD)
  }

  # Find the best fit solution
  MLE<-which.max(vapply(1:Method$Np,funcy,numeric(1)))
  # Save into ODD object
  ODD@data$Disp<-Disp[,MLE]*sum(ODD@gmax$gmax)/mean(sum(Disp[,MLE])) %>% round()
  # I know this is kind of repeating code, but I want the other function as fast as possible
  ODD@predictDisp<-merge(ODD@gmax,funcy(MLE,LLout=F),by="iso3")
  
  return(ODD)
  
})

nothingness<-function(ODDy,BDy,bbox=NULL){
  
  if(is.null(bbox)) bbox<-ODDy@bbox
  
  mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain",zoom=9)
  q<-ggmap(mad_map)
  p1<-q+ geom_raster(data=as.data.frame(ODDy),aes(Longitude,Latitude,fill=Disp),alpha=0.5,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
    scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
                         breaks=c(0,1,10,100),
                         na.value = "transparent")+#,limits=c(0,100)) +
    # labs(fill = "Displaced Pop. / Total Pop.")+xlab("Longitude") + ylab("Latitude");p
    labs(fill = "IDP Stock")+xlab("Longitude") + ylab("Latitude");p1
  
  mad_map <- get_stamenmap(bbox,source = "stamen",maptype = "terrain",zoom=9)
  q2<-ggmap(mad_map,base_layer = ggplot(as.data.frame(BDy),aes(Longitude,Latitude,group=grading)))
  p2<-q2+ geom_point(mapping = aes(colour=grading),size=0.5) + coord_cartesian() +
    # p2<-q2+ geom_jitter(mapping = aes(colour=grading),size=0.5,width = 0.008,height=0.008) + coord_cartesian() +
    xlab("Longitude") + ylab("Latitude");p2
  
  gridExtra::grid.arrange(p1,p2,nrow=1)
}

plotODDy<-function(ODDy){
  # ODDy@data$Disp<-rowMeans(predictDisp)*ODDy@gmax$gmax/mean(colSums(predictDisp))  
  
  mad_map <- get_stamenmap(ODDy@bbox,source = "stamen",maptype = "terrain",zoom=7)
  p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
  
  p+geom_contour_filled(data = as.data.frame(ODDy),
                        mapping = aes(Longitude,Latitude,z=Disp),
                        # breaks = c(0.1,1,5,10,50,300),alpha=0.5)+
                        breaks = c(0.1,1,5,10),alpha=0.5)+
    labs(fill = "Number of Displaced People")
    # geom_contour(data = as.data.frame(ODDy),
    #              mapping = aes(Longitude,Latitude,z=hazMean1,colour=..level..),alpha=1,bins = 5) +
    # scale_colour_gradient(low = "grey",high = "red",na.value = "transparent") + 
    # labs(colour = "Hazard Intensity")
  return(p)
  
}

plotODDyBG<-function(ODDy){
  # ODDy@data$Disp<-rowMeans(predictDisp)*ODDy@gmax$gmax/mean(colSums(predictDisp))  
  
  mad_map <- get_stamenmap(ODDy@bbox,source = "stamen",maptype = "terrain",zoom=7)
  p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
  
  p<-p+geom_contour_filled(data = as.data.frame(ODDy),
                           mapping = aes(Longitude,Latitude,z=Population),
                           alpha=0.5)+ 
    labs(fill = "Number of Displaced People")
  if(!is.null(ODDy@data$hazMean1)){
    p<-p+geom_contour(data = as.data.frame(ODDy),
                 mapping = aes(Longitude,Latitude,z=hazMean1,colour=..level..),alpha=0.8) +
    scale_colour_gradient(low = "transparent",high = "red",na.value = "transparent") + 
    labs(colour = "Hazard Intensity")
  }
  
  return(p)
  
}

MakeODDPlots<-function(ODDy, input){
  
    ODDy@data$Disp[ODDy@data$Disp<1]<-0
  
    mad_map <- get_stamenmap(ODDy@bbox,source = "stamen",maptype = "toner",zoom=9)
    q<-ggmap(mad_map)
    p1<-q+ geom_raster(data=as.data.frame(ODDy),aes(Longitude,Latitude,fill=Disp),
                       alpha=0.5,interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
      scale_fill_gradient2(low = "blue",mid="blue",high = "red",trans = "log",
                           # breaks=c(0,1,10,100),
                           na.value = "transparent")+#,limits=c(0,100)) +
      # labs(fill = "Displaced Pop. / Total Pop.")+xlab("Longitude") + ylab("Latitude");p
      labs(fill = "IDP Stock")+xlab("Longitude") + ylab("Latitude");p1

    p2<-q+ geom_raster(data = as.data.frame(ODDy),
                               mapping = aes(Longitude,Latitude,fill=hazMean1),
                               alpha=0.5,  inherit.aes = FALSE) + coord_cartesian() +
      scale_fill_gradient2(low = "darkorchid2",mid="olivedrab4",high = "yellow",
                           midpoint = 6,breaks=c(4.5,5,5.5,6,6.5,7,7.5,8),
                           na.value = "transparent")+#,limits=c(0,100)) +
      labs(fill = "Hazard Intensity")+xlab("Longitude") + ylab("Latitude");p2    
    # p2<-q+ geom_contour_filled(data = as.data.frame(ODDy),
    #                   mapping = aes(Longitude,Latitude,z=hazMean1,fill=..level..),
    #                   alpha=0.5, inherit.aes = FALSE,bins = 20,) + coord_cartesian() +
    #   scale_fill_discrete(low = "transparent",high = "purple",na.value = "transparent") + 
    #   labs(fill = "Hazard Intensity")+xlab("Longitude") + ylab("Latitude");p2
    
    plotty<-gridExtra::grid.arrange(p1,p2,nrow=1); plotty
        
    ggsave(paste0(input$datadir,input$plotdir,ODDy@hazard,input$sdate,input$iso3,".png"), plotty)
    
    return(plotty)
  
}

setGeneric("transIso", function(ODD) 
  standardGeneric("transIso") )
setMethod("transIso", "ODD", function(ODD)
  return(countrycode::countrycode(sourcevar = as.character(unique(ODD@cIndies$iso3)),
                                  origin = "iso3c",
                                  destination = "country.name")))
