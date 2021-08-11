
ExtractParams<-function(haz="EQ"){
  if(haz=="EQ") return(list(I0=4.5,minmag=5))
  if(haz=="TC") return(list(I0=3))
  stop("Hazard type not recognised")
}

GetEarthquake<-function(input){
  # Extract standard (or user modified) EQ parameters
  EQparams<-ExtractParams("EQ")
  
  if(!is.null(input$USGSid)) {
    out<-GetUSGS(USGSid=input$USGSid,
                 I0=EQparams$I0,minmag=EQparams$minmag)
  } else {
    # Extract bounding box of affected countries
    bbox<-countriesbbox(input$iso3)
    # Extract all earthquakes from USGS:
    out<-GetUSGS(bbox=bbox,sdate=input$sdate,fdate = input$fdate,
                 I0=EQparams$I0,minmag=EQparams$minmag)
  }
  if(is.null(out)) stop("Earthquake not found, try giving a broader/narrower range of dates")
  return(out)
}
