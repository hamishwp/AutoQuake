##########################################################################
##########################################################################
############################ Model formulation ###########################
### There are three sections to this file:                             ###
### 1) Define model variables, parameterisation and link functions     ###
###    Choose your important variables here, e.g. unemployment rate    ###
### 2) Linear predictor: damage*exp(linearpredictor)                   ###
###    This acts to modify the damage based on the country/area values ###
###    For example, local GDP or pop density decrease expected damage  ###
### 3) Log-likelihood, prior and posterior distribution definitions    ###
###    When optimising the parameters of the model based on historical ###
###    data, this section defines what function to minimise (cost fn)  ###
##########################################################################
##########################################################################
##########################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Model variables and format, based on the specific hazard
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

Model<-list()

haz<-"EQ"
# haz<-"TC"
Model$haz<-haz

# Get linear predictor variables (currently the same for both building damage and human displacement)
if(haz%in%c("EQ","TC","FL")) {
  # Form a list of the different components used in the linear predictor
  INFORM_vars<-c("CC.INS.GOV.GE", # Government Effectiveness
                 "VU.SEV.AD", # Economic Dependency Vulnerability
                 "CC.INS.DRR", # Disaster Risk Reduction
                 "VU.SEV.PD", # Multi-dimensional Poverty
                 "CC.INF.PHY" # Physical Infrastructure
  )
  INFORM_vars%<>%c(grep(haz,c("HA.NAT.EQ","HA.NAT.TC","HA.NAT.FL"),value = T))
  # MUST BE SAME LENGTH AND CORRESPOND EXACTLY TO INFORM_vars + HAZ.NAT.EQ/TC/FL + Sinc/dollar,
  fIndies<-list(CC.INS.GOV.GE=returnX, # Government Effectiveness
                VU.SEV.AD=returnX, # Economic Dependency Vulnerability
                CC.INS.DRR=returnX, # Disaster Risk Reduction
                VU.SEV.PD=returnX, # Multi-dimensional Poverty
                CC.INF.PHY=returnX, # Physical Infrastructure
                dollar=returnX, # IncomeDistribution*GDP
                Pdens=returnX) # Population Density
  if(haz=="EQ") {fIndies%<>%c(HA.NAT.EQ=function(x) x*x) # Hazard Exposure)
  } else if (haz=="TC") {fIndies%<>%c(HA.NAT.TC=function(x) x*x) # Hazard Exposure)
  } else if (haz=="FL") {fIndies%<>%c(HA.NAT.EQ=function(x) x*x) }# Hazard Exposure)
  WID_perc<-   c("p10p100", # top 90% share of Income Distribution
                 "p20p100", # top 80% share of Income Distribution
                 "p30p100", # top 70% share of Income Distribution
                 "p40p100", # top 60% share of Income Distribution
                 "p50p100", # top 50% share of Income Distribution
                 "p60p100", # top 40% share of Income Distribution
                 "p70p100", # top 30% share of Income Distribution
                 "p80p100", # top 20% share of Income Distribution
                 "p90p100" # top 10% share of Income Distribution
  )
  
  Model%<>%c(list(INFORM_vars=INFORM_vars,WID_perc=WID_perc,fIndies=fIndies))
  
}

# Link functions (MUST BE SAME LENGTH AS OMEGA)
Model$links<-list(
  Lambda=list(kappa='exp',nu='exp',omega='exp'),
  zeta=list(k='exp',lambda='exp'), # zeta=list(k=2.5,lambda=1.6),
  beta=list(xxx='exp',CC.INS.GOV.GE='exp',VU.SEV.AD='exp',CC.INS.DRR='exp',VU.SEV.PD='exp',CC.INF.PHY='exp',dollar='negexp',Pdens='returnX'),
  theta=list(e='exp'), #list(e=0.25),
  # rho=list(A='exp',H='exp'),
  eps=list(eps='exp',xi='exp')
  # mu=list(muplus='exp',muminus='exp',sigplus='exp',sigminus='exp')
)
names(Model$links$beta)[1]<-paste0("HA.NAT.",haz)

# These should be the inverse of the link functions! (You should create a checker function)
InverseLinks<-function(Omega){
  Omega%<>%unlist()
  Omega[c(1:11,14:16)]%<>%log()
  Omega[12]<-log(-Omega[12])
  Omega%>%relist(skeleton = Model$skeleton)
}

# Using the parameterisation skeleton, apply link functions to the proposed parameters
Proposed2Physical<-function(proposed,Model,index=NULL){
  
  Model$links%<>%unlist()
  
  if(is.null(index)) index<-1:length(Model$links)
  # Link functions to convert values into useable/physical values
  for (i in index)  {
    proposed[i] <- match.fun(Model$links[[names(proposed)[i]]])(proposed[i])
  }
  # Reshape into desired structure
  proposed%>%relist(skeleton=Model$skeleton)
  
}

# Skeleton
Model$skeleton <- list(
  Lambda=list(kappa=NA,nu=NA,omega=NA),
  zeta=list(k=NA,lambda=NA), # zeta=list(k=2.5,lambda=1.6),
  beta=list(xxx=NA,CC.INS.GOV.GE=NA,VU.SEV.AD=NA,CC.INS.DRR=NA,VU.SEV.PD=NA,CC.INF.PHY=NA,dollar=NA,Pdens=NA),
  theta=list(e=NA), #list(e=0.25),
  # rho=list(A=NA,H=NA),
  eps=list(eps=NA,xi=NA)
  # mu=list(muplus=NA,muminus=NA,sigplus=NA,sigminus=NA)
)
names(Model$skeleton$beta)[1]<-paste0("HA.NAT.",haz)

# Get the binary regression function
Model$BinR<-"weibull" # "gompertz"

# Implement higher order Bayesian priors?
Model$higherpriors<-TRUE

# Get the building damage beta distribution parameters
Model$BD_params<-list(functions=list(
    destroyed="dbeta",
    severe="dbeta",
    moderate="dbeta",
    possible="dbeta",
    notaffected="dbeta"
    # damaged="alldamaged" # damaged is a combination of all categories except 'notaffected'
  ),
  Params=list(
    destroyed=list(shape1=40,shape2=2),
    severe=list(shape1=40,shape2=25),
    moderate=list(shape1=15,shape2=30),
    possible=list(shape1=3,shape2=100),
    notaffected=list(shape1=0.05,shape2=100)
    # damaged=list() # damaged is a combination of all categories except 'notaffected'
  )
)

Model$center<-ExtractCentering(dir,Model,haz,T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Linear predictor calculations (act to modify the expected damage values)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Linear predictor function - ugly, but very fast!
llinpred<-function(Params,beta,center,func) { # Note all of these are vectors, not single values
  
  linp<-list()
  for(j in 1:length(Params)) {
    vars<-Params[[j]]$var
    # This part uses the name of the variable to search for the following:
    # 1) link function - func, 2) centering value - center, 3) linear predictor parameterisation - beta
    linp[[names(Params[j])]]<-prod(exp(vapply(1:nrow(vars),
                                function(i) beta[[vars$variable[i]]]*(do.call(func[[vars$variable[i]]],
                                                                              list(vars$value[i]-center[[vars$variable[i]]] ))),
                                FUN.VALUE = numeric(1))))
  }
  return(linp)
}

# Quicker version of llinpred for when only one vars$variable exists
dlinpred<-function(vars,beta,center,func) { # Note all of these are vectors, not single values
  return(exp(beta[[vars$variable[1]]]*(do.call(func[[vars$variable[1]]],list(vars$value-center[[vars$variable[1]]] )))))
}

GDPlinp<-function(ODD,Sinc,beta,center,fIndies,notnans){
  
  iGDP<-as.numeric(factor(ODD$GDP,levels=unique(ODD$GDP)))
  dGDP<-data.frame(ind=iGDP[notnans],GDP=ODD$GDP[notnans],iso=ODD$ISO3C[notnans])
  dGDP%<>%group_by(ind)%>%summarise(value=log(unique(GDP)*(Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]),"value"]/
                                                             Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]) & Sinc$variable=="p50p100","value"])),
                                    income=Sinc[Sinc$iso3==unique(dGDP$iso[dGDP$ind==unique(ind)]),"variable"],
                                    variable="dollar",
                                    iso3c=unique(dGDP$iso[dGDP$ind==unique(ind)]),.groups = 'drop_last')
  dGDP$linp<-rep(NA_real_,nrow(dGDP))
  dGDP$linp[which(!is.na(dGDP$ind))]<-dlinpred(dGDP[!is.na(dGDP$ind),],beta,center,fIndies)
  
  return(list(dGDP=dplyr::select(dGDP,c(ind,income,linp)),iGDP=iGDP))
}

Plinpred<-function(Pdens,beta,center,fIndies,notnans){
  
  Plinp<-rep(NA_real_,length(Pdens))
  Plinp[notnans]<-dlinpred(data.frame(variable=rep("Pdens",length(Pdens[notnans])),value=log(Pdens[notnans]+1)),beta,center,fIndies)
  return(Plinp)
  
}

GetLP<-function(ODD,Omega,Params,Sinc,notnans){
  linp<-llinpred(Params[c(unique(ODD@cIndies$iso3))],Omega$beta,Params$center,Params$fIndies)
  # Calculate possible dollar (GDP*income_dist) linear predictor values
  GDP<-GDPlinp(ODD,Sinc,Omega$beta,Params$center,Params$fIndies,notnans)
  # Calculate population density linear predictor values
  Plinp<-Plinpred(ODD@data$Population,
                  Omega$beta,
                  Params$center,
                  Params$fIndies,
                  notnans)
  return(list(linp=linp,dGDP=GDP$dGDP,iGDP=GDP$iGDP,Plinp=Plinp))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Generalised Linear Models for building damage & displacement calcs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Binary regression function
if(Model$BinR=="weibull") {
  BinR<-function(x,zeta) pweibull(x,shape=zeta$k,scale=zeta$lambda)
} else if(Model$BinR=="gompertz") {
  BinR<-function(x,zeta) flexsurv::pgompertz(x,shape=zeta$varrho,rate=zeta$eta)
} else stop("Incorrect binary regression function name, try e.g. 'weibull'")

# Survival model, given the hazard intensity
fDamUnscaled<-function(I,Params,Omega,Ik){
    h_0(I,Params$I0,Omega$theta)* # Expected damage level
       stochastic(Params$Np,Ik,Omega$eps) # Add stochasticity
}

# Baseline hazard function h_0
h_0<-function(I,I0,theta){
  ind<-I>I0
  h<-rep(0,length(I))
  h[ind]<-exp( theta$e*(I[ind]-I0) ) -1
  return(h)
}
# Building damage baseline hazard function hBD_0
hBD<-function(Ab,Population,rho,center){
  exp(-rho$A*(log(Ab)-center$A) - rho$H*(log(Population)-center$H))
}
# Stochastic damage function process
stochastic<-function(n,Ik,eps){
  # return(1+rnorm(n = n,mean = 0, sd = (eps$eps+eps$xi*Ik) ))
  # return(rgammaM(n = n,mu = 1, sig_percent = (eps$eps+eps$xi*Ik) ))
  return(rgammaM(n = n,mu = 1, sig_percent = eps$eps ))
}
# Ik is no longer used: ignore this function
# What is the model for Ik? I'm trying quadratic, but we'll see!
CalcStochasticIk<-function(WID){
  Ik<-WID%>%group_by(iso3)%>%summarise(value=sd(value)^2,.groups = 'drop_last')
  cbind(Ik,data.frame(variable=rep("Ik",nrow(Ik))))
}
# For building damage assessment data
fDamageBuilding<-function(BD,I,Params,Omega,linp,Ik){
  BinR(fDamUnscaled(I,Params,Omega,Ik)*linp,#*hBD(BD$Ab,BD$Population,Omega$rho,Params$center[c("A","H")]),
       Omega$zeta)
}

qualifierDisp<-function(Disp,qualifier,mu) {
  if(qualifier%in%c("total","approximately")) return(Disp)
  else if(qualifier=="more than") {
    return(vapply(Disp,function(Disp) rgammaM(n=1,mu = mu$muplus*Disp,
                                              sig_percent = mu$sigplus),numeric(1)))
  } else if(qualifier=="less than") {
    return(vapply(Disp,function(Disp) rgammaM(n=1,mu = mu$muminus*Disp,
                                              sig_percent = mu$sigminus),numeric(1)))
  } else stop(paste0("qualifier is not recognised - ",qualifier))
}

# Binomial displacement calculator function
rbiny<-function(size,p) rbinom(n = 1,size,p); 
Fbdisp<-function(lPopS,Dprime) mapply(rbiny,lPopS,Dprime)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Log likelihood, posterior and prior distribution calculations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# These high-level priors are to include expert opinion in the 
# model calculations, not in the individual parameters themselves (for this, see priors)
Model$HighLevelPriors<-function(Omega,Model){
  
  if(Model$haz=="EQ"){
    
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 
    Dispfun<-function(I_ij) BinR(Dfun(I_ij)*Dfun(I_ij)*Omega$Lambda$kappa+Omega$Lambda$nu*Dfun(I_ij) + Omega$Lambda$omega,Omega$zeta)
    Damfun<-function(I_ij) BinR(Dfun(I_ij),Omega$zeta)
    
    # Add lower bound priors:
    adder<-sum(500*pweibull(c(Dispfun(4.6),Damfun(4.6)),3,0.001))
    # Middle range priors:
    adder<-adder+sum(150*pweibull(c(Dispfun(6),Damfun(6)),15,0.8) + 150*(1- pweibull(c(Dispfun(6),Damfun(6)),3,0.005)))
    # Upper bound priors:
    adder<-adder+sum(500*(1-pweibull(c(Dispfun(9),Damfun(9)),30,0.85)))
    
    return(-adder)
    
  } else if(Model$haz=="TC"){
    
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 3,theta = Omega$theta) 
    Dispfun<-function(I_ij) BinR(Dfun(I_ij)*Dfun(I_ij)*Omega$Lambda$kappa+Omega$Lambda$nu*Dfun(I_ij) + Omega$Lambda$omega,Omega$zeta)
    Damfun<-function(I_ij) BinR(Dfun(I_ij),Omega$zeta)
    
    # Add lower bound priors:
    adder<-sum(50*pweibull(c(Dispfun(3.05),Damfun(3.05)),3,0.001))
    # Middle range priors:
    # adder<-adder+sum(15*(1- pweibull(c(Dispfun(3.53),Damfun(3.53)),3,0.005)) +15*pweibull(c(Dispfun(3.53),Damfun(3.53)),15,0.8) )
    # Upper bound priors:
    adder<-adder+sum(50*(1-pweibull(c(Dispfun(45),Damfun(45)),30,0.85)))
    
    return(-adder)
    
  }
  
}

# Get the log-likelihood for the displacement data
LL_IDP<-function(Y){
   -log(1+(Y$gmax-Y$predictor)^2/Y$gmax)
  # -((Y$gmax-Y$predictor)^2/Y$gmax)
  # dnorm((log(max(Ystar,1))-log(Y)), mean = 0, sd = log(Y), log = T)
  # dgammaM(Ystar,Y,0.1,log = T)
}

# Preparation for building damage log-likelihood
LL_beta_apply<-function(b,value,BD_params) do.call(BD_params$functions[[value]],as.list(c(x=b,unlist(BD_params$Params[[value]]))))

# Building damage log-likelihood based on building classification in EMS-98 format
LL_BD<-function(b,classified,BD_params){
  
  lls<-array(dim=c(length(b),length(BD_params$Params)))
  for(i in 1:length(BD_params$Params)){
    value<-(names(BD_params$Params))[i]
    lls[,i]<-vapply(b,FUN = LL_beta_apply,FUN.VALUE = numeric(1),
                    value=value,BD_params=BD_params)
  }
  lls%<>%as.data.frame.array();colnames(lls)<-names(BD_params$Params)
  
  # 'Damaged' is a special classification, requiring special treatment 
  if(classified=="Damaged") {
    tmp<-rowSums(lls)
    # Sum of all rows for the classifications that predict at least some damage
    return(log((tmp-lls[["notaffected"]])/tmp))
  }
  
  # Log likelihood is based on relativising against the likelihood to be in other classifications
  out<-log(lls[[classified]]/rowSums(lls))
  # machine level precision ~ -300
  out[is.infinite(out)| out< -300]<--300
  return(out)
  
}

# Log-likelihood for displacement (ODD) objects
LL_ODD<-function(LL,dir,Model,proposed,AlgoParams,expLL=T){
  # Load ODD files
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/ODDobjects_count"),pattern=Model$haz,recursive = T,ignore.case = T)
  for(i in 1:length(ufiles)){
    # Extract the ODD object
    ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects_count/",ufiles[i]))
    # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
    ODDy@fIndies<-Model$fIndies
    ODDy@gmax%<>%as.data.frame.list()
    # Apply DispX
    tLL<-tryCatch(DispX(ODD = ODDy,Omega = proposed,center = Model$center,LL = T,Method = AlgoParams),
                  error=function(e) NA)
    # If all is good, add the LL to the total LL
    if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
    if(expLL) LL<-LL+max(log(mean(exp(tLL),na.rm=T)),AlgoParams$cap,na.rm = T)
    else LL<-LL+max(mean(tLL,na.rm=T),AlgoParams$cap,na.rm = T)
  } 
  
  return(LL)
}

# Log-likelihood for building damage (BD) objects
LL_BD<-function(LL,dir,Model,proposed,AlgoParams,expLL=T){
  # Load BD files
  ufiles<-list.files(path=paste0(dir,"IIDIPUS_Input/BDobjects_count"),pattern=Model$haz,recursive = T,ignore.case = T)
  for(i in 1:length(ufiles)){
    # Extract the BD object
    BDy<-readRDS(paste0(dir,"IIDIPUS_Input/BDobjects_count/",ufiles[i]))
    # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
    BDy@fIndies<-Model$fIndies
    # Apply BDX
    tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,center = Model$center,Method=AlgoParams),
                  error=function(e) NA)
    # If all is good, add the LL to the total LL
    if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
    if(expLL) LL<-LL+max(log(mean(exp(tLL),na.rm=T)),AlgoParams$cap,na.rm = T)
    else LL<-LL+max(mean(tLL,na.rm=T),AlgoParams$cap,na.rm = T)
  }
  
  return(LL)
}

# Bayesian Posterior distribution: This is for a group of ODD objects with observed data
logTarget<-function(dir,Model,proposed,AlgoParams,expLL=T){
  
  HP<-0
  # Apply higher order priors
  if(!is.null(Model$HighLevelPriors)){
    HP<-Model$HighLevelPriors(proposed,Model)
    print(paste0("Higher Level Priors = ",HP))
  }
  
  # Add the log-likelihood values from the ODD (displacement) objects
  LL%<>%LL_ODD(dir,Model,proposed,AlgoParams,expLL=T)
  print(paste0("LL Displacements = ",LL)) ; sLL<-LL
  
  # Add the log-likelihood values from the BD (building damage) objects
  LL%<>%LL_BD(dir,Model,proposed,AlgoParams,expLL=T)
  print(paste0("LL Building Damages = ",LL-sLL))
  
  posterior<-LL+HP
  # Add Bayesian priors
  if(!is.null(Model$priors)){
    posterior<-posterior+sum(Priors(proposed,Model$priors),na.rm=T)
  }
  print(paste0("Posterior = ",posterior))
  
  return(posterior)
  
}






# }
# vals<-seq(0.0001,0.99999,length.out = 50)
# plot(vals,dbeta(vals,Model$BD_params$Params$destroyed$shape1,
#                      Model$BD_params$Params$destroyed$shape2,log = T),
#      ylim=c(-50,5),col="black",ylab="Log Probability (beta)",xlab="Damage %",
#      type="l",lwd=3)
# lines(vals,dbeta(vals,Model$BD_params$Params$severe$shape1,
#                        Model$BD_params$Params$severe$shape2,log = T),col="red",
#       lwd=2,type="b",pch=0)
# lines(vals,dbeta(vals,Model$BD_params$Params$moderate$shape1,
#                   Model$BD_params$Params$moderate$shape2,log = T),col="orange",
#       lwd=2,type="b",pch=1,lty=2)
# lines(vals,dbeta(vals,Model$BD_params$Params$possible$shape1,
#                   Model$BD_params$Params$possible$shape2,log = T),col="green",
#       lwd=2,type="b",pch=2,lty=3)
# lines(vals,dbeta(vals,Model$BD_params$Params$notaffected$shape1,
#                   Model$BD_params$Params$notaffected$shape2,log = T),col="blue",
#       lwd=2,type="b",pch=3,lty=4)
# legend(x=0.45, y=-20, legend=c("Destroyed","Severe","Moderate","Possible","Unaffected"),
#        col=c("black","red","orange","green","blue"),lty=c(1,1,2,3,4),lwd = 2,pch = c(0,0:3))