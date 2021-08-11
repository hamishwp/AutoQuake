##########################################################################
##########################################################################
###################### Machine Learning Methodologies ####################
### This file contains the Markov Chain Monte Carlo algorithm          ###
### that is used to optimise the model parameters, given historical    ###
### events and observation data, including building damage evaluations ###
###                                                                    ###
### I chose to use the Generalised Adaptative Metropolis-Hastings      ###
### algorithm with Global Adaptative Scaling                           ###
### see DOI 10.1007/s11222-008-9110-y                                  ###
### By C. Andrieu & J. Thoms, Stat Comput (2008) 18: 343–373           ###
### 'A tutorial on adaptive MCMC'                                      ###
### or otherwise http://drops.dagstuhl.de/opus/volltexte/2010/2813     ###
### C.L. Muller, (ETH) 'Exploring the common concepts of adaptive      ###
### MCMC and covariance matrix adaptation schemes'                     ###
##########################################################################
##########################################################################
##########################################################################

# Methodology parameters required
AlgoParams<-list(Np=30, # Number of Monte Carlo particles
                 cores=12, # Number of parallelised threads per event
                 itermax=500, # How many iterations do we want?
                 cap=-300, # if log values are too low, then log(mean(exp(LL)))=-Inf
                 GreedyStart=20, # How sure are we of the initial covariance matrix for accepted parameters? (Larger -> more confident)
                 Pstar=0.234, # Adaptive metropolis acceptance rate
                 gamzy0=1, # How quickly do the rejected parameters start having an influence on the covariance? (like GreedyStart) 
                 epsilon=2, # Do we still want values at larger numbers of iterations to have an influence on the covariance?
                 minVar=1e-4 # Prevent certain parameters from being too sure of themselves
                 )

# Metropolis-Hastings proposal distribution, given old values and covariance matrix
multvarNormProp <- function(xt, propPars){
  # purpose : A multivariate Gaussian random walk proposal for Met-Hastings
  #           MCMC
  # inputs  : xt       - The value of the chain at the previous time step 
  #           propPars - The correlation structure of the proposal
  return(array(mvtnorm::rmvnorm(1, mean=xt, sigma=propPars),dimnames = list(names(xt))))
}

# Generate the Adaptive Metropolis Global Scaling Factor discount factor
GenerateGamzy<-function(AlgoParams){
  AlgoParams$gamzy0/(1:(AlgoParams$itermax+AlgoParams$cores))^
    (seq(from=1/(1+AlgoParams$epsilon),to=1,length.out=(AlgoParams$itermax+AlgoParams$cores)))
}

# Check that the proposed initial values, model and methodology parameters
checkLTargs<-function(Model,iVals,AlgoParams){
  if(length(unlist(Model$links))!=length(unlist(iVals$x0))) stop("Mismatching link functions and initial values of Omega space")
  if(length(unlist(iVals$x0))!=nrow(iVals$COV) | length(unlist(iVals$x0))!=ncol(iVals$COV)) stop("Mismatching initial values and initial covariance matrix")
}

# Optimisation of the model parameters to minimise the posterior distribution
# (proper developers would split this, but I need to see all the goods!)
Algorithm<-function(dir,Model,iVals,AlgoParams){
  # purpose : Generalised Adaptative Metropolis-Hastings Algorithm with Global Adaptative Scaling
  # Details : 'a tutorial on adaptive MCMC', Andrieu & Thoms, 2008
  #           Look at algorithm 4.
  #           Uses the multivariate Gaussian distribution as proposal, 
  #           whereby the proposal covariance is updated at every iteration.
  
  # Check no mistakes have been made in the model, methodology and initial values
  checkLTargs(Model,iVals,AlgoParams)
  
  ###### Initialisations ######
  Lgsf<-0
  xPrev<-propMu<-unlist(iVals$x0)
  if(!is.null(Model$links)) iVals$x0%<>%unlist()%>%Proposed2Physical(Model)
  propCOV<-iVals$COV
  n <- length(xPrev) + 1
  output <- matrix(NA, nrow=AlgoParams$itermax, ncol=n)
  xNew <- rep(NA,n-1)
  lTargNew<-alpha<-c()
  # Create file names with unique names for storage reasons
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  # Generate Adaptive Metropolis Global Scaling Factor iterative vector
  gamzy <- GenerateGamzy(AlgoParams)
  ##############################
  
  # Find first log-target value using initial values
  output[1, ] <- c(logTarget(dir,Model,iVals$x0,AlgoParams[c("Np","cores")]),unlist(xPrev))
  print(output[1,])
  # Start the iterations!
  it <- 2
  while (it <= AlgoParams$itermax){
  
    # Parameter proposal
    xNew<-multvarNormProp(xt=xPrev, propPars=exp(2*Lgsf)*propCOV/sum(propCOV*t(propCOV)))
    # Convert parameters to physical/useable values
    if(!is.null(Model$links)) xProp<-xNew%>%Proposed2Physical(Model)
    # Calculate log-target value
    lTargNew <- tryCatch(logTarget(dir,Model,xProp,AlgoParams[c("Np","cores")]), error=function(e) NA)
    print(lTargNew)
    # Check if we have a NaN
    if(is.na(lTargNew)|is.infinite(lTargNew)) {
      output[it,] <- output[it-1,]
      it <- it + 1
      next
    }
    # Prepare for acceptance
    lTargOld <- output[it-1, 1]
    u <- runif(1)
    # Acceptance probability
    alpha <- min(c(exp(lTargNew - lTargOld),1))
    # Metropolis Acceptance Algorithm
    if (alpha>=u) { # Accepted!
      output[it,] <- c(lTargNew, xNew)
    } else {  # Rejected
      output[it,] <- output[it-1,]
      if(it<=AlgoParams$GreedyStart) {it <- it + 1; next}
    }
    # Store this for next time!
    xPrev<-xNew
    # Global Scaling Factor (GSF), mean & covariance update
    Lgsf <- Lgsf + gamzy[it]*(alpha-AlgoParams$Pstar)
    print(Lgsf)
    propMu <- propMu + gamzy[it]*(xNew - propMu)
    propCOV <- propCOV + gamzy[it]*((xNew - propMu)%*%t(xNew - propMu) - propCOV)
    propCOV[is.na(propCOV)]<-0
    
    print(" ")
    
    it <- it + 1
  }
  
  return(output)
  
}
