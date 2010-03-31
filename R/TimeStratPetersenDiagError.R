# 2009-12-05 CJS added title to argument list
# 2009-12-01 CJS (added WinBugs/OpenBugs directory to the argument list

TimeStratPetersenDiagError <- function(title, prefix, time, n1, m2, u2,
                     jump.after=NULL,
                     logitP.cov,
                     n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
                     tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
                     mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),
                     tau_xiP=1/var(logit((m2+.5)/(n1+1)),na.rm=TRUE),
                     tauP.alpha=.001, tauP.beta=.001, 
                     debug=FALSE, debug2=FALSE, openbugs=TRUE,
                     InitialSeed, 
                     OPENBUGS.directory, WINBUGS.directory){

#
#  Fit the smoothed time-Stratified Petersen estimator with Diagonal recoveries (i.e. no recoveries
#  outside stratum of release) and error in the smoothed U curve
#
#  Packages Required - must be installed BEFORE calling this functin
#
#    R2WinBugs  - needed to call WinBugs to fit the model
#    BRugs
#    Coda
#    actuar
#    splines
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususall
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Of the n1 fish released, m2 fish are recaptured in the same stratum (e.g. week) of release.
#     There is a related function that allows fish to be recaptured in subsequent weeks.
#  At the same tine, u2 other (unmarked) fish are newly captured in stratum i. 
#     These EXCLUDE recaptures of marked fish. These are the fish that are "expanded"
#     to estimate the population size of fish in stratum i.  
# 
#  Input
#      prefix - prefix for file name for initial plot of U's
#      time- the stratum number
#      n1  - vector of number of fish released in stratum i
#      m2  - vector of number of fish recovered in stratum i (EXCLUDING recaps)
#      u2  - vector of number of unmarked fish captured in stratum i
#      jump.after - points after which the spline is allowed to jump. Specify as a list of integers in the 
#              range of 1:Nstrata. If jump.after[i]=k, then the spline is split between strata k and k+1
#      logitP.cov - covariates for logit(P)=X beta.logitP


#  This routine makes a call to WinBugs to fit the model and then gets back the 
#  coda files for the posteriour distribution.

library("R2WinBUGS")  # Make sure that all the packages needed are available
library("coda")       # used for convergence diagnostics
library("splines")    # used to compute the Bspline design matrix
library("BRugs")

# Save the WinBugs progam to the model.txt file
#
sink("model.txt")  # NOTE: NO " allowed in model as this confuses the cat command
cat("
model TimeStratPetersenDiagError{
# Time Stratified Petersen with Diagonal recapture (no spillover in subsequent weeks or marked fish)
#    and allowing for error in the smoothed U curve.

# Refer to Bonner (2008) Ph.D. thesis from Simon Fraser University available at
#     http://www.stat.sfu.ca/people/alumni/Theses/Bonner-2008.pdf
# The model is in Appendix B. The discussion of the model is in Chapter 2.

#  Data input:
#      Nstrata - number of strata
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#      u2         - number of unmarked fish captured (To be expanded to population).
#      logitP.cov   - covariates for logitP
#      NlogitP.cov  - number of logitP covariates
#      SplineDesign- spline design matrix of size [Nstrata, maxelement of n.b.notflat]
#                   This is set up prior to the call.
#      b.flat   - vector of strata indices where the prior for the b's will be flat.
#                 this is normally the first two of each spline segment
#      n.b.flat - number of b coefficients that have a flat prior
#      b.notflat- vector of strata indices where difference in coefficients is modelled
#      n.b.notflat- number of b coefficients that do not have a flat prior
#      tauU.alpha, tauU.beta - parameters for prior on tauU
#      taueU.alpha, taueU.beta - parameters for prior on taueU
#      mu_xiP, tau_xiP  - parameters for prior on mean logit(P)'s [The intercept term]
#                         mu_xiP = mean; tau_xiP = 1/variance
#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#
#  Parameters of the model are:
#      p[i]
#       logitP[i]  = logit(p[i]) = logitP.cov*beta.logitP
#         The first beta.logitP has a prior from N(xiP, tauP)
#            and xiP and tauP are currently set internally
#         The remaining beta's are assigned a wider prior N(mu=0, var=1000).
#      U[i]
#       etaU[i]  = log(U[i])
#         which comes from spline with parameters bU[1... Knots+q]
#         + error term eU[i]        

   ##### Fit the spline and specify hierarchial model for the logit(P)'s ######
   for(i in 1:Nstrata){
        logUne[i] <- inprod(SplineDesign[i,1:n.bU],bU[1:n.bU])  # spline design matrix * spline coeff 
        etaU[i] ~ dnorm(logUne[i], taueU)              # add random error
        eU[i] <- etaU[i] - logUne[i]
        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov]) 
        logitP[i] ~ dnorm(mu.logitP[i],tauP)
   }
   ##### Hyperpriors #####
   ## Run size - flat priors
   for(i in 1:n.b.flat){
      bU[b.flat[i]] ~ dflat()
   }
   ## Run size - priors on the difference
   for(i in 1:n.b.notflat){
      xiU[b.notflat[i]] <- 2*bU[b.notflat[i]-1] - bU[b.notflat[i]-2]
      bU [b.notflat[i]] ~ dnorm(xiU[b.notflat[i]],tauU)
   }
   tauU ~ dgamma(tauU.alpha,tauU.beta)  # Notice reduction from .0005 (in thesis) to .05
   sigmaU <- 1/sqrt(tauU)
   taueU ~ dgamma(taueU.alpha,taueU.beta) # dgamma(100,.05) # Notice reduction from .0005 (in thesis) to .05
   sigmaeU <- 1/sqrt(taueU)

   ## Capture probabilities. The logit(p[i]) are n(logitP.cov*beta.logitP.cov, sigmaP**2)
   beta.logitP[1] ~ dnorm(mu_xiP,tau_xiP) # first term is usually an intercept
   for(i in 2:NlogitP.cov){
      beta.logitP[i] ~ dnorm(0, .01)   # rest of beta terms are normal 0 and a large variance
   }
   beta.logitP[NlogitP.cov+1] ~ dnorm(0, .01) # dummy so that covariates of length 1 function properly
   tauP ~ dgamma(tauP.alpha,tauP.beta)
   sigmaP <- 1/sqrt(tauP)

   ##### Likelihood contributions #####
   for(i in 1:Nstrata){
      logit(p[i]) <- logitP[i]       # convert from logit scale
      U[i]   <- round(exp(etaU[i]))       # convert from log scale
      m2[i] ~ dbin(p[i],n1[i])     # recovery of marked fish
      u2[i] ~ dbin(p[i],U [i])      # capture of newly unmarked fish
   }

   ##### Derived Parameters #####
   Utot <- sum( U[1:Nstrata])          # Total number of unmarked fish
   Ntot <- sum(n1[1:Nstrata]) + Utot  # Total population size including those fish marked and released
} # end of model

", fill=TRUE)
sink()  # End of saving the WinBugs program


# Now to create the initial values, and the data prior to call to WinBugs

Nstrata <- length(n1)
avgP <- sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)
Uguess <- pmax((u2+1)*(n1+2)/(m2+1), u2/avgP, na.rm=TRUE)  # try and keep Uguess larger than observed values


# create the B-spline design matrix
# Each set of strata separated at the jump.after[i] points forms a separate spline with a separate basis
# We need to keep track of the breaks as the first two spline coefficients will have a flat
# prior and the others are then related to the previous values.

ext.jump <- c(0, jump.after, Nstrata)  # add the first and last breakpoints to the jump sets
SplineDesign <- matrix(0, nrow=0, ncol=0)
SplineDegree <- 3           # Degree of spline between occasions
b.flat <- NULL              # index of spline coefficients with a flat prior distribution -first two of each segment
b.notflat <- NULL           # index of spline coefficients where difference is modelled
all.knots <- NULL    
for (i in 1:(length(ext.jump)-1)){
  nstrata.in.set <- ext.jump[i+1]-ext.jump[i]
  if(nstrata.in.set > 7)
    { knots   <- seq(5,nstrata.in.set-1,4)/(nstrata.in.set+1) # a knot roughly every 4th stratum
    } else{
      knots   <- .5       # a knot roughly every 4th stratum
    }
  all.knots <- c(all.knots, knots)
  # compute the design matrix for this set of strata
  z <- bs((1:nstrata.in.set)/(nstrata.in.set+1), knots=knots, degree=SplineDegree, 
             intercept=TRUE, Boundary.knots=c(0,1))
  # first two elements of b coeffients have a flat prior
  b.flat <- c(b.flat, ncol(SplineDesign)+(1:2))
  b.notflat <- c(b.notflat, ncol(SplineDesign)+3:(ncol(z)))
  # add to the full design matrix which is block diagonal
  SplineDesign <- cbind(SplineDesign, matrix(0, nrow=nrow(SplineDesign), ncol=ncol(z)))
  SplineDesign <- rbind(SplineDesign,
                         cbind( matrix(0,nrow=nrow(z),ncol=ncol(SplineDesign)-ncol(z)), z)  )
  } # end of for loop
n.b.flat <- length(b.flat)
n.b.notflat <- length(b.notflat)
n.bU <- n.b.flat + n.b.notflat


# get the logitP=logit(P) covariate matrix ready 
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))

datalist <- list("Nstrata", "n1", "m2", "u2", "logitP.cov", "NlogitP.cov",
                 "SplineDesign",
                 "b.flat", "n.b.flat", "b.notflat", "n.b.notflat", "n.bU",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")


# get the initial values for the parameters of the model
  
avgP <- sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE) 
Uguess <- pmax((u2+1)*(n1+2)/(m2+1), u2/avgP, 100, na.rm=TRUE)  # try and keep Uguess larger than observed values
init.bU   <- lm(log(Uguess+1) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
if(debug2) {
   cat("compute init.bU \n")
   browser()  # Stop here to examine the spline design matrix function
}

logitPguess <- logit( (m2+1)/(n1+1))
init.beta.logitP <- as.vector(lm( logitPguess ~ logitP.cov-1)$coefficients)
if(debug2) {
   cat(" obtained initial values of beta.logitP\n")
   browser()
}



# create an initial plot of the fit
pdf(file=paste(prefix,"-initialU.pdf",sep=""))
plot(time, log(Uguess), 
    main=paste(title,"\nInitial spline fit to estimated U[i]"),
    ylab="log(U[i])", xlab='Stratum')  # initial points on log scale.
lines(time, SplineDesign %*% init.bU)  # add smoothed spline through points
dev.off()


parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU", "tauU", "sigmaU", 
                "eU", "taueU", "sigmaeU", 
                "Ntot", "Utot", "logUne", "etaU", "U")
if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2))) {parameters <- c(parameters,"u2")}
                 
init.vals <- function(){
   init.logitP <- logit((m2+1)/(n1+2))         # initial capture rates based on observed recaptures
   init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
   init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
   init.beta.logitP[init.beta.logitP=NA] <- 0 
   init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a 
                                             # vector in the init files etc.
   init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

   init.bU   <- lm(log(Uguess+1) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
   init.eU   <- as.vector(log(Uguess)-SplineDesign%*%init.bU)  # error terms set as differ between obs and pred
   init.etaU <- log(Uguess)

   # variance of spline difference
   sigmaU <- sd( init.bU[b.notflat]-2*init.bU[b.notflat-1]+init.bU[b.notflat-2], na.rm=TRUE)
   init.tauU <- 1/sigmaU^2

   # variance of error in the U' over and above the spline fit
   sigmaeU <- sd(init.eU, na.rm=TRUE)
   init.taueU <- 1/sigmaeU^2

   # initialize the u2 where missing
   init.u2    <- u2
   init.u2[ is.na(u2)] <- 100
   init.u2[!is.na(u2)] <- NA
   
   list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP, 
        bU=init.bU,  tauU=init.tauU, taueU=init.taueU, etaU=init.etaU)
}

#browser()


# set up for the call to WinBugs or OpenBugs

working.directory <- getwd()          # store all data and init files in the current directory
if(openbugs){
   cat("\n\n*** Start of call to OpenBugs \n")
   bugs.directory = OPENBUGS.directory 
   # the following call sequence is mainly based on the "openbugs" function in the R2WinBugs package
   modelFile <- "model.txt"  # this was written in the current directory earlier
   numChains <- n.chains
   nBurnin   <- n.burnin
   nIterPostBurnin <- n.iter - n.burnin
   nThin     <- round(nIterPostBurnin/n.sims)  # want certain # of final samples so decide upon the thining rate
   over.relax <- FALSE
   parameters.to.save <- c(parameters, "deviance")  # always get the deviance
   parametersToSave <- parameters.to.save

   cat("OpenBugs files created in:", working.directory, "\n")
   BRugs::modelCheck(modelFile)    # check the model

   # get the data and write to the data.txt file in the current directory
   datafileName = file.path(working.directory, "data.txt")
   data.list <- list()
   for(i in 1:length(datalist)){
      data.list[[length(data.list)+1]] <-get(datalist[[i]])
   }
   names(data.list) <- as.list(datalist)
   temp <- BRugs::bugsData(data.list, fileName = datafileName, digits = 5)
   cat("Data files saved in ", temp, "\n")

   BRugs::modelData(datafileName)
   cat("Data loaded into model\n")
   BRugs::modelCompile(numChains)

   # set the random number generator seed
   BRugs::modelSetSeed(InitialSeed)
   cat("Random seed initialized with :", InitialSeed, "\n")

   # generate the files to save the initial values
   initfileNames <- file.path(working.directory, paste("inits", 1:numChains,".txt", sep = ""))
   inits <- BRugs::bugsInits(inits = init.vals, numChains = numChains, fileName=initfileNames)
   cat("Initial values generated in ",paste(inits,"\n"), "\n")
   BRugs::modelInits(inits)
   BRugs::modelGenInits()     # generate the initial values for any uninitialized variables
   cat("Initial values loaded into model\n")

   # now to generate the burnin sample
   cat("Burnin sampling has been started for ", nBurnin, " iterations.... \n")
   flush.console()
   BRugs::modelUpdate(nBurnin, overRelax = over.relax)
   cat("Burnin sampling completed \n")

   # generate the non-burnin samples
   BRugs::dicSet()      # turn on DIC computations
   on.exit(BRugs::dicClear(), add = TRUE) 
   cat("DIC collection set \n")  
   BRugs::samplesSet(parametersToSave)
   cat("Nodes to monitor set\n")
   cat("Starting sampling after burnin for ", n.chains," chain each with  a further ", 
        nIterPostBurnin, " iterations. \n A thining rate of ", nThin, 
        "will give about ", round(nIterPostBurnin/nThin), " posterior values in each chain... \n")
   BRugs::modelUpdate(round(nIterPostBurnin/nThin), thin=nThin, overRelax = over.relax) # we do the thining on the fly
   cat("Finished sampling after burnin and thining \n")
   FinalSeed <- BRugs::modelGetSeed(i=1)
   cat("Random seed ended with :", FinalSeed, "\n")


 
   # Now to extract the sampled values and create the bugs array
   cat("Extracting the sampled values\n")
   params <- BRugs::samplesMonitors("*")
   samples <- sapply(params, BRugs::samplesSample)
   n.saved.per.chain <- nrow(samples)/numChains
   samples.array <- array(samples, c(n.saved.per.chain, numChains, ncol(samples)))    
   dimnames(samples.array)[[3]] <- dimnames(samples)[[2]]
   DICOutput <- BRugs::dicStats()

   # save the information. The simulation has already been thinned, so no need to thin again
   results<- as.bugs.array(sims.array = samples.array, 
        model.file = modelFile, program = "OpenBUGS", DIC = TRUE, 
        DICOutput = DICOutput, n.iter = n.iter, n.burnin = n.burnin, 
        n.thin = nThin)
   results$Seed.initial <- InitialSeed
   results$Seed.final   <- FinalSeed
   cat("Final dimension of saved simulation output is ", dim(results$sims.array), "\n")

   # save the information to the coda files
   BRugs::samplesCoda("*", stem=paste(working.directory,"/",sep=""))  # write out code files
   cat("Coda file created \n")

   cat("\n\n*** Finished OpenBugs ***\n\n")
   results
 }

else {
   bugs.directory = WINBUGS.directory  
   results <- bugs( 
      datalist,   # notice this is a list of external variables 
      inits=init.vals, # function to generate initial values for each chain
      parameters,
      model.file="model.txt",
      n.chains=n.chains,
      n.iter=n.iter,    # includes burn in
      n.burnin=n.burnin,
      n.sims=n.sims,      # (approx) number of final simulations to save
      bugs.directory=bugs.directory,
      working.directory=working.directory,
      debug=debug     )
   results
   } # end of else clause


} # end of TimeStratPetersenDiagError function

