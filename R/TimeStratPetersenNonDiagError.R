# 2010-04-26 CJS fixed problem with init.logitP when n1=m2=k and you get +infinity which craps out lm()
# 2010-03-03 CJS allowed some logitP[j] top be fixed at arbitrary values (on the logit scale)
#                added definition of storage.class(logitP) to deal with no fixed values where the program bombed
# 2009-12-07 CJS changed etaP to logitP
# 2009-12-05 CJS added title to argument list
# 2009-12-01 CJS added openbugs/winbugs directory to argument list

TimeStratPetersenNonDiagError <- function(title, prefix, time, n1, m2, u2,
                     jump.after=NULL,
                     logitP.cov=rep(1,length(u2)), logitP.fixed=rep(NA,length(u2)),
                     n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
                     tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
                     mu_xiP=-2, tau_xiP=.6666, 
                     tauP.alpha=.001, tauP.beta=.001, 
                     debug=FALSE, debug2=FALSE, openbugs=TRUE, InitialSeed,
                     OPENBUGS.directory, WINBUGS.directory ){
#
#  Fit the smoothed time-Stratified Petersen estimator with NON-Diagonal recoveries 
#  This model allows recoveries outside the stratum of release and error in the smoothed U curve
#
#  Packages Required - must be installed BEFORE calling this functin
#
#    R2WinBugs  - needed to call WinBugs to fit the model
#    coda
#    splines
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususally
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Of the n1 fish released, some are recapturd in this stratum of release (column 1 of m2) or in
#     subsequent strata (subsequent columns of m2). No fish are assumed to be available for capture
#     outside the range of strata considered in the matrix of m2
#  At the same tine, u2 other (unmarked) fish are newly captured in stratum i. 
#     These EXCLUDE recaptures of marked fish. These are the fish that are "expanded"
#     to estimate the population size of fish in stratum i.  
# 
#  Input
#      prefix - prefix for file name for initial plot of U's
#      time- the stratum number
#      n1  - vector of number of fish released in stratum i
#      m2  - matrix of number of fish recovered who were released in stratum i and recovered in stratum j
#      u2  - vector of number of unmarked fish captured in stratum i
#      jump.after - points after which the spline is allowed to jump. Specify as a list of integers in the 
#              range of 1:Nstrata. If jump.after[i]=k, then the spline is split between strata k and k+1
#      logitP.cov - covariates for logit(P)=X beta.logitP.cov
#                 - specify anything you want for fixed logitP's as the covariate values are simply ignored. 
#                 - recommend that you specify 1 for the intercept and 0's for everything else
#      logitP.fixed- values for logitP that are fixed in advance. Use NA if corresponding value is not fixed, 
#                    otherwise specify the logitP value.


#  This routine makes a call to WinBugs to fit the model and then gets back the 
#  coda files for the posteriour distribution.

library("R2WinBUGS")  # Make sure that all the packages needed are available
library("coda")       # used for convergence diagnostics
library("splines")    # used to compute the Bspline design matrix

# Save the WinBugs progam to the model.txt file
#
sink("model.txt")  # NOTE: NO " allowed in model as this confuses the cat command
cat("
model TimeStratPetersenNonDiagError{
# Time Stratified Petersen with NON Diagonal recapture and allowing for error in the smoothed U curve.

# Refer to Bonner (2008) Ph.D. thesis from Simon Fraser University available at
#     http://www.stat.sfu.ca/people/alumni/Theses/Bonner-2008.pdf
# The model is in Appendix B. The discussion of the model is in Chapter 2.

#  Data input:
#      Nstrata.rel - number of strata where fish are releases
#      Nstrata.cap - number of (future strata) where fish are recaptured.
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#                   This is a matrix of size Nstrata.rel x (Nstrata.cap+1)
#                   with entries m2[i,j] = number of fish released in i and recaptured in j
#                   Entries in the last column are the number of fish NEVER recaptured from those
#                   released
#      u2         - number of unmarked fish captured (To be expanded to population).
#      logitP     - the recapture rates. Use NA if these are modelled, otherwise specify the logit(fixed value, e.g. -10 for 0).
#      Nfree.logitP - number of free logitP parameters
#      free.logitP.index - vector of length(Nfree.logitP) for the free logitP parameters
#      logitP.cov - covariates for logitP
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
#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#      xiMu, tauMu  - mean and precision (1/variance) for prior on mean(log travel-times)
#      siSd, tauSd  - mean and precision (1/variance) for prior on sd(log travel times)
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
#
#      muLogTT[i] = mean log(travel time) assuming a log-normal distribution for travel time
#      sdLogTT[i] = sd   log(travel time) assuming a log-normal distribution for travel time
#                        Note that the etasdLogTT=log(sdLogTT) is modelled to keep the sd positive
#       

   ##### Fit the spline for the U's and specify hierarchial model for the logit(P)'s ######
   for(i in 1:Nstrata.cap){
        logUne[i] <- inprod(SplineDesign[i,1:n.bU],bU[1:n.bU])  # spline design matrix * spline coeff 
        etaU[i] ~ dnorm(logUne[i], taueU)              # add random error
        eU[i] <- etaU[i] - logUne[i]
   }
   for(i in 1:Nfree.logitP){   # model the free capture rates using covariates
        mu.logitP[free.logitP.index[i]] <- inprod(logitP.cov[free.logitP.index[i],1:NlogitP.cov], beta.logitP[1:NlogitP.cov]) 
        logitP[free.logitP.index[i]] ~ dnorm(mu.logitP[free.logitP.index[i]],tauP)
   }


   ##### Hyperpriors #####
   ## Mean and sd of log travel-times
   for(i in 1:Nstrata.rel){
     muLogTT[i]    ~ dnorm(xiMu,tauMu)
     etasdLogTT[i] ~ dnorm(xiSd,tauSd)
   }

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

   ## prior on the  Mean and log(sd) of log travel times
   xiMu  ~ dnorm(0,.0625)
   tauMu ~ dgamma(1,.01)
   sigmaMu <- 1/sqrt(tauMu)

   xiSd  ~ dnorm(0,.0625)
   tauSd ~ dgamma(1,.1)
   sigmaSd <- 1/sqrt(tauSd)

   ##### Compute derived parameters #####
   ## Get the sd of the log(travel times)
   for(i in 1:Nstrata.rel){
      log(sdLogTT[i]) <- etasdLogTT[i]
   }

   ## Transition probabilities
   for(i in 1:Nstrata.rel){
     # Probability of transition in 0 days (T<1 days)
     Theta[i,i] <- phi((log(1)-muLogTT[i])/sdLogTT[i])
     for(j in (i+1):Nstrata.cap){
       # Probability of transition in j days (j-1<T<j)
       Theta[i,j] <- phi((log(j-i+1)-muLogTT[i])/sdLogTT[i])- phi((log(j-i)-muLogTT[i])/sdLogTT[i])
     }
     Theta[i,Nstrata.cap+1] <- 1-sum(Theta[i,i:Nstrata.cap])  # fish never seen again
   }

   ##### Likelihood contributions #####
   ## marked fish ##
   for(i in 1:Nstrata.rel){
      # Compute cell probabilities
      for(j in i:Nstrata.cap){
        Pmarked[i,j] <- Theta[i,j] * p[j]
      }
      Pmarked[i,Nstrata.cap+1] <- 1- sum(Pmarked[i,i:Nstrata.cap])

      # Likelihood contribution
      m2[i,i:(Nstrata.cap+1)] ~ dmulti(Pmarked[i,i:(Nstrata.cap+1)],n1[i])
   }
   ## Capture probabilities and run size
   for(j in 1:Nstrata.cap){
      logit(p[j]) <- max(-10,min(10,logitP[j]))    # convert from logit scale; use limits to avoid over/underflow
      U[j]   <- round(exp(etaU[j]))       # convert from log scale
      u2[j] ~ dbin(p[j],U[j])      # capture of newly unmarked fish
   }

   ##### Derived Parameters #####
   Utot <- sum( U[1:Nstrata.cap])          # Total number of unmarked fish
   Ntot <- sum(n1[1:Nstrata.rel]) + Utot  # Total population size including those fish marked and released
} # end of model

", fill=TRUE)
sink()  # End of saving the WinBugs program


# Now to create the initial values, and the data prior to call to WinBugs

Nstrata.rel <- length(n1)
Nstrata.cap <- ncol(m2)-1  # remember last column of m2 has the number of fish NOT recovered


Uguess <- pmax(c((u2[1:Nstrata.rel]+1)*(n1+2)/(apply(m2[,1:Nstrata.cap],1,sum)+1),rep(0,Nstrata.cap-Nstrata.rel)),
          (u2+1)/expit(mu_xiP), na.rm=TRUE)  # try and keep Uguess larger than observed values


# create the B-spline design matrix
# Each set of strata separated at the jump.after[i] points forms a separate spline with a separate basis
# We need to keep track of the breaks as the first two spline coefficients will have a flat
# prior and the others are then related to the previous values.

ext.jump <- c(0, jump.after, Nstrata.cap)  # add the first and last breakpoints to the jump sets
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


# get the logitP covariate matrix ready 
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))

# get the logitP's ready to allow for fixed values
logitP <- as.numeric(logitP.fixed)
storage.mode(logitP) <- "double" # if there are no fixed logits, the default class will be logical which bombs
free.logitP.index <- (1:Nstrata.cap)[ is.na(logitP.fixed)]  # free values are those where NA is specifed
Nfree.logitP <- length(free.logitP.index)


datalist <- list("Nstrata.rel", "Nstrata.cap", "n1", "m2", "u2", 
                 "logitP", "Nfree.logitP", "free.logitP.index",   # those indices that are fixed and free to vary
                 "logitP.cov", "NlogitP.cov",
                 "SplineDesign",
                 "b.flat", "n.b.flat", "b.notflat", "n.b.notflat", "n.bU",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")


# get the initial values for the parameters of the model


Uguess <- pmax(c((u2[1:Nstrata.rel]+1)*(n1+2)/(apply(m2[,1:Nstrata.cap],1,sum)+1),rep(0,Nstrata.cap-Nstrata.rel)),
          (u2+1)/expit(mu_xiP), na.rm=TRUE)  # try and keep Uguess larger than observed values
init.bU   <- lm(log(Uguess) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
if(debug2) {
   cat("compute init.bU \n")
   browser()  # Stop here to examine the spline design matrix function
}

logitPguess <- c(logit((apply(m2[,1:Nstrata.cap],1,sum)+1)/(n1+1)),rep(mu_xiP,Nstrata.cap-Nstrata.rel))
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
                "Ntot", "Utot", "logUne", "etaU", "U",
                 "muLogTT", "sdLogTT")      #  mean and sd of log(travel times)
if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2))) {parameters <- c(parameters,"u2")}
                 
init.vals <- function(){
   init.logitP <- c(logit((apply(m2[,1:Nstrata.cap],1,sum)+1)/(n1+1)),rep(mu_xiP,Nstrata.cap-Nstrata.rel))         # initial capture rates based on observed recaptures
   init.logitP <- pmin(10,pmax(-10,init.logitP))
   init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
   init.logitP[!is.na(logitP.fixed)]<- NA        # no need to initialize the fixed values
   init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
   init.beta.logitP[is.na(init.beta.logitP)] <- 0 
   init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a 
                                             # vector in the init files etc.
   init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

   init.bU   <- lm(log(Uguess) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
   init.eU   <- as.vector(log(Uguess)-SplineDesign%*%init.bU)  # error terms set as differ between obs and pred
   init.etaU <- log(Uguess)

   # variance of spline difference
   sigmaU <- sd( init.bU[b.notflat]-2*init.bU[b.notflat-1]+init.bU[b.notflat-2], na.rm=TRUE)
   init.tauU <- 1/sigmaU^2

   # variance of error in the U over and above the spline fit
   sigmaeU <- sd(init.eU, na.rm=TRUE)
   init.taueU <- 1/sigmaeU^2

   # initialize the u2 where missing
   init.u2    <- u2
   init.u2[ is.na(u2)] <- 100
   init.u2[!is.na(u2)] <- NA

   # mean log(travel time) and sd log(travel time)
   init.muLogTT <- as.vector(log(pmax(
                   (m2[,1:Nstrata.cap] %*% 1:Nstrata.cap)/(1+apply(m2[,1:Nstrata.cap],1,sum)) - 1:Nstrata.rel,
                    1, na.rm=TRUE)))
   init.muLogTT <- as.vector(init.muLogTT)
   init.etasdLogTT <- log(rep(.5,Nstrata.rel))  # note that log (sd(log travel time)) is being modelled
   
   list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP, 
        bU=init.bU,  tauU=init.tauU, taueU=init.taueU, etaU=init.etaU, 
        muLogTT=init.muLogTT, etasdLogTT=init.etasdLogTT)
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
   cat("Initial values generated in ",inits, "\n")
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
   BRugs::samplesCoda("*", stem=paste(working.directory,"/",sep=""))  # write out coda files
   cat("Coda file created \n")

   cat("\n\n*** Finished OpenBugs ***\n\n")
   results
   } else {
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


} # end of TimeStratPetersenNonDiagError function
