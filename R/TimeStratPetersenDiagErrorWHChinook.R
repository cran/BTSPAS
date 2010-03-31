# 2009-12-05 CJS added title to argument list
# 2009-12-01 CJS Added open/win bugs path names to argument list

TimeStratPetersenDiagErrorWHChinook <-
       function(title, prefix, time, n1, m2, u2.A, u2.N,
                hatch.after=NULL, clip.frac.H=.25, 
                logitP.cov,
                n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
                tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
                mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)), 
                tau_xiP=1/var(logit((m2+.5)/(n1+1)), na.rm=TRUE), 
                tauP.alpha=.001, tauP.beta=.001, 
                debug=FALSE, debug2=FALSE, openbugs=TRUE,
                InitialSeed,
                OPENBUGS.directory, WINBUGS.directory){

#
#  Fit the smoothed time-Stratified Petersen estimator with Diagonal recoveries (i.e. no recoveries
#  outside stratum of release), error in the smoothed U curve, and separating wild vs hatchery stocks
#
#  Packages Required - must be installed BEFORE calling this functin
#
#    R2WinBugs  - needed to call WinBugs to fit the model
#
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususally
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Of the n1 fish released, m2 fish are recaptured in the same stratum (e.g. week) of release.
#     There is a related function that allows fish to be recaptured in subsequent weeks.
#
#  At the same tine, u2.A (adipose fin clipped) and u2.N (not clipped)
#     other (unmarked) fish are newly captured in stratum i. 
#     These EXCLUDE recaptures of marked fish. These are the fish that are "expanded"
#     to estimate the population size of fish in stratum i.
#     All wild fish are NOT ad-clipped.
#     Only a fraction of hatchery fish are ad-clipped. It is assumed that the fraction of ad-clipped
#     hatchery fish is constant over the life of the run.       
# 
#  Input
#      prefix - prefix for file name for initial plot of U's
#      time  - the stratum number
#      n1    - vector of number of fish released in stratum i
#      m2    - vector of number of fish recovered in stratum i (EXCLUDING recaps)
#      u2.A  - vector of number of adclipped unmarked fish captured in stratum i
#      u2.N  - vector of number of non-clipped unmarked fish captured in stratum i
#      hatch.after - point AFTER which the hatchery fish are released.
#      clip.frac.H - what fraction of hatchery fish are clipped
#      logitP.cov - covariates for logit(P)=X beta.logitP


#  This routine makes a call to WinBugs to fit the model and then gets back the 
#  coda files for the posteriour distribution.

library("R2WinBUGS")  # Make sure that all the packages needed are available
library("coda")       # used for convergence diagnostics
library("splines")    # used to compute the Bspline design matrix

# Save the WinBugs progam to the model.txt file
#
sink("model.txt")  # NOTE: NO " allowed in model as this confuses the cat command
cat("
model TimeStratPetersenDiagErrorWHChinook{
# Time Stratified Petersen with Diagonal recapture (no spillover in subsequent weeks or marked fish)
#    and allowing for error in the smoothed U curve with separation of wild and hatchery fish
# Each of the wild and hatchery populations are fit using a SINGLE spline curve as this should be flexible 
#    enough to capture the individual behaviours

#  Data input:
#      Nstrata - number of strata
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#      u2.A       - number of adclipped unmarked fish captured (must be hatchery fish).
#      u2.N       - number of non-clipped unmarked fish captured (wild + hatchery fish)
#      clip.frac.H- what fraction of hatchery fish are clipped
#      logitP.cov   - covariates for logitP
#      NlogitP.cov  - number of logitP covariates
#      SplineDesign.W- wildfish spline design matrix of size [Nstrata, maxelement of n.b.notflat.W]
#      SplineDesign.H- hatchery spline design matrix of size [Nstrata, maxelement of n.b.notflat.H]
#                   This is set up prior to the call.
#      b.flat.W   - vector of strata indices where the prior for the b's will be flat for wild     fish
#      b.flat.H   - vector of strata indices where the prior for the b's will be flat for hatchery fish
#                 this is normally the first two weeks of each spline segment
#      n.b.flat.W - number of b coefficients that have a flat prior - wild     fish
#      n.b.flat.H - number of b coefficients that have a flat prior - hatchery fish
#      b.notflat.W - vector of strata indices where difference in coefficients is modelled - wild     fish
#      b.notflat.H - vector of strata indices where difference in coefficients is modelled - hatchery fish
#      n.b.notflat.W - number of b coefficients that do not have a flat prior - wild     fish
#      n.b.notflat.H - number of b coefficients that do not have a flat prior - hatchery fish
#      tauU.alpha, tauU.beta   - parameters for prior on tauU
#      taueU.alpha, taueU.beta - parameters for prior on taueU
#      mu_xiP, tau_xiP  - parameters for prior on mean logit(P)'s [The intercept term]
#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#      clip.frac.H    - what fraction of hatchery fish are clipped (KNOWN in advance)
#
#  Parameters of the model are:
#      p[i]
#       logitP[i]  = logit(p[i]) = logitP.cov*beta.logitP
#         The first beta.logitP has a prior from N(xiP, tauP)
#            and xiP and tauP are currently set internally
#         The remaining beta's are assigned a wider prior N(mu=0, var=1000).
#      U.W[i] - number of unmarked wild     fish passing stratam i in population
#      U.H[i] - number of unmarked hatchery fish passing stratum i in population
#       etaU.W[i]  = log(U.W[i])
#       etaU.H[i]  = log(U.H[i])
#         which comes from spline with parameters bU.W[1... Knots+q] or bU.H[1... knots+q]
#         + error term eU.W[i] or eu.H[i]        

   ##### Fit the spline for wildfish - this covers the entire experiment ######
   for(i in 1:Nstrata){
        logUne.W[i] <- inprod(SplineDesign.W[i,1:n.bU.W],bU.W[1:n.bU.W])  # spline design matrix * spline coeff 
        etaU.W[i] ~ dnorm(logUne.W[i], taueU)              # add random error
        eU.W[i] <- etaU.W[i] - logUne.W[i]
   }
   ##### Fit the spline for hatchery fish - these fish only enter AFTER hatch.after ######
   for(i in (hatch.after+1):Nstrata){
        logUne.H[i] <- inprod(SplineDesign.H[i,1:n.bU.H],bU.H[1:n.bU.H])  # spline design matrix * spline coeff 
        etaU.H[i] ~ dnorm(logUne.H[i], taueU)              # add random error
        eU.H[i] <- etaU.H[i] - logUne.H[i]
   }

   ##### Model the capture probabilities #####
   for(i in 1:Nstrata){
        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov]) 
        logitP[i] ~ dnorm(mu.logitP[i],tauP)
   }


   ##### Hyperpriors #####
   ## Run size - wild and hatchery fish - flat priors
   for(i in 1:n.b.flat.W){
      bU.W[b.flat.W[i]] ~ dflat()
   }
   for(i in 1:n.b.flat.H){
      bU.H[b.flat.H[i]] ~ dflat()
   }

   ## Run size - priors on the difference for wild and hatchery fish
   for(i in 1:n.b.notflat.W){
      xiU.W[b.notflat.W[i]] <- 2*bU.W[b.notflat.W[i]-1] - bU.W[b.notflat.W[i]-2]
      bU.W [b.notflat.W[i]] ~ dnorm(xiU.W[b.notflat.W[i]],tauU)
   }
   for(i in 1:n.b.notflat.H){
      xiU.H[b.notflat.H[i]] <- 2*bU.H[b.notflat.H[i]-1] - bU.H[b.notflat.H[i]-2]
      bU.H [b.notflat.H[i]] ~ dnorm(xiU.H[b.notflat.H[i]],tauU)
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
   ## Number of marked fish recovered ##
   for(i in 1:Nstrata){
      logit(p[i]) <- logitP[i]       # convert from logit scale
      m2[i] ~ dbin(p[i],n1[i])     # recovery of marked fish
   }

   ## captures of wild (unclipped fish) - these are the only fish available upto (and including) hatch.after
   for(i in 1:hatch.after){
      U.W[i] <- round(exp(etaU.W[i]))       # convert from log scale
      u2.N[i] ~ dbin(p[i],U.W[i])
   }

   ## captures of hatchery (clipped fish) - these can only occur AFTER hatch.after
   for(i in (hatch.after+1):Nstrata){
      U.W[i] <- round(exp(etaU.W[i]))      # convert from log scale
      U.H[i] <- round(exp(etaU.H[i]))      # convert from log scale
      U.clip[i] ~ dbin(clip.frac.H, U.H[i])
      #u2.A[i] ~ dbin(p[i], U.clip[i])
      p.temp[i] <- p[i]*clip.frac.H
      u2.A[i] ~ dbin(p.temp[i], U.H[i]) # must be hatchery and clipped
   }

   ## captures of wild+hatchery unclipped fish - these can only occur AFTER hatch.after
   for(i in (hatch.after+1):Nstrata){
      U.noclip[i] <- U.W[i] + U.H[i] - U.clip[i]
      u2.N[i] ~ dbin(p[i], U.noclip[i])
   } 

   ##### Derived Parameters #####
   Utot.W <- sum( U.W[1:Nstrata])              # Total number of unmarked fish - wild
   Utot.H <- sum( U.H[(hatch.after+1):Nstrata])# Total number of unmarked fish - hatchery
   Utot   <- Utot.W + Utot.H                   # Grand total number of fish
}  # end of model

", fill=TRUE)
sink()  # End of saving the WinBugs program


datalist <- list("Nstrata", "n1", "m2", "u2.A", "u2.N", "hatch.after", "clip.frac.H",
                 "logitP.cov", "NlogitP.cov",
                 "SplineDesign.W",
                 "b.flat.W", "n.b.flat.W", "b.notflat.W", "n.b.notflat.W", "n.bU.W",
                 "SplineDesign.H",
                 "b.flat.H", "n.b.flat.H", "b.notflat.H", "n.b.notflat.H", "n.bU.H",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")

parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU.W", "bU.H", "tauU", "sigmaU", 
                "eU.W", "eU.H", "taueU", "sigmaeU", 
                "Utot.W", "Utot.H", "Utot", "logUne.W", "logUne.H",
                "etaU.W", "etaU.H", "U.W", "U.H")
if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2.A))) {parameters <- c(parameters,"u2.A")}
if( any(is.na(u2.N))) {parameters <- c(parameters,"u2.N")}


# Now to create the initial values, and the data prior to call to WinBugs

Nstrata <- length(n1)

# Estimate number of wild and hatchery fish based on clip rate
u2.H <- u2.A/clip.frac.H  # only a portion of the hatchery fish are clipped
u2.W <- pmax(u2.N - u2.H*(1-clip.frac.H),0) # subtract the questimated number of hatchery fish
u2.H[is.na(u2.H)] <- 1  # in case of missing values
u2.W[is.na(u2.W)] <- 1  # in case of missing values

avg.P <- sum(m2,na.rm=TRUE)/sum(n1, na.rM=TRUE)
Uguess.W <- pmax((u2.W+1)*(n1+2)/(m2+1), u2.W/avg.P, na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess.H <- pmax((u2.H+1)*(n1+2)/(m2+1), u2.H/avg.P, na.rm=TRUE)
Uguess.H[1:hatch.after] <- 0   # no hatchery fish prior to release from hatchery

# create the B-spline design matrix for wild and hatchery fish
# The design matrix for hatchery fish will still have rows corresponding to entries PRIOR to 
#    the hatchery release but these are never used in the winbugs fitting routines
# There is a separate (single) spline for hatchery and wild fish with NO breakpoints
# The first two coefficient have a flat prior and the rest of the coefficients are modelled using
#    differences between the succesive coefficients

# Wild fish. This covers the entire experiment.
SplineDegree <- 3           # Degree of spline between occasions
knots <- seq(4,Nstrata,4)/(Nstrata+1) # a knot roughly every 4th stratum
SplineDesign.W <- bs((1:Nstrata)/(Nstrata+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.W      <- c(1,2)
b.notflat.W   <- 3:(ncol(SplineDesign.W))
n.b.flat.W    <- length(b.flat.W)
n.b.notflat.W <- length(b.notflat.W)
n.bU.W        <- n.b.flat.W + n.b.notflat.W
init.bU.W   <- lm(log(Uguess.W+1) ~ SplineDesign.W-1)$coefficients  # initial values for spline coefficients

# hatchery fish. Notice they can only enter AFTER hatch.after, The spline design matrix still has rows
# of zero for 1:hatch.after to make it easier in WinBugs
SplineDegree <- 3           # Degree of spline between occasions
knots <- (seq((hatch.after+4),Nstrata-1,4)-hatch.after)/(Nstrata-hatch.after+1) # a knot roughly every 4th stratum
SplineDesign.H <- bs((1:(Nstrata-hatch.after))/(Nstrata-hatch.after+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.H      <- c(1,2)
b.notflat.H   <- 3:(ncol(SplineDesign.H))
n.b.flat.H    <- length(b.flat.H)
n.b.notflat.H <- length(b.notflat.H)
n.bU.H        <- n.b.flat.H + n.b.notflat.H
init.bU.H   <- lm(log(Uguess.H[(hatch.after+1):Nstrata]+1) ~ SplineDesign.H-1)$coefficients  # initial values for spline coefficients
# patch up the initial rows of the spline design matrix
SplineDesign.H <- rbind(matrix(0,nrow=hatch.after, ncol=ncol(SplineDesign.H)), SplineDesign.H)

# create an initial plot of the fit to the number of unmarked fish
pdf(file=paste(prefix,"-initialU.pdf",sep=""))
plot(time, log(Uguess.H+1), 
    main=paste(title,"\nInitial spline fit to estimated U.W[i] and U.H[i]"),
    ylab="log(U[i])", xlab='Stratum', pch="H")  # initial points on log scale.
points(time, log(Uguess.W+1), pch="W")
lines(time, SplineDesign.W %*% init.bU.W)  # add smoothed spline through points
lines(time, SplineDesign.H %*% init.bU.H)  # add smoothed spline through points
dev.off()

#browser()

# get the logitP=logit(P) covariate matrix ready 
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))


#  initial values for the parameters of the model
                
init.vals <- function(){
   # Initial values for the probability of capture
   init.logitP <- logit((m2+1)/(n1+2))         # initial capture rates based on observed recaptures
   init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
   init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
   init.beta.logitP[init.beta.logitP=NA] <- 0 
   init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a 
                                             # vector in the init files etc.
   init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

   # inital values for the spline coefficients
   init.bU.W   <- lm(log(Uguess.W+1) ~ SplineDesign.W-1)$coefficients  # initial values for spline coefficients
   init.bU.H   <- lm(log(Uguess.H[(hatch.after+1):Nstrata]+1) ~ SplineDesign.H[(hatch.after+1):Nstrata,]-1)$coefficients  # initial values for spline coefficients

   init.eU.W   <- as.vector(log(Uguess.W+1)-SplineDesign.W%*%init.bU.W)  # error terms set as differ between obs and pred
   init.etaU.W <- log(Uguess.W+1)
   init.eU.H   <- as.vector(log(Uguess.H+1)-SplineDesign.H%*%init.bU.H)  # error terms set as differ between obs and pred
   init.etaU.H <- log(Uguess.H+1)
   init.etaU.H[1:hatch.after] <- NA  # these are never used.

   # variance of spline difference (use only the wild fish to initialize)
   sigmaU <- sd( init.bU.W[b.notflat.W]-2*init.bU.W[b.notflat.W-1]+init.bU.W[b.notflat.W-2], na.rm=TRUE)
   init.tauU <- 1/sigmaU^2

   # variance of error in the U' over and above the spline fit (use only the wild fish to initialize)
   sigmaeU <- sd(init.eU.W, na.rm=TRUE)
   init.taueU <- 1/sigmaeU^2

   # initialize the u2.A and u2.N where missing
   init.u2.A    <- u2.A
   init.u2.A[ is.na(u2.A)] <- 100
   init.u2.A[!is.na(u2.A)] <- NA

   init.u2.N    <- u2.N
   init.u2.N[ is.na(u2.N)] <- 100
   init.u2.N[!is.na(u2.N)] <- NA

   list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP, 
        bU.W=init.bU.W, bU.H=init.bU.H,  tauU=init.tauU, taueU=init.taueU, 
        etaU.W=init.etaU.W, etaU.H=init.etaU.H)
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
   cat("\nInitial values loaded into model\n")

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


} # end of TimeStratPetersenDiagError function
