## 2020-11-07 CJS Allow user to specify prior for beta parameters for covariates on logitP
# 2018-12-06 CJS converted initial plot to ggplot2
# 2018-11-25 CJS Remove all references to OpenBugs
# 2014-09-01 CJS converstion to JAGS
#                - no model name
#                - C(,20) -> T(,20)
#                - dflat() to dnorm(0, 1E-6)
#                - added u2copy to improve mixing based on Matt S. suggestion
# 2011-05-15 CJS limited etaU to 20 to prevent overflow in binomial computations
# 2011-03-09 CJS added prior to muTT (mean.muTT and sd.muTT) with defaults same a previously
# 2011-02-17 CJS limited initial Pguess to between .01 and .99 to avoid taking logit of 0 or 1
# 2011-02-17 CJS fixed initial values for theta; added as.matrix to deal with case of Delta.max=1
# 2011-01-24 SB  added call to run.windows.openbugs and run.windows.winbugs
# 2010-11-25 CJS output on progress of burnin and post-burnin phases
# 2010-04-25 CJS fixed problems of init.logitP=+infinity if n1=m2=k which crapped out the lm() call
# 2010-03-11 Added fixed values of logitP[j] to be provided by user.
# 2010-03-02 SJB Created file.

#' @keywords internal
#' @importFrom stats lm spline sd


TimeStratPetersenNonDiagErrorNP <- function(title,
                                            prefix,
                                            time,
                                            n1,
                                            m2,
                                            u2,
                                            jump.after=NULL,
                                            logitP.cov=as.matrix(rep(1,length(u2))),
                                            logitP.fixed=rep(NA,length(u2)),
                                            n.chains=3,
                                            n.iter=200000,
                                            n.burnin=100000,
                                            n.sims=2000,
                                            tauU.alpha=1,
                                            tauU.beta=.05,
                                            taueU.alpha=1,
                                            taueU.beta=.05,
                                            Delta.max,
                                            mean.muTT=rep(0,Delta.max),
                                            sd.muTT=rep(sqrt(.666),Delta.max),
                                            tauTT.alpha=.1,
                                            tauTT.beta=.1,
                                            prior.beta.logitP.mean = c(logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),rep(0,  ncol(as.matrix(logitP.cov))-1)),
                                            prior.beta.logitP.sd   = c(2,                                           rep(10, ncol(as.matrix(logitP.cov))-1)), 
                                            tauP.alpha=.001,
                                            tauP.beta=.001,
                                            debug=FALSE,
                                            debug2=FALSE,
                                            InitialSeed,
                                            save.output.to.files=TRUE){

#   browser()
set.seed(InitialSeed)  # set prior to initial value computations

#
#  Fit the smoothed time-Stratified Petersen estimator with NON-Diagonal recoveries.
#  This model allows recoveries outside the stratum of release and error in the smoothed U curve.
#  The travel time model is based on the continuation ratio and makes no parametric assumptions.
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
#      logitP.fixed - values for logitP that are fixed in advance. Use NA if corresponding value is not fixed,
#                    otherwise specify the logitP value.


#  This routine makes a call to the MCMC sample to fit the model and then gets back the
#  coda files for the posteriour distribution.

## Set working directory to current directory (we should allow users to select this)
working.directory <- getwd()

## Define paths for the model, data, and initial value files
model.file <- file.path(working.directory, "model.txt")
data.file <- file.path(working.directory,"data.txt")
init.files <- file.path(working.directory,
                       paste("inits", 1:n.chains,".txt", sep = ""))

# Save the Bugs progam to the model.txt file
#
sink(model.file)  # NOTE: NO " allowed in model as this confuses the cat command
cat("

model {

# Time Stratified Petersen with NON Diagonal recapture and allowing for error in the smoothed U curve.
# Non-parametric estimateion of travel times for marked individuals.
#

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
#      prior.beta.logitP.mean, prior.beta.logitP.sd  - parameters for prior of coefficient of covariates for logitP
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#
#  Parameters of the model are:
#      p[i]
#       logitP[i]  = logit(p[i]) = logitP.cov*beta.logitP
#         The beta coefficients have a prior that is N(mean= prior.beta.logitP.mean, sd= prior.beta.logitP.sd)
#      U[i]
#       etaU[i]  = log(U[i])
#         which comes from spline with parameters bU[1... Knots+q]
#         + error term eU[i]
#
#      muTT[j] = mean(logit(delta[i,i+j-1])), j=1,...,Delta.max
#      sdTT = stats::sd(logit(delta[i,i+j-1])), j=1,....,Delta.max
#      delta[i,i+j-1]=Theta[i,i+j-1]/(1-Theta[i,i]-...-Theta[i,i+j-2])
#

   ##### Fit the spline for the U's and specify hierarchial model for the logit(P)'s ######
   for(i in 1:(Nstrata.cap)){
        ## Model for U's
        logUne[i] <- inprod(SplineDesign[i,1:n.bU],bU[1:n.bU])  # spline design matrix * spline coeff
        etaU[i] ~ dnorm(logUne[i], taueU)T(,20)              # add random error but keep from getting too large
        eU[i] <- etaU[i] - logUne[i]
   }

   for(i in 1:Nfree.logitP){   # model the free capture rates using covariates
        mu.logitP[free.logitP.index[i]] <- inprod(logitP.cov[free.logitP.index[i],1:NlogitP.cov], beta.logitP[1:NlogitP.cov])
        ## logitP[free.logitP.index[i]] ~ dnorm(mu.logitP[free.logitP.index[i]],tauP)
        mu.epsilon[free.logitP.index[i]] <- mu.logitP[free.logitP.index[i]] - log(u2copy[free.logitP.index[i]] + 1) + etaU[free.logitP.index[i]]
        epsilon[free.logitP.index[i]] ~ dnorm(mu.epsilon[free.logitP.index[i]],tauP)

        logitP[free.logitP.index[i]] <- log(u2copy[free.logitP.index[i]] + 1) - etaU[free.logitP.index[i]] + epsilon[free.logitP.index[i]]
   }
   # define the last epsilon (including the extra needed for m2)
   for(i in Extra.strata.cap){
      epsilon[ Nstrata.cap + i] <- 0  # forces definition of epsilon1 ...epsilon[Nstrata.cap -> Extra.strata.cap]
   }

   ##### Priors and hyperpriors #####
   ## Transition probabilities -- continuation ratio model
   for(i in 1:Nstrata.rel){
     ## delta[i,j] is the probability that a marked fish released on day i passes the second trap
     ## on day i+j-1 given that it does not pass the on days i,...,i+j-2. r[i,j]=logit(delta[i,j])
     ## is assumed to have a normal distribution with mean muTT[j] and precision tauTT.
     r[i,1] ~ dnorm(muTT[1],tauTT)

     logit(Theta[i,1] ) <- r[i,1]

     for(j in 2:Delta.max){
       r[i,j] ~ dnorm(muTT[j],tauTT)
       logit(delta[i,j]) <- r[i,j]
       Theta[i,j] <- delta[i,j] * (1 - sum(Theta[i,1:(j-1)]))
     }
     Theta[i,Delta.max+1] <- 1- sum(Theta[i,1:Delta.max])
   }

#  derived parameters on actual movement probabilities
   logit(movep[1]) <- muTT[1]
   for(j in 2:Delta.max){
      movep[j] <- ilogit(muTT[j]) *(1- sum(movep[1:(j-1)]))
   }
   movep[Delta.max+1] <- 1- sum(movep[1:Delta.max])

#  prior on the movement rates. These are specified using the make.muTT.prior function
   for(j in 1:Delta.max){
     muTT[j] ~ dnorm(mean.muTT[j],tau.muTT[j])
   }
   tauTT~ dgamma(tauTT.alpha,tauTT.beta)
   sdTT <- 1/sqrt(tauTT)

   ## Run size - flat priors
   for(i in 1:n.b.flat){
      bU[b.flat[i]] ~ dnorm(0, 1E-6)
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

   ## Capture probabilities covariates
   for(i in 1:NlogitP.cov){
      beta.logitP[i] ~ dnorm(prior.beta.logitP.mean[i], 1/prior.beta.logitP.sd[i]^2)  # rest of beta terms are normal 0 and a large variance
   }
   beta.logitP[NlogitP.cov+1] ~ dnorm(0, .01) # dummy so that covariates of length 1 function properly
   tauP ~ dgamma(tauP.alpha,tauP.beta)
   sigmaP <- 1/sqrt(tauP)

   ##### Likelihood contributions #####
   ## marked fish ##
   for(i in 1:Nstrata.rel){
      # Compute cell probabilities
      for(j in 1:(Delta.max+1)){
        Pmarked[i,j] <- Theta[i,j] * p[i+j-1]
      }
      Pmarked[i,Delta.max+2] <- 1- sum(Pmarked[i,1:(Delta.max+1)])

      # Likelihood contribution
      m2[i,1:(Delta.max+2)] ~ dmulti(Pmarked[i,],n1[i])
   }

   ## Capture probabilities and run size
   for(j in 1:(Nstrata.cap + Extra.strata.cap)){
      logit(p[j]) <- logitP[j]       # convert from logit scale
   }
   for(j in 1:Nstrata.cap){
      U[j]   <- round(exp(etaU[j]))       # convert from log scale
      u2[j] ~ dbin(p[j],U[j])      # capture of newly unmarked fish
   }

   ##### Derived Parameters #####
   Utot <- sum( U[1:Nstrata.cap])          # Total number of unmarked fish
   Ntot <- sum(n1[1:Nstrata.rel]) + Utot  # Total population size including those fish marked and released
} # end of model
", fill=TRUE)

sink()  # End of saving the Bugs program

# Now to create the initial values, and the data prior to call to MCMC sampler

Nstrata.rel <- length(n1)
Nstrata.cap <- length(u2)

## Count extra columns that will have to be added to account for Delta.max
Extra.strata.cap <- max(0,Nstrata.rel + ncol(m2) - Nstrata.cap -1)

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
logitP <- c(as.numeric(logitP.fixed),rep(-10,Extra.strata.cap))
storage.mode(logitP) <- "double"  # force the storage class to be correct if there are no fixed values
free.logitP.index <- (1:Nstrata.cap)[ is.na(logitP.fixed)]  # free values are those where NA is specifed
Nfree.logitP <- length(free.logitP.index)

tau.muTT <- 1/sd.muTT**2  # convert from sd to precision = 1/variance

# make a copy of u2 to improve mixing
u2copy <- exp(stats::spline(x = 1:length(u2), y = log(u2+1), xout = 1:length(u2))$y)-1 # on log scale to avoid negative values
u2copy <- pmax(0,round(u2copy)) # round to integers

datalist <- list("Nstrata.rel", "Nstrata.cap","Extra.strata.cap",
                 "Delta.max","n1", "m2", "u2", "u2copy",
                 "logitP", "Nfree.logitP", "free.logitP.index",
                 "logitP.cov", "NlogitP.cov",
                 "SplineDesign",
                 "b.flat", "n.b.flat", "b.notflat", "n.b.notflat", "n.bU",
                 "mean.muTT", "tau.muTT", # priors on muTT
                 "tauTT.alpha","tauTT.beta",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "prior.beta.logitP.mean", "prior.beta.logitP.sd", 
                 "tauP.alpha", "tauP.beta")


## Generate the initial values for the parameters of the model

## 1) U and spline coefficients
Uguess <- pmax((u2+1)/expit(prior.beta.logitP.mean[1]),1)  # try and keep Uguess larger than observed values
Uguess[which(is.na(Uguess))] <- mean(Uguess,na.rm=TRUE)

init.bU   <- stats::lm(log(Uguess) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients

if(debug2) {
   cat("compute init.bU \n")
   browser()  # Stop here to examine the spline design matrix function
}

## 2) Capture probabilities
logitPguess <- c(logit(pmin(.99,pmax(.01,(apply(m2[,1:(Delta.max+1)],1,sum)+1)/(n1+1)))),
                 rep(prior.beta.logitP.mean[1],Nstrata.cap-Nstrata.rel))
#browser()
init.beta.logitP <- as.vector(stats::lm( logitPguess ~ logitP.cov-1)$coefficients)
if(debug2) {
   cat(" obtained initial values of beta.logitP\n")
   browser()
}


# create an initial plot of the fit
plot.data <- data.frame(time=time, 
                        logUguess=log(Uguess[1:Nstrata.cap]),
                        spline=SplineDesign %*% init.bU, stringsAsFactors=FALSE)
init.plot <- ggplot(data=plot.data, aes_(x=~time, y=~logUguess))+
  ggtitle(title, subtitle="Initial spline fit to estimated log U[i]")+
  geom_point()+
  geom_line(aes_(y=~spline))+
  xlab("Stratum")+ylab("log(U[i])")+
  scale_x_continuous(breaks=seq(min(plot.data$time, na.rm=TRUE),max(plot.data$time,na.rm=TRUE),2))

if(save.output.to.files)ggsave(init.plot, filename=paste(prefix,"-initialU.pdf",sep=""), height=4, width=6, units="in")
#results$plots$plot.init <- init.plot  # do this after running the MCMC chain (see end of function)


parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU", "tauU", "sigmaU",
                "eU", "taueU", "sigmaeU",
                "Ntot", "Utot", "logUne", "etaU", "U",
                 "muTT","sdTT","Theta","movep")

if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2))) {parameters <- c(parameters,"u2")}

init.vals <- genInitVals("TSPNDENP",
                         n1,m2,u2,
                         Delta.max=Delta.max,
                         logitP.cov=logitP.cov,
                         logitP.fixed=logitP.fixed,
                         SplineDesign=SplineDesign,
                         n.chains=n.chains)
#browser()
## Generate data list
data.list <- list()
for(i in 1:length(datalist)){
  data.list[[length(data.list)+1]] <-get(datalist[[i]])
}
names(data.list) <- as.list(datalist)

# Set up for the call to the MCMC sampler

results <- run.MCMC(modelFile=model.file,
                        dataFile=data.file,
                        dataList=data.list,
                        initFiles=init.files,
                        initVals=init.vals,
                        parameters=parameters,
                        nChains=n.chains,
                        nIter=n.iter,
                        nBurnin=n.burnin,
                        nSims=n.sims,
                        overRelax=FALSE,
                        initialSeed=InitialSeed,
                        working.directory=working.directory,
                        debug=debug)

results$plots$plot.init <- init.plot  # save initial plot to results object

return(results)
}
