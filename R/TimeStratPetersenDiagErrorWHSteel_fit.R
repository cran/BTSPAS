# 2024-05-09 CJS TestIfPool only if at least 2 valid release groups
# 2021-10-24 CJS Added trunc.logitP to fix plotting problems with extreme values of logitP
# 2020-12-15 CJS Remove sampfrac from code
# 2020-11-07 CJS Allowed user to specify prior for beta coefficient for logitP
# 2018-12-19 CJS deprecated use of sampling fraction
# 2018-12-06 CJS converted report to textConnection
# 2018-12-06 CJS converted initial plot to ggplots
# 2018-12-02 CJS converted trace plots to ggplots
# 2018-12-01 CJS converted posterior plots to ggplot2
# 2018-11-30 CJS converted acf to ggplot2
# 2018-11-29 CJS fixed issue with printing of large results got cutoff
# 2018-11-28 CJS remove reference to OpenBugs
# 2015-06-10 CJS convert gof plot to ggplot. Bug fix
# 2014-09-01 CJS conversion to jags
# 2012-08-30 CJS fixed problem with NAs in any() and all() in error checking
# 2011-06-13 CJS added p-values to results
# 2010-11-25 CJS pretty printing of final population estimates
# 2010-09-06 CJS forced input vectors to be vectors
# 2010-08-06 CJS added trace plots of logitP and logU
# 2010-08-03 CJS added version/date to final results
# 2010-03-12 CJS added n.sims etc to calling arguments so users can set
# 2009-12-01 CJS added openbugs/winbugs directory to argument list; added some basic error checking of arguments



#' Wrapper (*_fit) and function to call the Time Stratified Petersen Estimator
#' with Diagonal Entries and separating Wild from Hatchery Steelhead function.
#' 
#' Takes the number of marked fish released, the number of recaptures, and the
#' number of unmarked fish and uses Bayesian methods to fit a fit a spline
#' through the population numbers and a hierarchical model for the trap
#' efficiencies over time.  The output is written to files and an MCMC object
#' is also created with samples from the posterior.
#' 
#' Normally, data is passed to the wrapper which then calls the fitting
#' function.
#' 
#' 
#' @aliases TimeStratPetersenDiagErrorWHSteel_fit 
#' @template title
#' @template prefix
#' @template time
#' @template n1
#' @param m2 A numeric vector of the number of marked fish from n1 that are
#' recaptured in each time stratum. All recaptures take place within the
#' stratum of release. Use the \code{\link{TimeStratPetersenNonDiagError_fit}}
#' function for cases where recaptures take place outside the stratum of
#' release.
#' @param u2.W.YoY A numeric vector of the number of unmarked wild Young-of-Year
#' fish captured in each stratum.
#' @param u2.W.1 A numeric vector of the number of unmarked wild age 1+ fish
#' captured in each stratum.
#' @param u2.H.1 A numeric vector of the number of unmarked hatchery age 1+ fish
#' (i.e. adipose fin clipped) captured in each stratum.
#' @template sampfrac
#' @param hatch.after A numeric vector with elements belonging to \code{time}.
#' At which point do hatchery fish arrive? They arrive in the immediate stratum
#' AFTER these entries.
#' @template bad.n1
#' @template bad.m2 
#' @param bad.u2.W.YoY A numeric vector with elements belonging to \code{time}.
#' In some cases, something goes wrong in the stratum, and the number of wild
#' unmarked Young-of-Year fish should be ignored.
#' @param bad.u2.W.1 A numeric vector with elements belonging to \code{time}.
#' In some cases, something goes wrong in the stratum, and the number of wild
#' unmarked age 1+ fish should be ignored.
#' @param bad.u2.H.1 A numeric vector with elements belonging to \code{time}.
#' In some cases, something goes wrong in the stratum, and the number of
#' hatchery unmarked (but adipose fin clipped) age 1+ fish should be ignored.
#' @template logitP.cov
#' @template mcmc-parms
#' @template tauU.alpha.beta
#' @template taueU.alpha.beta 
#' @template prior.beta.logitP.mean.sd
#' @template tauP.alpha.beta
#' @template run.prob
#' @template debug 
#' @template InitialSeed
#' @template save.output.to.files 
#' @template trunc.logitP 
#' 
#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' @template author
#' @template references
#' @keywords ~models ~smooth
#' @examples
#'  
#' ##---- See the vignettes for example on how to use this package.
#' 
#' @export TimeStratPetersenDiagErrorWHSteel_fit
#' @importFrom stats runif var sd

TimeStratPetersenDiagErrorWHSteel_fit <-
  function( title="TSPDE-WHSteel", prefix="TSPDE-WHSteel-", 
           time, n1, m2, u2.W.YoY, u2.W.1, u2.H.1, sampfrac=rep(1,length(u2.W.YoY)), 
           hatch.after=NULL, 
           bad.n1=c(), bad.m2=c(), bad.u2.W.YoY=c(), bad.u2.W.1=c(), bad.u2.H.1=c(),
           logitP.cov=as.matrix(rep(1,length(n1))),
           n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
           tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
           prior.beta.logitP.mean = c(logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),rep(0,  ncol(as.matrix(logitP.cov))-1)),
           prior.beta.logitP.sd   = c(stats::sd(logit((m2+.5)/(n1+1)),na.rm=TRUE),        rep(10, ncol(as.matrix(logitP.cov))-1)), 
           tauP.alpha=.001, tauP.beta=.001,
           run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
           debug=FALSE, debug2=FALSE,
           InitialSeed=ceiling(stats::runif(1,min=0, 1000000)),
           save.output.to.files=TRUE,
           trunc.logitP=15) {
# Fit a Time Stratified Petersen model with diagonal entries and with smoothing on U allowing for random error,
# covariates for the the capture probabilities, and separating the wild vs hatchery fish for STEELHEAD releases
# The steelhead are nice because 100% of hatchery fish are adipose fin clipped and no wild fish are adipose fin clipped
# The "diagonal entries" implies that no marked fish are recaptured outside the (time) stratum of release
#
    version <- '2024-05-09'
    options(width=200)

# Input parameters are
#    prefix - prefix used for files created with the analysis results
#             this should be in standard Window's format, eg. JC-2002-ST-TSPDE
#             to which is appended various suffixes for plots etc
#    time   - vector of stratum numbers. For example, 9:38 would indicate that the
#             Trinity River system sampled weeks 9 to 38. If some values are omitted
#             e.g. time=10 not present, this indicates sampling did not take place this
#             week. The data are expanded and interpolation for the missing week takes place
#    n1, m2   - the input data consisting of fish marked and released and then recaptured.
#                 The n1 and m2 are used to calibrate the trap
#    u2.W.YoY - number of wild YoY fish (no clips)
#    u2.W.1   - number of wild age 1+ fish (no clips)
#    u2.H.1   - number of hatchery age 1+ fish (ad fin clipped). 100% of hatchery production is fin clipped
#    sampfrac - DEPRECATED. DO NOT USE ANYMORE. 
#    hatch.after - julian week AFTER which hatchery fish are released 
#    bad.m2  - list of julian numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to JAGS
#    bad.u2.W.YoY
#    bad.u2.W.1
#    bad.u2.H.1 - list of julian weeks where the values of u2.* are suspect. 
#                 These are set to NA prior to the call to JAGS
#    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
#               for the first element of the covariance matrix and 0 for the rest of the covariates.
#               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
#               if you have any covariates. The default, if not specified, is a constant (the mean logit)
#    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
#    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline 
#    prior.beta.logitP.mean, prior.beta.logitP.sd   - parameters for the prior on mean logit(P)'s [The intercept term]
#                              The other covariates are assigned priors of a mean of 0 and a sd of 30
#    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
#    run.prob  - percentiles of run timing wanted 
#    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data
#             and JAGS will run and stop waiting for your to exit and complete

# force input vectors to be vectors
time      <- as.vector(time)
n1        <- as.vector(n1)
m2        <- as.vector(m2)
u2.W.YoY  <- as.vector(u2.W.YoY)
u2.W.1    <- as.vector(u2.W.1)
u2.H.1    <- as.vector(u2.H.1)
sampfrac  <- as.vector(sampfrac)

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time all match
if(stats::var(c(length(n1),length(m2),length(u2.W.YoY),length(u2.W.1),length(u2.H.1), length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of n1, m2, u2.W.YoY, u2.W.1, u2.H.1, sampfrac, time must all be equal. They are:",
        length(n1)," ",length(m2)," ",length(u2.W.YoY)," ",length(u2.W.1)," ",length(u2.H.1)," ", length(sampfrac)," ",length(time),"\n")
   return()}
if(!is.numeric(n1)){
   cat("***** ERROR ***** n1 must be numeric. You have:",
        paste(n1,collapse=", "),"\n")
   return()} 
if(any(is.na(n1))){
  cat("***** ERROR ***** All values of n1 must not be missing. You have: ",
        paste(n1,collapse=", "),"\n")
   return()}
if(any(n1 < 0, na.rm=TRUE)){
  cat("***** ERROR ***** All values of n1 must be non-negative. You have: ",
        paste(n1,collapse=", "),"\n")
   return()}

if(length(logitP.cov) %% length(n1) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of n1 etc They are:",
        length(n1)," ",length(logitP.cov)," ",dim(logitP.cov),"\n")
   return()}
#  2. Check that m2<= n1
if(any(m2>n1, na.rm=TRUE)){
   cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:",
       paste(n1,collapse=","),"\n m2:",
       paste(m2,collapse=","),"\n")
   return()}

#  3. Elements of bad.m2, bad.u2.W.YoY, bad.u2.W.1, bad.u2.H.1, and hatch.after must belong to time
if(!all(bad.m2 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",
       paste(bad.m2,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,  collapse=","), "\n")
   return()}
if(!all(bad.u2.W.YoY %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.W.YoY must be elements of strata identifiers. You entered \n bad.u2.W.YoY:",
       paste(bad.u2.W.YoY,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,        collapse=","), "\n")
   return()}
if(!all(bad.u2.W.1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.W.1 must be elements of strata identifiers. You entered \n bad.u2.W.1:",
       paste(bad.u2.W.1,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,      collapse=","), "\n")
   return()}
if(!all(bad.u2.H.1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.H.1 must be elements of strata identifiers. You entered \n bad.u2.H.1:",
       paste(bad.u2.H.1,collapse=","),"\n Strata identifiers are \n time:",
       paste(time      ,collapse=","), "\n")
   return()}
if(!all(hatch.after %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** hatch.after must be elements of strata identifiers. You entered \n hatch.after:",
       paste(hatch.after,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,       collapse=","), "\n")
   return()}

# Check that that the prior.beta.logitP.mean and prior.beta.logitP.sd length=number of columns of covariates
logitP.cov <- as.matrix(logitP.cov)
if(!is.vector(prior.beta.logitP.mean) | !is.vector(prior.beta.logitP.sd)){
   stop("prior.beta.logitP.mean and prior.beta.logitP.sd must be vectors")
}
if(!is.numeric(prior.beta.logitP.mean) | !is.numeric(prior.beta.logitP.sd)){
   stop("prior.beta.logitP.mean and prior.beta.logitP.sd must be numeric")
}
if(length(prior.beta.logitP.mean) != ncol(logitP.cov) | length(prior.beta.logitP.sd) != ncol(logitP.cov)){
   stop("prior.beta.logitP.mean and prior.beta.logitP.sd must be same length as number columns in covariate matrix")
}

# Deprecation of sampling fraction.
if(any(sampfrac != 1)){
  cat("***** ERROR ***** Sampling fraction is deprecated for any values other than 1. DO NOT USE ANYMORE. ")
  return()
}


   results.filename <- paste(prefix,"-results.txt",sep="")   

stdout <- vector('character')
report <- textConnection('stdout', 'wr', local = TRUE)
sink(report)

cat(paste("Time Stratified Petersen with Diagonal recaptures, error in smoothed U, separating wild and hatchery fish, STEELHEAD ONLY - ", date()))
cat("\nVersion: ", version)

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** \n")
temp<- cbind(time, n1, m2, u2.W.YoY, u2.W.1, u2.H.1, logitP.cov)
colnames(temp)<- c('time', 'n1','m2','u2.W.YoY', 'u2.W.1', 'u2.H.1', paste("logitPcov[", 1:ncol(as.matrix(logitP.cov)),"]",sep="") )
print(temp) 
cat("\n\n")
cat("Hatchery fish are released AFTER strata: ", hatch.after,"\n\n")
cat("The following strata had m2       set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.W.YoY set to missing: ", 
     if(length(bad.u2.W.YoY)>0){bad.u2.W.YoY} else {" NONE"}, "\n")
cat("The following strata had u2.W.1   set to missing: ", 
     if(length(bad.u2.W.1)>0){bad.u2.W.1} else {" NONE"}, "\n")
cat("The following strata had u2.H.1   set to missing: ", 
     if(length(bad.u2.H.1)>0){bad.u2.H.1} else {" NONE"}, "\n")


# Pooled Petersen estimator over ALL of the data including when no releases take place, bad m2, bad.u2.* values.
cat("\n\n*** Pooled Petersen Estimates based on pooling over ALL strata ***\n\n")
cat("W.YoY  Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",sum(u2.W.YoY, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.W.YoY, na.rm=TRUE))
cat("W.YoY  Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("W.YoY  Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

cat("W.1    Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=", sum(u2.W.1, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.W.1, na.rm=TRUE))
cat("W.1    Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("W.1    Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

cat("H.1    Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",sum(u2.H.1, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.H.1, na.rm=TRUE))
cat("H.1    Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("H.1    Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")


# Obtain the Pooled Petersen estimator for each group excluding bad.m2, bad.u2.* values
select <- (n1>0) & (!is.na(n1)) & (!is.na(m2)) & (!is.na(u2.W.YoY)) & (!is.na(u2.W.1)) & (!is.na(u2.H.1))
cat("\n\n*** Pooled Petersen Estimate excluding bad m2, u2.* values adjusting for sampling fractions ***\n\n")
cat("The following strata are excluded because n1=0 or NA values in m2, u2.W.YoY, u2.W.1, u2.H.1 :", time[!select],"\n\n")

temp.n1       <- n1      [select]
temp.m2       <- m2      [select]
temp.u2.W.YoY <- u2.W.YoY[select]
temp.u2.W.1   <- u2.W.1  [select]
temp.u2.H.1   <- u2.H.1  [select]

cat("W.YoY  Total n1=", sum(temp.n1, na.rm=TRUE),";  m2=",sum(temp.m2, na.rm=TRUE),";  u2=",sum(temp.u2.W.YoY, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(temp.n1, na.rm=TRUE), sum(temp.m2, na.rm=TRUE), sum(temp.u2.W.YoY, na.rm=TRUE))
cat("W.YoY  Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("W.YoY  Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

cat("W.1    Total n1=", sum(temp.n1, na.rm=TRUE),";  m2=",sum(temp.m2, na.rm=TRUE),";  u2=", sum(temp.u2.W.1, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(temp.n1, na.rm=TRUE), sum(temp.m2, na.rm=TRUE), sum(temp.u2.W.1, na.rm=TRUE))
cat("W.1    Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("W.1    Est U(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

cat("H.1    Total n1=", sum(temp.n1, na.rm=TRUE),";  m2=",sum(temp.m2, na.rm=TRUE),";  u2=",sum(temp.u2.H.1, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(temp.n1, na.rm=TRUE), sum(temp.m2, na.rm=TRUE), sum(temp.u2.H.1, na.rm=TRUE))
cat("H.1    Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("H.1    Est U(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")


# Obtain the Pooled Petersen estimator after fixup of bad.m2, bad.u2.* values after adjusting for sampling fractions
temp.m2 <- m2
index.bad.m2 <- as.vector((1:length(time)) %*% outer(time,bad.m2,"=="))
temp.m2[index.bad.m2] <- NA

temp.u2.W.YoY      <- u2.W.YoY
index.bad.u2.W.YoY <- as.vector((1:length(time)) %*% outer(time,bad.u2.W.YoY,"=="))
temp.u2.W.YoY[index.bad.u2.W.YoY] <- NA

temp.u2.W.1      <- u2.W.1
index.bad.u2.W.1  <- as.vector((1:length(time)) %*% outer(time,bad.u2.W.1,"=="))
temp.u2.W.1[index.bad.u2.W.1] <- NA

temp.u2.HJ.1      <- u2.H.1
index.bad.u2.H.1  <- as.vector((1:length(time)) %*% outer(time,bad.u2.H.1,"=="))
temp.u2.H.1[index.bad.u2.H.1] <- NA


select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2)) & (!is.na(temp.u2.W.YoY) & (!is.na(temp.u2.W.1)) & (!is.na(temp.u2.H.1)) )
cat("\n\n*** Pooled Petersen Estimate after removing bad m2, u2.* values adjusting for sampling fraction  ***\n\n")
cat("The following strata had m2       set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.W.YoY set to missing: ", 
     if(length(bad.u2.W.YoY)>0){bad.u2.W.YoY} else {" NONE"}, "\n")
cat("The following strata had u2.W.1   set to missing: ", 
     if(length(bad.u2.W.1)>0){bad.u2.W.1} else {" NONE"}, "\n")
cat("The following strata had u2.H.1   set to missing: ", 
     if(length(bad.u2.H.1)>0){bad.u2.H.1} else {" NONE"}, "\n")
cat("\nThe following strata are excluded because n1=0 or NA values in m2, u2.*:", time[!select],"\n\n")

temp.n1        <- n1      [select]
temp.m2        <- m2      [select]
temp.u2.W.YoY  <- u2.W.YoY[select]
temp.u2.W.1    <- u2.W.1  [select]
temp.u2.H.1    <- u2.H.1  [select]

cat("W.YoY  Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.W.YoY),"\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.W.YoY))
cat("W.YoY  Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("W.YoY  Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

cat("W.1    Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.W.1),"\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.W.1))
cat("W.1    Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("W.1    Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

cat("H.1   Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.H.1),"\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.H.1))
cat("H.1   Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("H.1   Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")



# Obtain Petersen estimator for each stratum prior to removing bad m2 values
cat("*** Stratified Petersen Estimator for each stratum PRIOR to removing bad m2 values adjusting for sampling fractions ***\n\n")
cat("W.YoY  raw data\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.W.YoY
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.YoY', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.YoY  Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
                  "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")


cat("W.1   raw data\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.W.1
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.1    Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
                  "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("H.1   raw data\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.H.1
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("H.1    Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
                  "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")




# Obtain Petersen estimator for each stratum after removing bad m2, u2.* values
cat("*** Stratified Petersen Estimator for each stratum AFTER removing bad m2, u2.* values ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.m2[index.bad.m2] <- NA

temp.u2.W.YoY <- u2.W.YoY
temp.u2.W.YoY[index.bad.u2.W.YoY] <- NA

temp.u2.W.1 <- u2.W.1
temp.u2.W.1[index.bad.u2.W.1] <- NA

temp.u2.H.1 <- u2.H.1
temp.u2.H.1[index.bad.u2.H.1] <- NA

cat("\nW.YoY raw data\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.W.YoY)
temp <- cbind(time, temp.n1, temp.m2, temp.u2.W.YoY, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.YoY(adj)', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.YoY Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("\nW.1   raw data\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.W.1)
temp <- cbind(time, temp.n1, temp.m2, temp.u2.W.1, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1(adj)', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.1   Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("\nH.1   raw data\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.H.1)
temp <- cbind(time, temp.n1, temp.m2, temp.u2.H.1, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.H.1(adj)', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("H.1   Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")


# Test if pooling can be done
cat("*** Test if pooled Petersen is allowable. [Check if marked fractions are equal] ***\n\n")
select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2)) 
if(sum(select)<2){
  cat("Test for pooling not done because less than 2 release groups remaining\n")
}
if(sum(select)>=2){
  temp.n1 <- n1[select]
  temp.m2 <- m2[select]
  test <- TestIfPool( temp.n1, temp.m2)
  cat("(Large Sample) Chi-square test statistic ", test$chi$statistic," has p-value", test$chi$p.value,"\n\n")
  temp <- cbind(time[select],test$chi$observed, round(test$chi$expected,1), round(test$chi$residuals^2,1))
  colnames(temp) <- c('time','n1-m2','m2','E[n1-m2]','E[m2]','X2[n1-m2]','X2[m2]')
  print(temp)
  cat("\n Be cautious of using this test in cases of small expected values. \n\n")
}



# Fix up any data problems and prepare for the call.
# Notice that for strata entries that are missing any covariate values, only an intercept is added

# Expand the entries in case of missing time entries
new.n1         <- rep(0, max(time)-min(time)+1)
new.m2         <- rep(0, max(time)-min(time)+1)
new.u2.W.YoY   <- rep(0, max(time)-min(time)+1)
new.u2.W.1     <- rep(0, max(time)-min(time)+1)
new.u2.H.1     <- rep(0, max(time)-min(time)+1)
new.logitP.cov <- matrix(NA, nrow=max(time)-min(time)+1, ncol=ncol(as.matrix(logitP.cov)))
new.time       <- min(time):max(time)


new.n1  [time-min(time)+1]         <- n1
new.m2  [time-min(time)+1]         <- m2
new.m2  [bad.m2-min(time)+1]       <- NA    # wipe out strata where m2 is known to be bad

new.u2.W.YoY[time-min(time)+1]             <- u2.W.YoY
new.u2.W.YoY[bad.u2.W.YoY-min(time)+1]     <- NA    # wipe out strata where u2.W.YoY is known to be bad

new.u2.W.1  [time-min(time)+1]             <- u2.W.1
new.u2.W.1  [bad.u2.W.1-min(time)+1]       <- NA    # wipe out strata where u2.W.1 is known to be bad

new.u2.H.1  [time-min(time)+1]             <- u2.H.1
new.u2.H.1  [bad.u2.H.1-min(time)+1]       <- NA    # wipe out strata where u2.W.1 is known to be bad


new.logitP.cov[time-min(time)+1,]<- as.matrix(logitP.cov)
new.logitP.cov[ is.na(new.logitP.cov[,1]), 1] <- 1  # insert a 1 into first columns where not specified
new.logitP.cov[ is.na(new.logitP.cov)] <- 0         # other covariates are forced to zero not in column 1


# Check for and fix problems with the data
# If n1=m2=0, then set n1 to 1, and set m2<-NA
new.m2[new.n1==0] <- NA
new.n1[new.n1==0] <- 1

# Print out the revised data
hatch.indicator <- rep('   ', max(time)-min(time)+1)
hatch.indicator[hatch.after-min(time)+1]<- '***'

cat("\n\n*** Revised data *** \n")
temp<- data.frame(time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, 
       new.logitP.cov=new.logitP.cov, 
       hatch.indicator=hatch.indicator)
print(temp) 
cat("\n\n")

# Print out information on the prior distributions used
cat("\n\n*** Information on priors *** \n")
cat("   Parameters for prior on tauU (variance in spline coefficients: ", tauU.alpha, tauU.beta, 
    " which corresponds to a mean/std dev of 1/var of:",
    round(tauU.alpha/tauU.beta,2),round(sqrt(tauU.alpha/tauU.beta^2),2),"\n")
cat("   Parameters for prior on taueU (variance of log(U) about spline: ",taueU.alpha, taueU.beta, 
    " which corresponds to a mean/std dev of 1/var of:",
    round(taueU.alpha/taueU.beta,2),round(sqrt(taueU.alpha/taueU.beta^2),2),"\n")
cat("   Parameters for prior on beta.logitP[1] (intercept) (mean, sd): \n", cbind(round(prior.beta.logitP.mean,3), round(prior.beta.logitP.sd,5)),"\n")
cat("   Parameters for prior on tauP (residual variance of logit(P) after adjusting for covariates: ",tauP.alpha, tauP.beta, 
    " which corresponds to a mean/std dev of 1/var of:",
    round(tauP.alpha/tauP.beta,2),round(sqrt(tauP.alpha/tauP.beta^2),2),"\n")

cat("\n\nInitial seed for this run is: ",InitialSeed, "\n")

sink()

if (debug2) {
   cat("\nprior to formal call to TimeStratPetersenDiagErrorWH\n")
   browser()
}


if (debug) 
   {results <- TimeStratPetersenDiagErrorWHSteel(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, 
            hatch.after=hatch.after-min(time)+1, 
            logitP.cov=new.logitP.cov,
            n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,   # set to low value for debugging only
            prior.beta.logitP.mean=prior.beta.logitP.mean, 
            prior.beta.logitP.sd  =prior.beta.logitP.sd,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, debug2=debug2, InitialSeed=InitialSeed,
            save.output.to.files=save.output.to.files)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenDiagErrorWHSteel(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, 
            hatch.after=hatch.after-min(time)+1, 
            logitP.cov=new.logitP.cov,
            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims, 
            prior.beta.logitP.mean=prior.beta.logitP.mean, 
            prior.beta.logitP.sd  =prior.beta.logitP.sd,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, debug2=debug2, InitialSeed=InitialSeed,
            save.output.to.files=save.output.to.files)
  }

# Now to create the various summary tables of the results

#browser()
# A plot of the observered log(U) on the log scale, and the final mean log(U)
# in the diagonal case, all of n1, m2, u2 have the same length
   Nstrata <- length(n1)

  plot.df   <- data.frame(time =new.time)
 
  # adjust the u2 for the clipping fractions
  plot.df$n1       <- new.n1
  plot.df$m2       <- new.m2
  plot.df$u2.H.1   <- new.u2.H.1 
  plot.df$u2.W.1   <- new.u2.W.1
  plot.df$u2.W.YoY <- new.u2.W.YoY

  plot.df$u2.W.YoY[is.na(plot.df$u2.W.YoY)] <- 1  # in case of missing values
  plot.df$u2.H.1  [is.na(plot.df$u2.H.1)  ] <- 1  # in case of missing values
  plot.df$u2.W.1  [is.na(plot.df$u2.W.1)  ] <- 1  # in case of missing values

  get.est <- function(est.name, plot.df, hatch.after, results){
      # get the inital estimates, and extract from the results data structure and put into a data frame
      est.df <- data.frame(group=est.name, time=plot.df$time)
      avgP <- sum(plot.df$m2,na.rm=TRUE)/sum(plot.df$n1, na.rM=TRUE)
      #browser()
      # initial guess
      est.df$logUguess <- log(1+pmax( (plot.df[, paste("u2.",est.name,sep="")]+1)*(plot.df$n1+2)/(plot.df$m2+1), 
                                       plot.df[, paste("u2.",est.name,sep="")]/avgP, na.rm=TRUE))
      # extract estimates from results
      results.row.names <- rownames(results$summary)
      est.row.index    <- grep(paste("etaU.",est.name, sep=""), results.row.names)
      etaU <- results$summary[est.row.index,]
      est.df$logU    =etaU[,"mean"]
      est.df$logUlcl =etaU[,"2.5%"]
      est.df$logUucl =etaU[,"97.5%"]
      # extract the spline
      logUne.row.index <- grep(paste("logUne.",est.name,sep=""), results.row.names)
      est.df$spline    <- results$summary[logUne.row.index,"mean"]

      if(est.name=="H.1"){
         est.df$logUguess[1:(hatch.after-min(plot.df$time)+1)]<- NA
         est.df$logU     [1:(hatch.after-min(plot.df$time)+1)]<- NA
         est.df$logUlcl  [1:(hatch.after-min(plot.df$time)+1)]<- NA
         est.df$logUucl  [1:(hatch.after-min(plot.df$time)+1)]<- NA
         est.df$spline   [1:(hatch.after-min(plot.df$time)+1)]<- NA
      }
      est.df
  }
  plot.data <-rbind( get.est("H.1"  ,plot.df, hatch.after, results),
                     get.est("W.YoY",plot.df, hatch.after, results),
                     get.est("W.1"  ,plot.df, hatch.after, results))
  
  # add limits to the plot to avoid non-monotone secondary axis problems with extreme values
  plot.data$logUguess <- pmax(-10 , pmin(20, plot.data$logUguess))
  plot.data$logU      <- pmax(-10 , pmin(20, plot.data$logU ))
  plot.data$logUlcl   <- pmax(-10 , pmin(20, plot.data$logUlcl  ))
  plot.data$logUucl   <- pmax(-10 , pmin(20, plot.data$logUucl  ))
  plot.data$spline    <- pmax(-10 , pmin(20, plot.data$spline))

  fit.plot <- ggplot(data=plot.data, aes_(x=~time, color=~group))+
     ggtitle(title, subtitle="Fitted spline curve with 95% credible intervals")+
     geom_point(aes_(y=~logUguess), shape=16, position=position_dodge(width=.2))+  # guesses for population
     geom_point(aes_(y=~logU), shape=19, position=position_dodge(width=.2))+
     geom_line (aes_(y=~logU), position=position_dodge(width=.2))+
     geom_errorbar(aes_(ymin=~logUlcl, ymax=~logUucl), width=.1, position=position_dodge(width=.2))+
     geom_line(aes_(y=~spline),linetype="dashed", position=position_dodge(width=.2)) + 
     xlab("Time Index\nFitted/Smoothed/Raw values plotted for W(black) and H(blue)")+
     ylab("log(U[i]) + 95% credible interval")+
     theme(legend.justification = c(0, 0), legend.position = c(0, 0))+
     scale_color_discrete(name="Group")+
     scale_x_continuous(breaks=seq(min(plot.data$time, na.rm=TRUE),max(plot.data$time, na.rm=TRUE),2))+
     scale_y_continuous(sec.axis = sec_axis(~ exp(.), name="U + 95% credible interval",
                      breaks=c(1,10,20,50,
                                 100,200,500,
                                 1000,2000,5000,
                                 10000,20000, 50000,
                                 100000,200000, 500000,
                                 1000000,2000000,5000000,10000000),
                      labels = scales::comma))



if(save.output.to.files)ggsave(plot=fit.plot, filename=paste(prefix,"-fit.pdf",sep=""), height=6, width=10, units="in")
results$plots$fit.plot <- fit.plot


# plot the estimated logits over time
logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=new.m2, 
               u2=new.u2.W.YoY+new.u2.W.1+new.u2.H.1, logitP.cov=new.logitP.cov, results=results,
               trunc.logitP=trunc.logitP)
if(save.output.to.files)ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
results$plots$logitP.plot <- logitP.plot


# Look at autocorrelation function for Utot.W.YoY and Utot.W.1, U.tot.H.1
mcmc.sample1<- data.frame(parm="Utot.W.YoY", sample=results$sims.matrix[,"Utot.W.YoY"], stringsAsFactors=FALSE)
mcmc.sample2<- data.frame(parm="Utot.W.1",   sample=results$sims.matrix[,"Utot.W.1"], stringsAsFactors=FALSE)
mcmc.sample3<- data.frame(parm="Utot.H.1",   sample=results$sims.matrix[,"Utot.H.1"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2, mcmc.sample3)
acf.Utot.plot <- plot_acf(mcmc.sample)
if(save.output.to.files)ggsave(plot=acf.Utot.plot, filename=paste(prefix,"-Utot-acf.pdf",sep=""), height=4, width=6, units="in")
results$plots$acf.Utot.plot <- acf.Utot.plot


# Look at posterior plot for Utot.W.YoY and Utot.W.1, U.tot.H.1
mcmc.sample1<- data.frame(parm="Utot.W.YoY", sample=results$sims.matrix[,"Utot.W.YoY"], stringsAsFactors=FALSE)
mcmc.sample2<- data.frame(parm="Utot.W.1",   sample=results$sims.matrix[,"Utot.W.1"], stringsAsFactors=FALSE)
mcmc.sample3<- data.frame(parm="Utot.H.1",   sample=results$sims.matrix[,"Utot.H.1"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2, mcmc.sample3)
post.Utot.plot <- plot_posterior(mcmc.sample)
if(save.output.to.files)ggsave(plot=post.Utot.plot, filename=paste(prefix,"-Utot-posterior.pdf",sep=""), height=4, width=6, units="in")
results$plots$post.Utot.plot <- post.Utot.plot


#save the Bayesian predictive distribution (Bayesian p-value plots)
#browser()
discrep <-PredictivePosterior.TSPDE.WHSteel (time, new.n1, new.m2, new.u2.W.YoY, new.u2.W.1, new.u2.H.1, 
          expit(results$sims.list$logitP), 
          round(results$sims.list$U.W.YoY), 
          round(results$sims.list$U.W.1), 
          round(pmax(results$sims.list$U.H.1,0)), hatch.after) #don't forget that hatchery fish is 0 until hatch.after
gof <- PredictivePosteriorPlot.TSPDE.WHSteel(discrep)
if(save.output.to.files)ggsave(gof[[1]],filename=paste(prefix,"-GOF.pdf",sep=""),   height=8, width=8, units="in", dpi=300 )
results$plots$gof.plot <- gof


# create traceplots of logU, U, and logitP (along with R value) to look for non-convergence
# the plot_trace will return a list of plots (one for each page as needed)
varnames <- names(results$sims.array[1,1,])  # extract the names of the variables 

# Trace plots of logitP
trace.plot <- plot_trace(title=title, results=results, parms_to_plot=varnames[grep("^logitP", varnames)])
if(save.output.to.files){
   pdf(file=paste(prefix,"-trace-logitP.pdf",sep=""))
   plyr::l_ply(trace.plot, function(x){plot(x)})
   dev.off()
}
results$plots$trace.logitP.plot <- trace.plot

# now for the traceplots of logU (etaU), Utot, and Ntot
trace.plot <- plot_trace(title=title, results=results, parms_to_plot=varnames[c(grep("Utot",varnames), grep("Ntot",varnames), grep("^etaU", varnames))])
if(save.output.to.files){
   pdf(file=paste(prefix,"-trace-logU.pdf",sep=""))
   plyr::l_ply(trace.plot, function(x){plot(x)})
   dev.off()
}
results$plots$trace.logU.plot <- trace.plot



sink(report, append=TRUE)
# What was the initial seed
cat("\n\n*** Initial State for this run ***: ", results$Seed.initial,"\n")
cat("*** See help(modelSetRN) for details. ***\n")

# Global summary of results
cat("\n\n*** Summary of MCMC results *** \n\n")
  save.max.print <- getOption("max.print")
  options(max.print=.Machine$integer.max)
  
  print(results, digits.summary=3)#, max=.Machine$integer.max)
  
  options(max.print=save.max.print)
  
  
# Give an alternate computation of DIC based on the variance of the deviance
# Refer to http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/DIC-slides.pdf for derivation and why
# this alternate method may be superior to that automatically computed by WinBugs/OpenBugs

cat("\n\n*** Alternate DIC computation based on p_D = var(deviance)/2 \n")
results.row.names <- rownames(results$summary)
deviance.row.index<- grep("deviance", results.row.names)
deviance          <- results$summary[deviance.row.index,]
p.D <- deviance["sd"]^2/2
dic <- deviance["mean"]+p.D
cat("    D-bar: ", deviance["mean"],";  var(dev): ", deviance["sd"]^2,
    "; p.D: ", p.D, "; DIC: ", dic)

# Summary of population sizes. Add pretty printing
cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
cat("W.YoY \n")
temp <- results$summary[ grep("Utot.W.YoY", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nW.1   \n")
temp<- results$summary[ grep("Utot.W.1",   rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nH.1\n")
temp<- results$summary[ grep("Utot.H.1",   rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)


#browser()
time.H <- time>hatch.after
cat("\n\n\n\n*** Summary of Quantiles of Run Timing.Wild.YoY *** \n")
cat(    "    This is based on the sample weeks provided and the U.W.YoY[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W.YoY, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n*** Summary of Quantiles of Run Timing.Wild.1 *** \n")
cat(    "    This is based on the sample weeks provided and the U.W.1[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W.1, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n*** Summary of Quantiles of Run Timing.Hatchery.1 *** \n")
cat(    "    This is based on the sample weeks provided and the U.H.1[i] values \n") 
q <- RunTime(time=time[time.H], U=results$sims.list$U.H.1[,time.H], prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

# Add the runtiming to the output object
results$runTime <- temp


cat("\n\n")
cat(paste("*** end of fit *** ", date()))

sink()

# save the report to a files?
if(save.output.to.files)writeLines(stdout, results.filename)
results$report <- stdout


# add some of the raw data to the bugs object for simplicity in referencing it later
results$data <- list( time=time, n1=n1, m2=m2, u2.W.YoY=u2.W.YoY, u2.W.1=u2.W.1, u2.H.1=u2.H.1,
                      hatch.after=hatch.after, 
                      bad.m2=bad.m2, bad.u2.W.YoY=bad.u2.W.YoY, bad.u2.W.1=bad.u2.W.1, bad.u2.H.1=bad.u2.H.1, 
                      logitP.cov=logitP.cov,
                      version=version, date_run=date(),title=title)

return(results)
} # end of function
