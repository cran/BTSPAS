# 2024-05-09 CJS TestIfPool only if at least 2 valid release groups
# 2021-10-23 CJS Added trunc.logitP to deal with extreme values of logitP during plotting
# 2020-12-15 CJS Removed sampfrac from code
# 2020-12-14 CJS All bad.n1 or bad.m2 are set to 0. No NA allows for n1 or m2
#                Fixed problem with some u2=NA in the GOF computations and plots
# 2020-11-07 CJS Allowed user to specify prior for beta coefficient for logitP
# 2108-12-19 CJS deprecate of use of sampling fraction
# 2018-12-15 CJS converted the muTT and sdTT plots to ggplot objects
# 2018-12-15 CJS tested fixing logitP to certain values especially at the end of the sampling chain
# 2018-12-06 CJS converted report to textConnection() object
# 2018-12-03 CJS converted fit plot to ggplot2
# 2018-12-02 CJS converted trace plots to ggplot2
# 2018-12-01 CJS Converted acf, posterior plots to ggplot2
# 2018-11-28 CJS Fixed problem of large printing being cutoff in results
# 2018-11-25 CJS Removed all references to OpenBugs
# 2014-09-01 CJS converstion to jags
# 2012-08-30 CJS fixed NAs problem in any() and all() in error checking
# 2011-06-13 CJS added p-values to results
# 2010-02-21 CJS changed u2 to new.u2 in expanded.m2 section
# 2010-11-25 CJS pretty printing of final population estimates
# 2010-09-06 CJS forced input vectors to be vectors
# 2010-08-06 CJS produced traceplot
# 2010-08-04 CJS added version/date to final result
# 2010-03-12 CJS added n.chains etc to calling arguments
# 2010-03-03 CJS allowed the user for fix some logitPs to account for time when no sampling (or other events
#                with known probability of capture take place
# 2009-12-08 CJS added some basic error checking on arguments
# 2009-12-07 CJS added bad.n1, bad.u2 arguments and fixups
# 2009-12-01 CJS added open/winbugs directory to argument list




#' Wrapper (*_fit) to fit the Time Stratified Petersen Estimator
#' with NON Diagonal Entries function.
#' 
#' Takes the number of marked fish released, the number of recaptures, and the
#' number of unmarked fish and uses Bayesian methods to fit a fit a spline
#' through the population numbers and a hierarchical model for the trap
#' efficiencies over time.  The output is written to files and an MCMC object
#' is also created with samples from the posterior.
#' 
#' Normally the user makes a call to the *_fit function which then calls the
#' fitting function.
#' 
#' Use the \code{\link{TimeStratPetersenDiagError_fit}} function for cases
#' where recaptures take place ONLY in the stratum of release, i.e. the
#' diagonal case.
#' 
#' The non-diagonal case fits a log-normal distribution for the travel time.
#' The *NP functions fit a non-parametric distribution for the travel times.
#' The *MarkAvail functions extend the *NP functions to allow for reductions in
#' mark availability because of fall back, immediate tagging mortality, etc.
#' 
#' 
#' @template title
#' @template prefix
#' @template time
#' @template n1 
#' @param m2 A numeric matrix of the number of fish released in stratum [i] and
#' recovered in [j-1] strata later.  For example m2[3,5] is the number of
#' marked fish released in stratum 3 and recovered 4 strata later in stratum 7.
#' The first column is the number of marked fish recovered in the stratum of
#' release, i.e. 0 strata later.  Use the
#' \code{\link{TimeStratPetersenDiagError_fit}} function for cases where
#' recaptures take place ONLY in the stratum of release, i.e. the diagonal
#' case.
#' @template u2.ND
#' @template sampfrac
#' @template jump.after
#' @template bad.n1
#' @template bad.m2 
#' @template bad.u2 
#' @template logitP.cov
#' @template param.logitP.fixed
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

#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' @template author 
#' @template references
#' @keywords ~models ~smooth
#' @examples
#'  
#' ##---- See the vignettes for examples of how to use this package
#' 
#' @export TimeStratPetersenNonDiagError_fit
#' @importFrom stats runif var sd

TimeStratPetersenNonDiagError_fit <-
  function( title="TSPNDE",
           prefix="TSPNDE-", time, n1, m2, u2,
           sampfrac=rep(1,length(u2)), jump.after=NULL,
           bad.n1=c(), bad.m2=c(), bad.u2=c(),
           logitP.cov=as.matrix(rep(1,length(u2))),
           logitP.fixed=NULL, logitP.fixed.values=NULL, 
           n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
           tauU.alpha=1,tauU.beta=.05,
           taueU.alpha=1, taueU.beta=.05,
           prior.beta.logitP.mean = c(logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),rep(0,  ncol(as.matrix(logitP.cov))-1)),
           prior.beta.logitP.sd   = c(2,                                           rep(10, ncol(as.matrix(logitP.cov))-1)), 
           tauP.alpha=.001,
           tauP.beta=.001,
           run.prob=seq(0,1,.1), # what percentiles of run timing are wanted
           debug=FALSE, debug2=FALSE,
           InitialSeed=ceiling(stats::runif(1,min=0,1000000)),
           save.output.to.files=TRUE,
           trunc.logitP=15) {
# Fit a Time Stratified Petersen model with NON-diagonal entries and with smoothing on U allowing for random error
# This is the classical stratified Petersen model where the recoveries can take place for this and multiple
# strata later
#
    version <- as.character(packageDate("BTSPAS"))
    options(width=200)

# Input parameters are
#    title  - title for the analysis (character string)
#    prefix - prefix used for files created with the analysis results
#             this should be in standard Window's format, eg. JC-2002-ST-TSPNDE
#             to which is appended various suffixes for plots etc (character string)
#    time   - vector of stratum numbers. For example, 9:38 would indicate that the
#             Trinity River system sampled weeks 9 to 38.
#             These values are never used in the analysis and only serve as labels for the weeks and for plotting purposes.
#             They should be contiguous equally spaced values and be the same length as u2.
#    n1, m2, u2 - the input data consisting of fish marked and released, recapture, and unmarked captured
#             Note that m2 is a MATRIX. The first column are those fish recaptured in the stratum of release
#             and the subsequent columns are those recoveries in later strata.
#             This is expanded to the full matrix [i,j] for released released in stratum i and recovered in stratum j
#             The vector u2 should be long enough to account for any fish that are recaptured later on
#             from releases late in the season. The bottom right diagonal of m2 may be all zeros - that is ok
#             Notice that length(u2) can be longer than length(n1)+nrow(m2).
#    sampfrac - Deprecated. DO NOT USE ANYMORE. 
#    jump.after - in some cases, a single spline is still not flexible enough to cope with rapid
#                 changes in the run curve. For example, in the Trinity River project, a larger
#                 hatchery release occurs around stratum 14. This is a vector indicating the
#                 strata AFTER which the spline curve is allowed to jump.
#                 null or vector of arbitrary length.
#    bad.n1  - vector of stratum numbers where the value of n1 is suspect.
#    bad.m2  - vector of stratum numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to  JAGS
#    bad.u2  - vector of stratum numbers where the value of u2 is suspect.
#    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
#               for the first element of the covariance matrix and 0 for the rest of the covariates.
#               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
#               if you have any covariates. The default, if not specified, is a constant (the mean logit)
#    logitP.fixed - vector of time values where the logitP will be specified in advance. Typically this occurs
#                   when no sampling takes place and logitP <- logit(0) = -10
#    logitP.fixed.values - vector of fixed values (on the logit scale). Use -10 for p[i] <- 0, and 10 for p[i] <- 1
#    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
#    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline
#    prior.beta.logitP.mean, prior.beta.logitP.sd   - parameters for the prior on mean logit(P)'s [The intercept term]
#                              The other covariates are assigned priors of a mean of 0 and a sd of 30
#    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
#    run.prob  - percentiles of run timing wanted
#    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data
#             and JAGS will run and stop waiting for your to exit and complete

# force input vectors to be vectors as needed. Note that m2 is NOT a vector!
time     <- as.vector(time)
n1       <- as.vector(n1)
u2       <- as.vector(u2)
sampfrac <- as.vector(sampfrac)

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time are consistent with each other.
#  In the non-diagonal case, they don't have to match
if(length(n1)!=nrow(m2)){
   cat("***** ERROR ***** Length of n1 and number of rows of m2 must be equal. They are:",
        length(n1)," ",nrow(u2),"\n")
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
if(stats::var(c(length(u2),length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of u2, sampfrac, time must all be equal. They are:",
        length(u2)," ",length(sampfrac)," ",length(time),"\n")
   return()}
if(length(logitP.cov) %% length(u2) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of u2. They are:",
        length(u2)," ",length(logitP.cov)," ",dim(logitP.cov),"\n")
   return()}
#  2. Check that rowsum of m2<= n1
if(any(apply(m2,1,sum, na.rm=TRUE)>n1)){
   cat("***** ERROR ***** m2[i,+] must be <= n1[i]. The arguments are \n n1:",
       paste(n1,collapse=","),"\n m2:",
       paste(m2,collapse=","),"\n")
   return()}
#  3. Elements of bad.m2 and jump.after must belong to time
if(!all(bad.n1 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.n1 must be elements of strata identifiers. You entered \n bad.n1:",
       paste(bad.n1,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,  collapse=","), "\n")
   return()}
if(!all(bad.m2 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",
       paste(bad.m2,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,  collapse=","), "\n")
   return()}
if(!all(bad.u2 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2 must be elements of strata identifiers. You entered \n bad.u2:",
       paste(bad.u2,collapse=","),"\n Strata identifiers are \n time:",
       paste(time  ,collapse=","), "\n")
   return()}
if(!all(jump.after %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** jump.after must be elements of strata identifiers. You entered \n jump.after:",
       paste(jump.after,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,      collapse=","), "\n")
   return()}

#  4. check that index of logitP.fixed belong to time
if(!all(logitP.fixed %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** logitP.fixed must be elements of strata identifiers. You entered \n logitP.fixed:",
       paste(logitP.fixed,collapse=","),"\n Strata identifiers are \n time:",
       paste(time        ,collapse=","), "\n")
   return()}
if(length(logitP.fixed)!=length(logitP.fixed.values)){
   cat("***** ERROR ***** Lengths of logitP.fixed and logitP.fixed.values must all be equal. They are:",
        length(logitP.fixed)," ",length(logitP.fixed.values),"\n")
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

cat(paste("Time Stratified Petersen with Non-Diagonal recaptures, error in smoothed U, and log-normal distribution for travel time - ", date()))
cat("\nVersion: ", version)

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** (padded to match length of u2) \n")
jump.indicator <- rep('   ', length(u2))
jump.indicator[time %in% jump.after]<- '***'
ex.n1 <- c(n1, rep(NA, length(u2)-length(n1)))
ex.m2 <- rbind(m2,matrix(NA, nrow=length(u2)-length(n1), ncol=ncol(m2)))
temp<- data.frame(time=time, n1=ex.n1, m2=ex.m2, u2=u2, logitP.cov=logitP.cov, jump=jump.indicator)
print(temp)
cat("\n\n")
cat("Jump point are after strata: ", jump.after)
if(length(jump.after)==0) cat("none - A single spline is fit")
cat("\nFixed logitP indices are: ", logitP.fixed)
if(length(logitP.fixed)==0) cat("none - NO fixed values")
cat("\nFixed logitP values  are: ", logitP.fixed.values)
if(length(logitP.fixed)==0) cat("none - NO fixed values")


# Obtain the Pooled Petersen estimator prior to fixup of bad.n1, bad.m2, and bad.u2 values
cat("\n\n*** Pooled Petersen Estimate prior to fixing bad n1, m2, or u2 values ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2

cat("Total n1=", sum(temp.n1,na.rm=TRUE),";  m2=",sum(temp.m2,na.rm=TRUE),";  u2=",sum(temp.u2,na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(temp.n1,na.rm=TRUE), sum(temp.m2,na.rm=TRUE), sum(temp.u2,na.rm=TRUE))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")
.

# Obtain the Pooled Petersen estimator after removal of entries with bad.n1, m2, or u2 values
select.rel <- !(time[1:length(n1)] %in% bad.n1 | time[1:length(n1)] %in% bad.m2 )
select.rec <- ! time %in% bad.u2
cat("\n\n*** Pooled Petersen Estimate after removing release and recovery strata flagged as bad ***\n\n")
cat("The following release strata were excluded:",
     if(length(time[!select.rel])>0){time[!select.rel]} else {" NONE"}, "\n")
cat("The following recovery strata were excluded:",
     if(length(time[!select.rec])>0){time[!select.rec]} else {" NONE"}, "\n")

temp.n1 <- n1[select.rel]
temp.m2 <- m2[select.rel]
temp.u2 <- u2[select.rec]

cat("Total n1=", sum(temp.n1,na.rm=TRUE),";  m2=",sum(temp.m2,na.rm=TRUE),";  u2=",sum(temp.u2, na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(temp.n1,na.rm=TRUE), sum(temp.m2,na.rm=TRUE), sum(temp.u2,na.rm=TRUE))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")




# Test if pooling can be done
# We only do the release strata that are not flagged as bad and have no missing values

cat("*** Test if pooled Petersen is allowable. [Check if equal recovery from each stratum not flagged and without missing recoveries] ***\n\n")
#browser()
select <- select.rel  & (!is.na(apply(m2,1,sum)))
if(sum(select)<2){
  cat("Test for pooling not done because less than 2 release groups remaining\n")
}
if(sum(select)>=2){
  temp.n1 <- n1[select]
  temp.m2 <- m2[select,]
  test <- TestIfPool( temp.n1, apply(temp.m2,1,sum, na.rm=TRUE))
  cat("(Large Sample) Chi-square test statistic ", test$chi$statistic," has p-value", test$chi$p.value,"\n\n")
  temp <- cbind(time[select],test$chi$observed, round(test$chi$expected,1), round(test$chi$residuals^2,1))
  colnames(temp) <- c('time','n1-m2','m2','E[n1-m2]','E[m2]','X2[n1-m2]','X2[m2]')
  print(temp)
  cat("\n Be cautious of using this test in cases of small expected values. \n\n")
}


# Adjust the data for the explicitly bad values or other problems
new.time <- time
new.n1   <- n1
new.m2   <- m2
new.u2   <- u2
new.logitP.cov <- logitP.cov


# Set the bad n1 values to 0 for the number of fish released and corresponding values of m2 also to 0 recovered subsequently.
# Set any bad m2 values to 0 for the number of releases and subsequent recoveries as well.
# But we don't set bqd(u2) values to 0 as this would imply no catch. We set these to missing
new.n1[time[1:length(n1)] %in% c(bad.n1,bad.m2) ]  <- 0
new.m2[time[1:length(n1)] %in% c(bad.m2,bad.n1),]  <- 0
new.u2[time %in% bad.u2]                <- NA

# Print out the revised data
cat("\n\n*** Revised data after setting flagged strata to 0 (releases) or NA (recoveries) *** \n")
jump.indicator <- rep('   ', length(u2))
jump.indicator[time %in% jump.after]<- '***'
ex.n1 <- c(new.n1, rep(NA, length(new.u2)-length(new.n1)))
ex.m2 <- rbind(new.m2,matrix(NA, nrow=length(new.u2)-length(new.n1), ncol=ncol(new.m2)))
temp<- data.frame(time=new.time, n1=ex.n1, m2=ex.m2, u2=new.u2, logitP.cov=new.logitP.cov,
       jump.after=jump.indicator)
print(temp)
cat("\n\n")

# Need to expand the m2 matrix to the full Nstrata.rel x Nstrata.cap+1 matrix with
# the [i,j] entry for the number of fish released in stratum i and recaptured in stratum j
# The last column of the expanded m2 matrix is the number of fish released but never recaptured later.

expanded.m2 <- matrix(0, nrow=length(new.n1), ncol=length(new.u2)+1)
for(i in 1:length(new.n1)){
   expanded.m2[i,1:length(new.u2)] <- c(rep(0,i-1),new.m2[i,],rep(0,length(new.u2)))[1:length(new.u2)]
   expanded.m2[i,length(new.u2)+1] <- new.n1[i] - sum(new.m2[i,])
}

cat("*** Expanded m2 array ***\n\n")
save.max.print <- getOption("max.print")
options(max.print=.Machine$integer.max)
print(expanded.m2)
options(max.print=save.max.print)


# assign the logitP fixed values etc.
new.logitP.fixed <- rep(NA, length(new.u2))
new.logitP.fixed[match(logitP.fixed, time)] <- logitP.fixed.values


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
   cat("\nprior to formal call to TimeStratPetersenNonDiagError\n")
   browser()
}

if (debug)
  {results <- TimeStratPetersenNonDiagError(
        title=title, prefix=prefix,
        time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2,
        jump.after=(1:length(u2))[time %in% jump.after],
        logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
        n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,  # set to small values for debugging only
        prior.beta.logitP.mean=prior.beta.logitP.mean, 
        prior.beta.logitP.sd  =prior.beta.logitP.sd,
        tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
        debug=debug, debug2=debug2, InitialSeed=InitialSeed, save.output.to.files=save.output.to.files)
 }
else #notice R syntax requires { before the else
  {results <- TimeStratPetersenNonDiagError(
        title=title, prefix=prefix,
        time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2,
        jump.after=(1:length(u2))[time %in% jump.after],
        logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
        n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims,
        prior.beta.logitP.mean=prior.beta.logitP.mean, 
        prior.beta.logitP.sd  =prior.beta.logitP.sd,
        tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
        debug=debug, debug2=debug2, InitialSeed=InitialSeed, save.output.to.files=save.output.to.files)
 }

# Now to create the various summary tables of the results

Nstrata.rel <- length(n1)
Nstrata.cap <- ncol(expanded.m2) -1 # don't forget that last column of m2 is number of fish never seen

# A plot of the observed log(U) on the log scale, and the final mean log(U)
plot.df   <- data.frame(time =new.time)
plot.df$logUi <-log( c((new.u2[1:Nstrata.rel]+1)*(new.n1+2)/(apply(expanded.m2[,1:Nstrata.cap],1,sum)+1), rep(NA, length(u2)-Nstrata.rel)))

# extract the fitted U values
results.row.names <- rownames(results$summary)
etaU.row.index    <- grep("etaU", results.row.names)
etaU<- results$summary[etaU.row.index,]
plot.df$logU    =etaU[,"mean"]
plot.df$logUlcl =etaU[,"2.5%"]
plot.df$logUucl =etaU[,"97.5%"]

# extract the spline values
logUne.row.index <- grep("logUne", results.row.names)
logUne<- results$summary[logUne.row.index,"mean"]
plot.df$spline <- results$summary[logUne.row.index,"mean"]

# add limits to the plot to avoid non-monotone secondary axis problems with extreme values
   plot.df$logUi     <- pmax(-10 , pmin(20, plot.df$logUi))
   plot.df$logU      <- pmax(-10 , pmin(20, plot.df$logU ))
   plot.df$logUlcl   <- pmax(-10 , pmin(20, plot.df$logUlcl  ))
   plot.df$logUucl   <- pmax(-10 , pmin(20, plot.df$logUucl  ))
   plot.df$spline    <- pmax(-10 , pmin(20, plot.df$spline))

fit.plot <- ggplot(data=plot.df, aes_(x=~time))+
     ggtitle(title, subtitle="Fitted spline curve with 95% credible intervals for estimated log(U[i])")+
     geom_point(aes_(y=~logUi), color="red", shape=1)+  # open circle
     xlab("Time Index\nOpen/closed circles - initial and final estimates")+
     ylab("log(U[i]) + 95% credible interval")+
     geom_point(aes_(y=~logU), color="black", shape=19)+
     geom_line (aes_(y=~logU), color="black")+
     geom_errorbar(aes_(ymin=~logUlcl, ymax=~logUucl), width=.1)+
     geom_line(aes_(y=~spline),linetype="dashed")+
     scale_x_continuous(breaks=seq(min(plot.df$time, na.rm=TRUE),max(plot.df$time, na.rm=TRUE),2))+
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


# plot of the logit(p) over time
logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2, 
                           logitP.cov=new.logitP.cov, results=results,
                           trunc.logitP=trunc.logitP)
if(save.output.to.files)ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
results$plots$logitP.plot <- logitP.plot


# Look at autocorrelation function for Ntot
mcmc.sample <- data.frame(parm="Utot", sample=results$sims.matrix[,"Utot"], stringsAsFactors=FALSE)
acf.Utot.plot <- plot_acf(mcmc.sample)
if(save.output.to.files)ggsave(plot=acf.Utot.plot, filename=paste(prefix,"-Utot-acf.pdf",sep=""), height=4, width=6, units="in")
results$plots$acf.Utot.plot <- acf.Utot.plot


# plot the posterior plots
mcmc.sample1 <- data.frame(parm="Utot", sample=results$sims.matrix[,"Utot"], stringsAsFactors=FALSE)
mcmc.sample2 <- data.frame(parm="Ntot", sample=results$sims.matrix[,"Ntot"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2)
post.UNtot.plot <- plot_posterior(mcmc.sample)
post.UNtot.plot
if(save.output.to.files)ggsave(plot=post.UNtot.plot, filename=paste(prefix,"-UNtot-posterior.pdf",sep=""),
                               height=ifelse(length(unique(mcmc.sample$parm))<=2,4,6), width=6, units="in")
results$plots$post.UNtot.plot <- post.UNtot.plot


# plot the mean and sd log(travel times) (the muLogTT and the sdLogTT) vs release stratum number
#browser()
results.row.names <- rownames(results$summary)
muLogTT.row.index    <- grep("muLogTT", results.row.names)
muLogTT<- data.frame(results$summary[muLogTT.row.index,])
muLogTT$stratum <- new.time[1:Nstrata.rel]
muLogTT$source  <- "Mean log(travel time)"

results.row.names <- rownames(results$summary)
sdLogTT.row.index    <- grep("sdLogTT", results.row.names)
sdLogTT<- data.frame(results$summary[sdLogTT.row.index,])
sdLogTT$stratum <- new.time[1:Nstrata.rel]
sdLogTT$source  <- "SD log(travel time)"

plotdata <- rbind( muLogTT, sdLogTT)

musdLogTT.plot <- ggplot(data=plotdata, aes_(x=~stratum, y=~mean))+
   ggtitle(title)+
   geom_point()+
   geom_line()+
   geom_errorbar(aes_(ymin=~X2.5., ymax=~X97.5.), width=.1)+
   xlab("Stratum")+ylab("mean/sd log(travel time)")+
   facet_wrap(~source, ncol=1, scales="free_y")
musdLogTT.plot

if(save.output.to.files)ggsave(plot=musdLogTT.plot, filename=paste(prefix,"-musdLogTT.pdf",sep=""),
                               height=6, width=6, units="in")
results$plots$musdLogTTt <- musdLogTT.plot


## save the Bayesian predictive distribution (Bayesian p-value plots)

discrep <-PredictivePosterior.TSPNDE (new.n1, expanded.m2, new.u2,
                                      new.logitP.fixed,
                                      expit(results$sims.list$logitP),
                                      round(results$sims.list$U),
                                      results$sims.list$muLogTT,
                                      results$sims.list$sdLogTT)
#browser()
gof <- PredictivePosteriorPlot.TSPNDE (discrep)
if(save.output.to.files)ggsave(gof[[1]],filename=paste(prefix,"-GOF.pdf",sep=""),  height=12, width=8, units="in", dpi=300 )
results$plots$gof.plot <- gof

# create traceplots of logU, U, and logitP (along with R value) to look for non-convergence
# the plot_trace will return a list of plots (one for each page as needed)
varnames <- names(results$sims.array[1,1,])  # extract the names of the variables 

#browser()
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

#browser()
# Global summary of results
cat("\n\n*** Summary of MCMC results *** \n\n")
  save.max.print <- getOption("max.print")
  options(max.print=.Machine$integer.max)
  
  print(results, digits.summary=3)#, max=.Machine$integer.max)
  
  options(max.print=save.max.print)
  
  
cat("\n\n*** Alternate DIC computation based on p_D = var(deviance)/2 \n")
results.row.names <- rownames(results$summary)
deviance.row.index<- grep("deviance", results.row.names)
deviance          <- results$summary[deviance.row.index,]
p.D <- deviance["sd"]^2/2
dic <- deviance["mean"]+p.D
cat("    D-bar: ", deviance["mean"],";  var(dev): ", deviance["sd"]^2,
    "; p.D: ", p.D, "; DIC: ", dic)

# Summary of population sizes. Add pretty printing of results.
cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
temp<- results$summary[ grep("Utot", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\n*** Summary of Total Population Size *** \n")
temp<- results$summary[ grep("Ntot", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)



cat("\n\n\n\n*** Summary of Quantiles of Run Timing *** \n")
cat(    "    This is based on the sample weeks provided and the U[i] values \n")
q <- RunTime(time=time, U=results$sims.list$U, prob=run.prob)
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
results$data <- list( time=time, n1=n1, m2=m2, u2=u2,
                      jump.after=jump.after,
                      bad.n1=bad.n1, bad.m2=bad.m2, bad.u2=bad.u2,
                      logitP.cov=logitP.cov,
                      version=version, date_run=date(),title=title)

return(results)
} # end of function

