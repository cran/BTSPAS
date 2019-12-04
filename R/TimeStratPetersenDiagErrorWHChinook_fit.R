# 2018-12-19 CSJ deprecated use of sampling fraction
# 2018-12-06 CJS saved report to a text connections
# 2018-12-05 CJS converted final spline plot gggplot
# 2018-12-02 CJS converted trace plots to ggplot
# 2018-12-01 CJS changed posterior plot to ggplot
# 2018-11-30 CJS changed acf to ggplot
# 2018-11-28 CJS fixed issued with printing of results getting cut off
# 2018-11-27 CJS Remove refrence to OpenBugs
# 2015-06-10 CJS Fixed error in Bayesian p-value plots. Converted them to ggplot
# 2014-09-01 CJS converted to JAGS
# 2012-08-30 CJS fixed problem in any() and all() in error checking with NAs
# 2011-06-13 CJS add p-values to results
# 2010-11-25 CJS pretty printing of final population estimates
# 2010-09-06 CJS forced input vectors to be vectors
# 2010-08-06 CJS added creation of trace plots to output
# 2010-08-03 CJS added version/date to final data structure
# 2010-03-12 CJS added n.chains etc to calling arguments of the _fit function
# 2009-12-08 CJS added some basic error checking on arguments
# 2009-12-05 CJS added title to argument list
# 2009-12-01 CJS added openbugs/winbugs to argument list; some basic error checking



#' Wrapper (*_fit) to fit the Time Stratified Petersen Estimator
#' with Diagonal Entries and separating Wild from Hatchery Chinook function.
#' 
#' Takes the number of marked fish released, the number of recaptures, and the
#' number of unmarked fish and uses Bayesian methods to fit a fit a spline
#' through the population numbers and a hierarchical model for the trap
#' efficiencies over time.  The output is written to files and an MCMC object
#' is also created with samples from the posterior.
#' 
#' Normally use the *_fit to pass the data to the fitting function.
#' 
#' 
#' @aliases TimeStratPetersenDiagErrorWHChinook_fit TimeStratPetersenDiagErrorWHChinook2_fit
#' 
#' @param title A character string used for a title on reports and graphs
#' @param prefix A character string used as the prefix for created files. All
#' created graph files are of the form prefix-xxxxx.pdf.
#' @param time A numeric vector of time used to label the strata. For example,
#' this could be julian week for data stratified at a weekly level.
#' @param n1 A numeric vector of the number of marked fish released in each
#' time stratum.
#' @param m2 A numeric vector of the number of marked fish from n1 that are
#' recaptured in each time stratum. All recaptures take place within the
#' stratum of release. Use the \code{\link{TimeStratPetersenNonDiagError_fit}}
#' function for cases where recaptures take place outside the stratum of
#' release.
#' @param u2.A A numeric vector of the number of unmarked fish with adipose
#' clips captured in each stratum.
#' @param u2.N A numeric vector of the number of unmarked fish with NO-adipose
#' clips captured in each stratum.
#' @param u2.A.YoY,u2.N.YoY Number of YoY unmarked fish with/without adipose fin clips
#'               All YoY wild fish have NO adipose fin clips; however, hatchery fish are a mixture
#'               of fish with adipose fin clips (a known percentage are marked) and unmarked fish.
#'               So u2.A.YoY MUST be hatchery fish.
#'                  u2.N.YoY is a mixture of wild and hatchery fish.
#' @param u2.A.1,u2.N.1 Number of Age1 unmarked fish with/with out adipose fin clips
#'               All Age1 wild fish have NO adipose fin clips; however, hatchery fish are a mixture
#'               of fish with adipose fin clips (a known percentage are marked) and unmarked fish.
#'               So u2.A.1 MUST be hatchery fish.
#'                  u2.N.1 is a mixture of wild and hatchery fish.
#' @param clip.frac.H.YoY,clip.frac.H.1 Fraction of the YoY hatchery/Age1 (from last year's releases) hatchery fish are clipped?\ (between 0 and 1)
#' @param clip.frac.H A numeric value for the fraction of the hatchery fish
#' that have the adipose fin clipped (between 0 and 1).
#' @param sampfrac \strong{Deprecated}. DO NOT USE ANYMORE. A numeric vector with entries between 0 and 1 indicating
#' what fraction of the stratum was sampled. For example, if strata are
#' calendar weeks, and sampling occurred only on 3 of the 7 days, then the
#' value of \code{sampfrac} for that stratum would be 3/7.
#' @param hatch.after A numeric vector with elements belonging to \code{time}.
#' At which point do hatchery fish arrive? They arrive in the immediate stratum
#' AFTER these entries.
#' @param hatch.after.YoY A numeric vector with elements belonging to
#' \code{time}.  At which point do YoY hatchery fish arrive? They arrive in the
#' immediate stratum AFTER these entries.
#' @param bad.m2 A numeric vector with elements belonging to \code{time}.  In
#' some cases, something goes wrong in the stratum, and the number of recovered
#' fish should be ignored.  For example, poor handling is suspected to induce
#' handling induced mortality in the marked fish and so only very few are
#' recovered. The values of \code{m2} will be set to NA for these strata.
#' @param bad.u2.N A numeric vector with elements belonging to \code{time}.  In
#' some cases, something goes wrong in the stratum, and the number of unmarked
#' fish with NO adipose fin clip should be ignored.
#' @param bad.u2.A.YoY,bad.u2.N.YoY List of julian weeks where the value of u2.A.YoY/u2.N.YoY is suspect. 
#'               These are set to NA prior to the fit.
#' @param bad.u2.A A numeric vector with elements belonging to \code{time}.  In
#' some cases, something goes wrong in the stratum, and the number of unmarked
#' fish with an adipose fin clip should be ignored.
#' @param bad.u2.A.1,bad.u2.N.1   List of julian weeks where the value of u2.A.1/u2.N.1 is suspect. 
#'               These are set to NA prior to the fit.
#' @param logitP.cov A numeric matrix for covariates to fit the
#' logit(catchability). Default is a single intercept, i.e. all strata have the
#' same mean logit(catchability).
#' @param n.chains Number of parallel MCMC chains to fit.
#' @param n.iter Total number of MCMC iterations in each chain.
#' @param n.burnin Number of burn-in iterations.
#' @param n.sims Number of simulated values to keeps for posterior
#' distribution.
#' @param tauU.alpha One of the parameters along with \code{tauU.beta} for the
#' prior for the variance of the random noise for the smoothing spline.
#' @param tauU.beta One of the parameters along with \code{tauU.alpha} for the
#' prior for the variance of the random noise for the smoothing spline.
#' @param taueU.alpha One of the parameters along with \code{taueU.beta} for
#' the prior for the variance of noise around the spline.
#' @param taueU.beta One of the parameters along with \code{taueU.alpha} for
#' the prior for the variance of noise around the spline.
#' @param mu_xiP One of the parameters for the prior for the mean of the
#' logit(catchability) across strata
#' @param tau_xiP One of the parameter for the prior for the mean of the
#' logit(catchability) across strata
#' @param tauP.alpha One of the parameters for the prior for the variance in
#' logit(catchability) among strata
#' @param tauP.beta One of the parameters for the prior for the variance in
#' logit(catchability) among strata
#' @param run.prob Numeric vector indicating percentiles of run timing should
#' be computed.
#' @param debug Logical flag indicating if a debugging run should be made. In
#' the debugging run, the number of samples in the posterior is reduced
#' considerably for a quick turn around.
#' @param debug2 Logical flag indicated if additional debugging information is
#' produced. Normally the functions will halt at \code{browser()} calls to
#' allow the user to peek into the internal variables. Not useful except to
#' package developers.
#' @param InitialSeed Numeric value used to initialize the random numbers used
#' in the MCMC iterations.
#' @param save.output.to.files Should the plots and text output be save to the files
#' in addition to being stored in the MCMC object? 

#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' @template author 
#' @template references
#' @keywords ~models ~smooth
#' @examples
#'  
#' ##---- See the vignettes for examples on how to run this analysis.
#' 
#' @export TimeStratPetersenDiagErrorWHChinook_fit

TimeStratPetersenDiagErrorWHChinook_fit<- 
       function( title="TSPDE-WHChinook", prefix="TSPDE-WHChinook-", 
                 time, n1, m2, u2.A, u2.N, clip.frac.H, sampfrac=rep(1,length(u2.A)), 
                 hatch.after=NULL, 
                 bad.m2=c(), bad.u2.A=c(), bad.u2.N=c(),
                 logitP.cov=rep(1,length(n1)),
                 n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
                 tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
                 mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)), 
                 tau_xiP=1/var(logit((m2+.5)/(n1+1)), na.rm=TRUE), 
                 tauP.alpha=.001, tauP.beta=.001,
                 run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
                 debug=FALSE, debug2=FALSE,
                 InitialSeed=ceiling(runif(1,min=0, max=1000000)),
                 save.output.to.files=TRUE) {
# Fit a Time Stratified Petersen model with diagonal entries and with smoothing on U allowing for random error,
# covariates for the the capture probabilities, and separating the wild vs hatchery fish
# The "diagonal entries" implies that no marked fish are recaptured outside the (time) stratum of release
#
   version <- '2020-01-01'
   options(width=200)

# Input parameters are
#    title - title for the analysis
#    prefix - prefix used for files created with the analysis results
#             this should be in standard Window's format, eg. JC-2002-ST-TSPDE
#             to which is appended various suffixes for plots etc
#    time   - vector of stratum numbers. For example, 9:38 would indicate that the
#             Trinity River system sampled weeks 9 to 38. If some values are omitted
#             e.g. time=10 not present, this indicates sampling did not take place this
#             week. The data are expanded and interpolation for the missing week takes place
#    n1, m2 - the input data consisting of fish marked and released and then recaptured.
#                 The n1 and m2 are used to calibrate the trap
#    u2.A     - number of unmarked fish with adipose fin clips
#    u2.N     - number of unmarked fish with NO adipose fin clips
#               All wild fish have NO adipose fin clips; however, hatchery fish are a mixture
#               of fish with adipose fin clips (a known percentage are marked) unmarked fish.
#               So u2.A MUST be hatchery fish.
#                  u2.N is a mixture of wild and hatchery fish.
#    clip.frac.H - what fraction of the hatchery fish are clipped?
#    sampfrac - Deprecated. DO NOT USE ANYMORE. sampling fraction to adjust for how many days of the week was the trap operating
#              This is expressed as fraction i.e. 3 days out of 7 is expressed as 3/7=.42 etc.
#              If the trap was operating ALL days, then the SampFrac = 1. It is possible for the sampling
#              fraction to be > 1 (e.g. a mark is used for 8 days instead of 7. The data are adjusted
#              back to a 7 day week as well.
#    hatch.after - julian week AFTER which hatchery fish are released 
#    bad.m2  - list of julian numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to OpenBugs
#    bad.u2.A - list of julian weeks where the value of u2.A is suspect. 
#               These are set to NA prior to the call to OpenBugs
#    bad.u2.N - list of julian weeks where the value of u2.N is suspect.
#               These are set to NA prior to the call to OpenBugs
#    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
#               for the first element of the covariance matrix and 0 for the rest of the covariates.
#               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
#               if you have any covariates. The default, if not specified, is a constant (the mean logit)
#    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
#    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline 
#    mu_xiP, tau_xiP         - parameters for the prior on mean logit(P)'s [The intercept term]
#                              The other covariates are assigned priors of a mean of 0 and a variance of 1000
#    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
#    run.prob  - percentiles of run timing wanted 
#    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data
#             and OpenBUGS will run and stop waiting for your to exit and complete

# force the input vectors to be vectors
time      <- as.vector(time)
n1        <- as.vector(n1)
m2        <- as.vector(m2)
u2.A      <- as.vector(u2.A)
sampfrac  <- as.vector(sampfrac)

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time all match
if(var(c(length(n1),length(m2),length(u2.A),length(u2.N),length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of n1, m2, u2.A, u2.N, sampfrac, time must all be equal. They are:",
        length(n1)," ",length(m2)," ",length(u2.A)," ",length(u2.N)," ",length(sampfrac)," ",length(time),"\n")
   return()}
if(length(logitP.cov) %% length(n1) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of n1 etc They are:",
        length(n1)," ",length(logitP.cov)," ",paste(dim(logitP.cov),collapse=","),"\n")
   return()}

#  2. Check that m2<= n1
if(any(m2>n1, na.rm=TRUE)){
   cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:",
       paste(n1,collapse=","),"\n m2:",
       paste(m2,collapse=","),"\n")
   return()}

#  3. Elements of bad.m2, bad.u2.A, and bad.u2.N, and hatch.after must belong to time
if(!all(bad.m2 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",
       paste(bad.m2,collapse=","),"\n Strata identifiers are \n time:",
       paste(time  ,collapse=","), "\n")
   return()}
if(!all(bad.u2.A %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.A must be elements of strata identifiers. You entered \n bad.u2.A:",
       paste(bad.u2.A,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,    collapse=","), "\n")
   return()}
if(!all(bad.u2.N %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.N must be elements of strata identifiers. You entered \n bad.u2.N:",
       paste(bad.u2.N, collapse=","),"\n Strata identifiers are \n time:",
       paste(time,     collapse=","), "\n")
   return()}
if(!all(hatch.after %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** hatch.after must be elements of strata identifiers. You entered \n hatch.after:",
       paste(hatch.after,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,       collapse=","), "\n")
   return()}

# Deprecation of sampling fraction.
if(any(sampfrac != 1)){
  cat("***** ERROR ***** Sampling fraction is deprecated for any values other than 1. DO NOT USE ANYMORE. ")
  return()
}

results.filename <- paste(prefix,"-results.txt",sep="")   

stdout <- vector('character')
report <- textConnection('stdout', 'wr', local = TRUE)
sink(report)

cat(paste("Time Stratified Petersen with Diagonal recaptures, error in smoothed U, separating wild and hatchery fish - ", date()))
cat("\nVersion:", version)

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** \n")
temp<- cbind(time, n1, m2, u2.A, u2.N, round(sampfrac,digits=2), logitP.cov)
colnames(temp)<- c('time', 'n1','m2','u2.A', 'u2.N','SampFrac', paste("logitPcov[", 1:ncol(as.matrix(logitP.cov)),"]",sep="") )
print(temp) 
cat("\n\n")
cat("Hatchery fish are released AFTER strata: ", hatch.after,"\n\n")
cat("Hatchery fish are clipped at a rate of :", clip.frac.H,"\n\n")
cat("The following strata had m2   set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.A set to missing: ", 
     if(length(bad.u2.A)>0){bad.u2.A} else {" NONE"}, "\n")
cat("The following strata had u2.N set to missing: ", 
     if(length(bad.u2.N)>0){bad.u2.N} else {" NONE"}, "\n")



# Pooled Petersen estimator over ALL of the data including when no releases take place, bad m2, bad.u2.A or bad.u2.N values.
cat("\n\n*** Pooled Petersen Estimate based on pooling over ALL strata adjusting for sampling fraction***\n\n")
cat("Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",
     sum(u2.A/sampfrac, na.rm=TRUE)+sum(u2.N/sampfrac, na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.A/sampfrac, na.rm=TRUE)+sum(u2.N/sampfrac, na.rm=TRUE))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

# estimate for clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.A=",  sum(u2.A/sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum(u2.A/sampfrac, na.rm=TRUE))
cat("Est U.H(total) ", format(round(pp$U.est)/clip.frac.H,big.mark=","),
    "  (SE ",          format(round(pp$U.se) /clip.frac.H,big.mark=","), ")\n")
cat("Est N.H(total) ", format(round(pp$N.est)/clip.frac.H,big.mark=","),
    "  (SE ",          format(round(pp$N.se) /clip.frac.H,big.mark=","), ")\n\n\n")


# estimate for wild YoY fish found by subtraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.W=",  sum((u2.N+u2.A-u2.A/clip.frac.H)/sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum((u2.N+u2.A-u2.A/clip.frac.H)/sampfrac, na.rm=TRUE))
cat("Est U.W(total) ", format(round(pp$U.est),big.mark=","),
    "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE\n")
cat("Est N.W(total) ", format(round(pp$N.est),big.mark=","),
    "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE\n\n\n")



# Obtain the Pooled Petersen estimator without excluding bad.m2, bad.u2.A, or bad.u2.N values but after removing 0 or NA values
select <- (n1>0) & (!is.na(n1)) & (!is.na(m2)) & (!is.na(u2.A)) & (!is.na(u2.N))
cat("\n\n*** Pooled Petersen Estimate prior to excluding bad m2, u2.A, or u2.N values  ***\n\n")
cat("The following strata are excluded because n1=0 or NA values in m2, u2.A, u2.N :", time[!select],"\n\n")

temp.n1 <-      n1      [select]
temp.m2 <-      m2      [select]
temp.u2.A <-    u2.A    [select]
temp.u2.N <-    u2.N    [select]
temp.sampfrac<- sampfrac[select]

cat("Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.A/temp.sampfrac+temp.u2.N/temp.sampfrac),"\n\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.A/temp.sampfrac+temp.u2.N/temp.sampfrac))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

# estimate for clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.A=",  sum(temp.u2.A/temp.sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum(temp.u2.A/temp.sampfrac, na.rm=TRUE))
cat("Est U.H(total) ", format(round(pp$U.est)/clip.frac.H,big.mark=","),
    "  (SE ",          format(round(pp$U.se) /clip.frac.H,big.mark=","), ")\n")
cat("Est N.H(total) ", format(round(pp$N.est)/clip.frac.H,big.mark=","),
    "  (SE ",          format(round(pp$N.se) /clip.frac.H,big.mark=","), ")\n\n\n")


# estimate for wild YoY fish
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.W=",  sum((temp.u2.N+temp.u2.A-temp.u2.A/clip.frac.H)/temp.sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum((temp.u2.N+temp.u2.A-temp.u2.A/clip.frac.H)/temp.sampfrac, na.rm=TRUE))
cat("Est U.W(total) ", format(round(pp$U.est),big.mark=","),
    "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE \n")
cat("Est N.W(total) ", format(round(pp$N.est),big.mark=","),
    "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE \n\n\n")



# Obtain the Pooled Petersen estimator after fixup of bad.m2, bad.u2.A, and bad.u2.N values
temp.m2 <- m2
index.bad.m2 <- as.vector((1:length(time)) %*% outer(time,bad.m2,"=="))
temp.m2[index.bad.m2] <- NA
temp.u2.A <- u2.A
index.bad.u2.A <- as.vector((1:length(time)) %*% outer(time,bad.u2.A,"=="))
temp.u2.A[index.bad.u2.A] <- NA
temp.u2.N <- u2.A
index.bad.u2.N <- as.vector((1:length(time)) %*% outer(time,bad.u2.N,"=="))
temp.u2.N[index.bad.u2.N] <- NA

select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2)) & (!is.na(temp.u2.A) & (!is.na(temp.u2.N)) )
cat("\n\n*** Pooled Petersen Estimate after removing bad m2, u2.A, and u2.N values adjusting for sampling fraction  ***\n\n")
cat("The following strata had m2   set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.A set to missing: ", 
     if(length(bad.u2.A)>0){bad.u2.A} else {" NONE"}, "\n")
cat("The following strata had u2.N set to missing: ", 
     if(length(bad.u2.N)>0){bad.u2.N} else {" NONE"}, "\n")
cat("The following strata are excluded because n1=0 or NA values in m2, u2.A, or u2.N:", time[!select],"\n\n")

temp.n1       <- n1      [select]
temp.m2       <- m2      [select]
temp.u2.A     <- u2.A    [select]
temp.u2.N     <- u2.N    [select]
temp.sampfrac <- sampfrac[select]

cat("Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.A/temp.sampfrac+temp.u2.N/temp.sampfrac),"\n\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.A/temp.sampfrac+temp.u2.N/temp.sampfrac))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")


# estimate for clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.A=",  sum(temp.u2.A/temp.sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum(temp.u2.A/temp.sampfrac, na.rm=TRUE))
cat("Est U.H(total) ", format(round(pp$U.est)/clip.frac.H,big.mark=","),
    "  (SE ",          format(round(pp$U.se) /clip.frac.H,big.mark=","), ")\n")
cat("Est N.H(total) ", format(round(pp$N.est)/clip.frac.H,big.mark=","),
    "  (SE ",          format(round(pp$N.se) /clip.frac.H,big.mark=","), ")\n\n\n")

# estimate for wild YoY fish
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.W=",  sum((temp.u2.N+temp.u2.A-temp.u2.A/clip.frac.H)/temp.sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum((temp.u2.N+temp.u2.A-temp.u2.A/clip.frac.H)/temp.sampfrac, na.rm=TRUE))
cat("Est U.W(total) ", format(round(pp$U.est),big.mark=","),
    "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE\n")
cat("Est N.W(total) ", format(round(pp$N.est),big.mark=","),
    "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE\n\n\n")



# Obtain Stratified-Petersen estimator for each stratum prior to removing bad m2 values
cat("*** Stratified Petersen Estimator for each stratum PRIOR to removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- (u2.A + u2.N)/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','(u2.A+u2.N)*adj', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ",        format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Hatchery YoY PRIOR to removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.A/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.A*adj', 'U[i]', 'SE(U[i])')
print(temp)
cat("** Estimates not adjusted for clip fraction above \n")
cat("Est U.H(total) ", format(round(sum(sp$U.est, na.rm=TRUE)/clip.frac.H),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))/clip.frac.H), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Wild YoY  PRIOR to removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- pmax(0,(u2.N+u2.A-u2.A/clip.frac.H)/sampfrac)
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W-est', 'U[i]', 'SE(U[i])')
print(temp)
cat("Est U.W(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ") APPROXIMATE\n\n\n")





# Obtain Stratified-Petersen estimator for each stratum after removing bad m2, u2.A, or u2.N values
cat("*** Stratified Petersen Estimator for each stratum AFTER removing bad m2, u2.A, u2.N values adjusting for sampling fraction***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.m2[index.bad.m2] <- NA
temp.u2.A <- u2.A
temp.u2.A[index.bad.u2.A] <- NA
temp.u2.N <- u2.N
temp.u2.N[index.bad.u2.N] <- NA
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.A/sampfrac+temp.u2.N/sampfrac)
temp <- cbind(time, temp.n1, temp.m2, (temp.u2.A+temp.u2.N)/sampfrac, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','(u2.a+u2.N)*adj', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum YoY hatchery PRIOR to removing bad m2 values after adjusting for sampling fration ***\n\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.A/sampfrac)
temp <- cbind(time, temp.n1, temp.m2, round(temp.u2.A/sampfrac), round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.A*adj', 'U[i]', 'SE(U[i])')
print(temp)
cat("** Estimates not adjusted for clip fraction above \n")
cat("Est U.H(total) ", format(round(sum(sp$U.est, na.rm=TRUE)/clip.frac.H),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))/clip.frac.H), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Wild YoY  PRIOR to removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2.W <- pmax(0,(temp.u2.N+temp.u2.A-temp.u2.A/clip.frac.H)/sampfrac)
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.W)
temp <- cbind(time, temp.n1, temp.m2, round(temp.u2.W), round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W-est', 'U[i]', 'SE(U[i])')
print(temp)
cat("Est U.W(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ") APPROXIMATE\n\n\n")









# Test if pooling can be done
cat("*** Test if pooled Petersen is allowable. [Check if marked fractions are equal] ***\n\n")
select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2)) 
temp.n1 <- n1[select]
temp.m2 <- m2[select]
test <- TestIfPool( temp.n1, temp.m2)
cat("(Large Sample) Chi-square test statistic ", test$chi$statistic," has p-value", test$chi$p.value,"\n\n")
temp <- cbind(time[select],test$chi$observed, round(test$chi$expected,1), round(test$chi$residuals^2,1))
colnames(temp) <- c('time','n1-m2','m2','E[n1-m2]','E[m2]','X2[n1-m2]','X2[m2]')
print(temp)
cat("\n Be cautious of using this test in cases of small expected values. \n\n")



# Fix up any data problems and prepare for the call.
# Notice that for strata entries that are missing any covariate values, only an intercept is added

# Expand the entries in case of missing time entries
new.n1         <- rep(0, max(time)-min(time)+1)
new.m2         <- rep(0, max(time)-min(time)+1)
new.u2.A       <- rep(0, max(time)-min(time)+1)
new.u2.N       <- rep(0, max(time)-min(time)+1)
new.sampfrac   <- rep(0, max(time)-min(time)+1)
new.logitP.cov <- matrix(NA, nrow=max(time)-min(time)+1, ncol=ncol(as.matrix(logitP.cov)))
new.time       <- min(time):max(time)


new.n1  [time-min(time)+1]         <- n1
new.m2  [time-min(time)+1]         <- m2
new.m2  [bad.m2-min(time)+1]       <- NA    # wipe out strata where m2 is known to be bad
new.u2.A[time-min(time)+1]         <- u2.A
new.u2.A[bad.u2.A-min(time)+1]     <- NA    # wipe out strata where u2.A is known to be bad
new.u2.N[time-min(time)+1]         <- u2.N
new.u2.N[bad.u2.N-min(time)+1]     <- NA    # wipe out strata where u2.N is known to be bad
new.sampfrac[time-min(time)+1]   <- sampfrac
new.logitP.cov[time-min(time)+1,]<- as.matrix(logitP.cov)
new.logitP.cov[ is.na(new.logitP.cov[,1]), 1] <- 1  # insert a 1 into first columns where not specified
new.logitP.cov[ is.na(new.logitP.cov)] <- 0         # other covariates are forced to zero not in column 1


# Check for and fix problems with the data
# If n1=m2=0, then set n1 to 1, and set m2<-NA
new.m2[new.n1==0] <- NA
new.n1[new.n1==0] <- 1

# Adjust data when a stratum has less than 100% sampling fraction to "estimate" the number
# of unmarked fish that were captured. It is not necessary to adjust the n1 and m2 values 
# as these are used ONLY to estimate the capture efficiency. 
# In reality, there should be a slight adjustment
# to the precision to account for this change, but this is not done.
# Similarly, if the sampling fraction is more than 1, the adjustment forces the total unmarked catch back to a single week.
new.u2.A <- round(new.u2.A/new.sampfrac)
new.u2.N <- round(new.u2.N/new.sampfrac)

# Adjust for strata where sampling fraction=0. On these strata
# u2.A and u2.N is set to NA so that there is NO information on U2 for this stratum
new.u2.A[new.sampfrac<.001] <- NA
new.u2.N[new.sampfrac<.001] <- NA

# Print out the revised data
hatch.indicator <- rep('   ', max(time)-min(time)+1)
hatch.indicator[hatch.after-min(time)+1]<- '***'

cat("\n\n*** Revised data *** \n")
temp<- data.frame(time=new.time, n1=new.n1, m2=new.m2, u2.A=new.u2.A, u2.N=new.u2.N, 
       sampfrac=round(new.sampfrac,digits=2), new.logitP.cov=new.logitP.cov, 
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
cat("   Parameters for prior on beta.logitP[1] (intercept) (mean, 1/var):", round(mu_xiP,3), round(tau_xiP,5),
    " which corresponds to a median P of ", round(expit(mu_xiP),3), "\n")
cat("   Parameters for prior on tauP (residual variance of logit(P) after adjusting for covariates: ",tauP.alpha, tauP.beta, 
    " which corresponds to a mean/std dev of 1/var of:",
    round(tauP.alpha/tauP.beta,2),round(sqrt(tauP.alpha/tauP.beta^2),2),"\n")

cat("\n\nInitial seed for this run is: ",InitialSeed, "\n")

sink()

if (debug2) {
   cat("\nprior to formal call to TimeStratPetersenDiagErrorWHChinook\n")
   browser()
}


if (debug) 
   {results <- TimeStratPetersenDiagErrorWHChinook(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2.A=new.u2.A, u2.N=new.u2.N, 
            hatch.after=hatch.after-min(time)+1, clip.frac.H=clip.frac.H,
            logitP.cov=new.logitP.cov,
            n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,  # set to low values for debugging only
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug,InitialSeed=InitialSeed,
            save.output.to.files=save.output.to.files)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenDiagErrorWHChinook(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2.A=new.u2.A, u2.N=new.u2.N, 
            hatch.after=hatch.after-min(time)+1, clip.frac.H=clip.frac.H,
            logitP.cov=new.logitP.cov,
            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            InitialSeed=InitialSeed,
            save.output.to.files=save.output.to.files)
   }

# Now to create the various summary tables of the results

if (debug2) {
   cat("\nAfter formal call to TimeStratPetersenDiagErrorWHChinook\n")
   browser()
}

  # A plot of the observered log(U) on the log scale, and the final mean log(U)
  # A plot of the observered log(U) on the log scale, and the final mean log(U)
  # Create the data frame needed for ggplot. 
  # In the diagonal case, time, n1, m2, u2 are the same length
  Nstrata <- length(n1)
  plot.df   <- data.frame(time =new.time)
  
  # adjust the u2 for the clipping fractions
  plot.df$n1   <- new.n1
  plot.df$m2   <- new.m2
  plot.df$u2.H <- new.u2.A/clip.frac.H  # only a portion of the hatchery fish are clipped
  plot.df$u2.N <- new.u2.N
  plot.df$u2.W <- pmax(plot.df$u2.N - plot.df$u2.H*(1-clip.frac.H),0) # subtract the guestimated number of hatchery fish
  plot.df$u2.H[is.na(plot.df$u2.H)] <- 1  # in case of missing values
  plot.df$u2.W[is.na(plot.df$u2.W)] <- 1  # in case of missing values

  avg.P <- sum(plot.df$m2,na.rm=TRUE)/sum(plot.df$n1, na.rM=TRUE)
  plot.df$logUguess.W <- log(pmax((plot.df$u2.W+1)*(plot.df$n1+2)/(plot.df$m2+1), plot.df$u2.W/avg.P, na.rm=TRUE) ) # try and keep Uguess larger than observed values
  plot.df$logUguess.H <- log(pmax((plot.df$u2.H+1)*(plot.df$n1+2)/(plot.df$m2+1), plot.df$u2.H/avg.P, na.rm=TRUE) )
  plot.df$logUguess.H[1:(hatch.after-min(time))] <- NA   # no hatchery fish prior to release from hatchery
  
# extract the fitted U values for W (wild) and H (hatchery)
  results.row.names <- rownames(results$summary)
  etaU.W.row.index    <- grep("etaU.W", results.row.names)
  etaU.W <- results$summary[etaU.W.row.index,]
  plot.df$logU.W =etaU.W[,"mean"]
  plot.df$lcl.W  =etaU.W[,"2.5%"]
  plot.df$ucl.W  =etaU.W[,"97.5%"]
  
  etaU.H.row.index    <- grep("etaU.H", results.row.names)
  etaU.H <- results$summary[etaU.H.row.index,]
  plot.df$logU.H =etaU.H[,"mean"]
  plot.df$lcl.H  =etaU.H[,"2.5%"]
  plot.df$ucl.H  =etaU.H[,"97.5%"]
  plot.df$logU.H [1:(hatch.after - min(time)+1)] <- NA # no hatchery fish until release at hatch.after
  plot.df$lcl.H  [1:(hatch.after - min(time)+1)] <- NA
  plot.df$ucl.H  [1:(hatch.after - min(time)+1)] <- NA

# extract the spline values for W (wild) and H (hatchery) fish
  logUne.W.row.index <- grep("logUne.W", results.row.names)
  plot.df$spline.W  <- results$summary[logUne.W.row.index,"mean"]
  logUne.H.row.index <- grep("logUne.H", results.row.names)
  plot.df$spline.H  <- results$summary[logUne.H.row.index,"mean"]
  plot.df$spline.H [1:(hatch.after - min(time)+1)] <- NA # no hatchery fish until release at hatch.after

fit.plot <- ggplot(data=plot.df, aes_(x=~new.time))+
   ggtitle(title, subtitle="Fitted spline curve to raw U.W[i] U.H[i] with 95% credible intervals")+
   geom_point(aes_(y=~logUguess.W), color="red",  shape="w")+  # guesses for wild file
   geom_point(aes_(y=~logUguess.H), color="green", shape="h")+  # guesses for hatchery fish
   xlab("Time Index\nFitted/Smoothed/Raw values plotted for W(black) and H(blue)")+ylab("log(U[i])")+
   geom_point(aes_(y=~logU.W), color="black", shape=19)+
   geom_line (aes_(y=~logU.W), color="black")+
   geom_errorbar(aes_(ymin=~lcl.W, ymax=~ucl.W), width=.1)+
   geom_line(aes_(y=~spline.W),linetype="dashed") +  
   geom_point(aes_(y=~logU.H), color="blue", shape=19)+
   geom_line (aes_(y=~logU.H), color="blue")+
   geom_errorbar(aes_(ymin=~lcl.H, ymax=~ucl.H), width=.1, color="blue")+
   geom_line(aes_(y=~spline.H),linetype="dashed",color="blue")+
   ylim(c(-2,NA))
 
if(save.output.to.files)ggsave(plot=fit.plot, filename=paste(prefix,"-fit.pdf",sep=""), height=6, width=10, units="in")
results$plots$fit.plot <- fit.plot

# Plot the logitP over time
logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=new.m2, u2=new.u2.A+new.u2.N,  logitP.cov=new.logitP.cov, results=results)
if(save.output.to.files)ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in", dpi=300)
results$plots$logitP.plot <- logitP.plot


# Look at autocorrelation function for Utot.W and Utot.H
mcmc.sample1<- data.frame(parm="Utot.W", sample=results$sims.matrix[,"Utot.W"], stringsAsFactors=FALSE)
mcmc.sample2<- data.frame(parm="Utot.H", sample=results$sims.matrix[,"Utot.H"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2)
acf.Utot.plot <- plot_acf(mcmc.sample)
if(save.output.to.files)ggsave(plot=acf.Utot.plot, filename=paste(prefix,"-Utot-acf.pdf",sep=""), height=4, width=6, units="in")
results$plots$acf.Utot.plot <- acf.Utot.plot


# Look at the shape of the posterior distribution
mcmc.sample1<- data.frame(parm="Utot.W", sample=results$sims.matrix[,"Utot.W"], stringsAsFactors=FALSE)
mcmc.sample2<- data.frame(parm="Utot.H", sample=results$sims.matrix[,"Utot.H"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2)
post.Utot.plot <- plot_posterior(mcmc.sample)
if(save.output.to.files)ggsave(plot=post.Utot.plot, filename=paste(prefix,"-Utot-posterior.pdf",sep=""), height=4, width=6, units="in")
results$plots$post.Utot.plot <- post.Utot.plot


#save the Bayesian predictive distribution (Bayesian p-value plots)
#browser()
discrep <-PredictivePosterior.TSPDE.WHCH (time, new.n1, new.m2, new.u2.A, new.u2.N, clip.frac.H, 
          expit(results$sims.list$logitP), round(results$sims.list$U.W), 
          round(pmax(results$sims.list$U.H,0)),
          hatch.after)
gof <- PredictivePosteriorPlot.TSPDE.WHCH (discrep)
if(save.output.to.files)ggsave(gof[[1]],filename=paste(prefix,"-GOF.pdf",sep=""),   height=8, width=8, units="in", dpi=300 )
results$plots$gof.plot <- gof

# create traceplots of logU, U, and logitP (along with R value) to look for non-convergence
# the plot_trace will return a list of plots (one for each page as needed)
varnames <- names(results$sims.array[1,1,])  # extract the names of the variables 

# Trace plots of logitP
trace.plot <- plot_trace(title=title, results=results, parms_to_plot=varnames[grep("^logitP", varnames)])
if(save.output.to.files){
   pdf(file=paste(prefix,"-trace-logitP.pdf",sep=""))
   l_ply(trace.plot, function(x){plot(x)})
   dev.off()
}
results$plots$trace.logitP.plot <- trace.plot

# now for the traceplots of logU (etaU), Utot, and Ntot
trace.plot <- plot_trace(title=title, results=results, parms_to_plot=varnames[c(grep("Utot",varnames), grep("Ntot",varnames), grep("^etaU", varnames))])
if(save.output.to.files){
   pdf(file=paste(prefix,"-trace-logU.pdf",sep=""))
   l_ply(trace.plot, function(x){plot(x)})
   dev.off()
}
results$plots$trace.logU.plot <- trace.plot


sink(report, append=TRUE)
# What was the initial seed
cat("\n\n*** Initial Seed for this run ***: ", results$Seed.initial,"\n")

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

# Summary of population sizes. Add pretty printing.
cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
cat("Wild\n")
temp<- results$summary[ grep("Utot.W", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nHatchery\n")
temp<- results$summary[ grep("Utot.H", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nTotal\n")
temp<- results$summary[ rownames(results$summary) == "Utot",]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)



#browser()
cat("\n\n\n\n*** Summary of Quantiles of Run Timing.Wild *** \n")
cat(    "    This is based on the sample weeks provided and the U.W[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n*** Summary of Quantiles of Run Timing.Hatchery *** \n")
cat(    "    This is based on the sample weeks provided and the U.H[i] values \n") 
q <- RunTime(time=time[time>hatch.after], U=results$sims.list$U.H[,time>hatch.after], prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n")
cat(paste("*** end of fit *** ", date()))

sink()

# save the report to a files?
if(save.output.to.files)writeLines(stdout, results.filename)
results$report <- stdout


# add some of the raw data to the bugs object for simplicity in referencing it later
results$data <- list( time=time, n1=n1, m2=m2, u2.A=u2.A, u2.N=u2.N, clip.frac.H=clip.frac.H,
                      sampfrac=sampfrac, 
                      hatch.after=hatch.after, 
                      bad.m2=bad.m2, bad.u2.A=bad.u2.A, bad.u2.N=bad.u2.N, 
                      logitP.cov=logitP.cov,
                      version=version, date_run=date(),
                      title=title)

return(results)
} # end of function
