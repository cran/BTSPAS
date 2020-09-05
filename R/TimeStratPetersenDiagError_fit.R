# 2018-12-19 CJS Deprecated use of sampling.fraction
# 2018-12-15 CJS Added ability to fix some logitP values
# 2018-12-02 CJS Convert trace plots to ggplot
# 2018-12-01 CJS COnvert posterior plots to ggplot
# 2018-11-30 CJS Convert acf plot to ggplot
# 2018-11-28 CJS Fixed problem where printing results got cutoff
# 2018-11-25 CJS Removed all OpenBugs stuff
# 2015-06-10 CJS Change gof plot to ggplot()
# 2014-09-01 CJS Converted to JAGS engine from OpenBugs
# 2013-09-04 CJS Initialized n1, m2, u2 that are NA to sensible values.
#                Removed references to WinBugs
# 2012-08-30 CJS fixed problem in error checking in any() function that includes missing values
# 2011-06-13 CJS inserted the bayesian p-values in the results 
# 2010-11-29 CJS added code for bad.n1 to the call.
# 2010-11-25 CJS added code for bad.u2 to the call. Simplified the two pooled and simple Petersen estimates
# 2010-11-25 CJS pretty printing for estimates of Utot, Ntot
# 2010-09-06 CJS forced time, n1, m2, u2, sampfrac to be vectors
# 2010-08-04 CJS added traceplots of logitP, logU, Utot, and Ntot to help diagnose non-mixing
# 2010-08-04 CJS added version/date to results data structure
# 2009-12005 CJS added title to argument list
# 2009-12-01 CJS Added some basic error checking; added OPENBUGS/WINBUGS to argument list


#' Wrapper (*_fit) to call the Time Stratified Petersen Estimator
#' with Diagonal Entries function.
#' 
#' Takes the number of marked fish released, the number of recaptures, and the
#' number of unmarked fish and uses Bayesian methods to fit a fit a spline
#' through the population numbers and a hierarchical model for the trap
#' efficiencies over time.  The output is written to files and an MCMC object
#' is also created with samples from the posterior.
#' 
#' Normally, the wrapper (*_fit) function is called which then calls the
#' fitting routine.
#' 
#' Use the \code{\link{TimeStratPetersenNonDiagError_fit}} function for cases
#' where recaptures take place outside the stratum of release.
#' 
#' 
#' @aliases TimeStratPetersenDiagError_fit 
#' @param title A character string used for a title on reports and graphs
#' @param prefix A character string used as the prefix for created files. All
#' created graph files are of the form prefix-xxxxx.pdf.
#' @param time A numeric vector of time used to label the strata. For example,
#' this could be julian week for data stratified at a weekly level.
#' @param n1 A numeric vector of the number of marked fish released in each
#' time stratum.
#' @param m2 A numeric vector of the number of marked fish from n1 that are
#' recaptured in each time stratum. All recaptures take place within the
#' stratum of release.
#' @param u2 A numeric vector of the number of unmarked fish captured in each
#' stratum. These will be expanded by the capture efficiency to estimate the
#' population size in each stratum.
#' @param sampfrac \strong{Deprecated} because it really doesn't work as intended.
#' A numeric vector with entries between 0 and 1 indicating
#' what fraction of the stratum was sampled. For example, if strata are
#' calendar weeks, and sampling occurred only on 3 of the 7 days, then the
#' value of \code{sampfrac} for that stratum would be 3/7. 
#' @param jump.after A numeric vector with elements belonging to \code{time}.
#' In some cases, the spline fitting the population numbers should be allowed
#' to jump. For example, the population size could take a jump when hatchery
#' released fish suddenly arrive at the trap. The jumps occur AFTER the strata
#' listed in this argument.
#' @param bad.n1 A numeric vector with elements belonging to \code{time}. In
#' some cases, something goes wrong in the stratum, and the number of marked
#' fish released should be ignored.  The values of \code{m2} for this stratum
#' will also be set to NA for these strata.
#' @param bad.m2 A numeric vector with elements belonging to \code{time}. In
#' some cases, something goes wrong in the stratum, and the number of recovered
#' fish should be ignored. For example, poor handling is suspected to induce
#' handling induced mortality in the marked fish and so only very few are
#' recovered. The values of \code{m2} will be set to NA for these strata.
#' @param bad.u2 A numeric vector with elements belonging to \code{time}. In
#' some cases, something goes wrong in the stratum, and the number of unmarked
#' fish should be ignored. For example, the trap didn't work properly in this
#' stratum.  The values of \code{u2} will be set to NA for these strata.
#' @param logitP.cov A numeric matrix for covariates to fit the
#' logit(catchability). Default is a single intercept, i.e. all strata have the
#' same mean logit(catchability).
#' @param logitP.fixed A numeric vector (could be null) of the time strata
#' where the logit(P) would be fixed. Typically, this is used when the capture
#' rates for some strata are 0 and logit(P) is set to -10 for these strata. The
#' fixed values are given in \code{logitP.fixed.values}
#' @param logitP.fixed.values A numerical vector (could be null) of the fixed
#' values for logit(P) at strata given by logitP.fixed. Typically this is used
#' when certain strata have a 0 capture rate and the fixed value is set to -10
#' which on the logit scale gives p[i] essentially 0. Don't specify values such
#' as -50 because numerical problems could occur in WinBugs/OpenBugs.

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

#' 
#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' @template author
#' @template references
#' @keywords ~models ~smooth
#' @export TimeStratPetersenDiagError_fit
#' 
#' 

TimeStratPetersenDiagError_fit <-
  function( title="TSDPE", prefix="TSPDE-", 
           time, n1, m2, u2, sampfrac=rep(1,length(u2)), 
           jump.after=NULL, 
           bad.n1=c(), bad.m2=c(), bad.u2=c(),
           logitP.cov=rep(1,length(n1)),
           logitP.fixed=NULL, logitP.fixed.values=NULL,
           n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
           tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
           mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),
           tau_xiP=1/var(logit((m2+.5)/(n1+1)),na.rm=TRUE), 
           tauP.alpha=.001, tauP.beta=.001,
           run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
           debug=FALSE, debug2=FALSE,
           InitialSeed=ceiling(runif(1,min=0, max=1000000)),
           save.output.to.files=TRUE) {
    
# Fit a Time Stratified Petersen model with diagonal entries and with smoothing on U allowing for random error
# The "diagonal entries" implies that no marked fish are recaptured outside the (time) stratum of release
#
   version <- '2020-09-01'
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
#    n1, m2, u2 - the input data consisting of fish marked and released, recapture, and unmarked captured
#    sampfrac - sampling fraction to adjust for how many days of the week was the trap operating
#              This is expressed as fraction i.e. 3 days out of 7 is expressed as 3/7=.42 etc.
#              If the trap was operating ALL days, then the SampFrac = 1. It is possible for the sampling
#              fraction to be > 1 (e.g. a mark is used for 8 days instead of 7. The data are adjusted
#              back to a 7 day week as well.
#    jump.after - in some cases, a single spline is still not flexible enough to cope with rapid
#                 changes in the run curve. For example, in the Trinity River project, a larger
#                 hatchery release occurs around stratum 14. This is a vector indicating the
#                 strata AFTER which the spline curve is allowed to jump.
#    bad.n1  - list of stratum numbers where the value of n1 is suspect.
#              Note that if the value of n1 is suspect, the value of m2 is also likely suspect.
#              These are replaced by the value of (1,0). We need to specify a value of 1 for bad.n1 values
#              because OpenBugs gets upset with n1=0 or n1=NA.
#    bad.m2  - list of stratum numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to OpenBugs
#    bad.u2  - list of stratum numbers where the value of u2 is suspect.
#              For example, the trap may not be operating completely for some strata, or there was miss counting?
#              These are set to NA prior to the call to OpenBugs
#    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
#               for the first element of the covariance matrix and 0 for the rest of the covariates.
#               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
#               if you have any covariates. The default, if not specified, is a constant (the mean logit)
#    logitP.fixed, logitP.fixed.values - if you are fixing any of the logit P and at what values. 
#    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
#    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline 
#    mu_xiP, tau_xiP         - parameters for the prior on mean logit(P)'s [The intercept term]
#                              The other covariates are assigned priors of a mean of 0 and a variance of 1000
#    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
#    run.prob  - percentiles of run timing wanted 
#    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data
#             and OpenBUGS will run and stop waiting for your to exit and complete

# force input parameters to be vectors
time     <- as.vector(time)
n1       <- as.vector(n1)
m2       <- as.vector(m2)
u2       <- as.vector(u2)
sampfrac <- as.vector(sampfrac)


#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time all match
if(var(c(length(n1),length(m2),length(u2),length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of n1, m2, u2, sampfrac, time must all be equal. They are:",
        length(n1),length(m2),length(u2),length(sampfrac),length(time),"\n")
   return()}
if(length(logitP.cov) %% length(n1) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of n1 etc They are:",
        length(n1),length(logitP.cov),dim(logitP.cov),"\n")
   return()}
#  2. Check that m2<= n1
if(any(m2>n1, na.rm=TRUE)){
   cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:",n1,"\n m2:",m2,"\n")
   return()}
#  3. Elements of bad.n1, bad.m2, bad.u2, and jump.after must belong to time
if(!all(bad.n1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.n1 must be elements of strata identifiers. You entered \n bad.n1:",bad.n1,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.m2 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",bad.m2,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2 must be elements of strata identifiers. You entered \n bad.u2:",bad.u2,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(jump.after %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** jump.after must be elements of strata identifiers. You entered \n jump.after:",jump.after,"\n Strata identifiers are \n time:",time, "\n")
   return()}

#  4. check that index of logitP.fixed belong to time
if(!all(logitP.fixed %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** logitP.fixed must be elements of strata identifiers. You entered \n logitP.fixed:",logitP.fixed,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(length(logitP.fixed)!=length(logitP.fixed.values)){
   cat("***** ERROR ***** Lengths of logitP.fixed and logitP.fixed.values must all be equal. They are:",
        length(logitP.fixed),length(logitP.fixed.values),"\n")
   return()}

# Deprecation of sampling fraction.
if(any(sampfrac != 1)){
   cat("***** ERROR ***** Sampling fraction is deprecated for any values other than 1. DO NOT USE ANYMORE. ")
   return()
}


results.filename <- paste(prefix,"-results.txt",sep="")   

# Create the report
stdout <- vector('character')
report <- textConnection('stdout', 'wr', local = TRUE)
 
sink(report)
#sink(results.filename)
cat(paste("Time Stratified Petersen with Diagonal recaptures and error in smoothed U - ", date()))
cat("\nVersion:", version, "\n\n")

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** \n")
temp<- cbind(time, n1, m2, u2, round(sampfrac,digits=2), logitP.cov)
colnames(temp)<- c('time', 'n1','m2','u2','SampFrac', paste("logitPcov[", 1:ncol(as.matrix(logitP.cov)),"]",sep="") )
print(temp) 
cat("\n\n")
cat("Jump point are after strata: ", jump.after)
if(length(jump.after)==0) cat("none - A single spline is fit")

cat("\nFixed logitP indices are:", logitP.fixed)
if(length(logitP.fixed)==0) cat("none - NO fixed values")
cat("\nFixed logitP values  are:", logitP.fixed.values)
if(length(logitP.fixed)==0) cat("none - NO fixed values")

cat("\n\nValues of bad.n1 are : ", bad.n1, ". The value of n1 will be set to 1 and m2 to NA for these strata")
if(length(bad.n1)==0) cat("none.")
cat(  "\nValues of bad.m2 are : ", bad.m2, ". The value of m2 will be set to NA for these strata")
if(length(bad.m2)==0) cat("none.")
cat(  "\nValues of bad.u2 are : ", bad.u2, ". The value of u2 will be set to NA for these strata")
if(length(bad.u2)==0) cat("none.")

# Pooled Petersen estimator over ALL of the data including when no releases take place, bad.n1, bad.m2, bad.u2 and missing values.
cat("\n\n*** Pooled Petersen Estimate based on pooling over ALL strata")
cat("\nValues of u2 are adjusting for sampling fraction \n\n")
cat("Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",sum(u2/sampfrac, na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2/sampfrac, na.rm=TRUE))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")


# Obtain the Pooled Petersen estimator after EXCLUDING strata with missing data or strata that are flagged as having bad.n1, bad.m2, and bad.u2  values
temp.n1 <- n1
temp.n1[match(bad.n1,time)] <- NA
temp.m2 <- m2
temp.m2[match(bad.m2,time)] <- NA
temp.u2 <- u2
temp.u2[match(bad.u2,time)] <- NA

select <- (n1>0) & (!is.na(temp.n1)) & (!is.na(temp.m2)) & (!is.na(temp.u2))
cat("\n\n*** Pooled Petersen Estimate after EXCLUDING strata with missing value or flagged as bad.n1, bad.m2 or bad.m2. ")
cat("\nValues of u2 are adjusted for sampling fraction\n\n")
cat("The following strata are excluded because n1=0, NA values in m2 or u2, or flagged by bad.n1, bad.m2 or bad.u2:", time[!select],"\n\n")

temp.n1 <-       n1[select]
temp.m2 <-       m2[select]
temp.u2 <-       u2[select]
temp.sampfrac <- sampfrac[select]

cat("Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2/temp.sampfrac),"\n\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2/temp.sampfrac))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

   
# Obtain Petersen estimator for each stratum prior to excluding any strata flagged as bad values
cat(  "*** Stratified Petersen Estimator for each stratum PRIOR to removing strata with bad.n1, bad.m2, or bad.u2 values.")
cat("\n    Values of u2 are adjusted for sampling fraction ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")


# Obtain Petersen estimator for each stratum after excluding strata where n1=0, or flagged by bad.n1,  bad.m2, or bad.u2
cat(  "*** Stratified Petersen Estimator for each stratum EXCLUDING strata with n1=0, NA values, or flagged by bad.n1, bad.m2, or bad.u2 values ***")
cat("\n    Values of u2 are adjusted for sampling fraction ***\n\n")
temp.n1 <- n1
temp.n1[match(bad.n1,time)] <- NA  # if any value is bad, exclude this entire stratum
temp.n1[match(bad.m2,time)] <- NA
temp.n1[match(bad.u2,time)] <- NA
temp.n1[temp.n1==0]         <- NA  # if n1 is zero, then there is no estimate of capture probability for this stratum
temp.m2 <- m2
temp.m2[match(bad.n1,time)] <- NA  # if any value is bad, exlude this entire stratum
temp.m2[match(bad.m2,time)] <- NA
temp.m2[match(bad.u2,time)] <- NA
temp.u2 <- u2
temp.u2[match(bad.n1,time)] <- NA  # if any value is bad, exclude this entire stratum
temp.u2[match(bad.m2,time)] <- NA
temp.u2[match(bad.u2,time)] <- NA
temp.u2 <- temp.u2/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("Est U(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

# Test if pooling can be done
cat("*** Test if pooled Petersen is allowable on strata without problems in n1 or m2. [Check if marked fractions are equal] ***\n\n")
select <- (temp.n1>0) & (!is.na(temp.n1)) & (!is.na(temp.m2)) 
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
new.u2         <- rep(0, max(time)-min(time)+1)
new.sampfrac   <- rep(0, max(time)-min(time)+1)
new.logitP.cov <- matrix(NA, nrow=max(time)-min(time)+1, ncol=ncol(as.matrix(logitP.cov)))
new.time       <- min(time):max(time)

new.n1[time-min(new.time)+1]         <- n1
new.m2[time-min(new.time)+1]         <- m2
new.m2[match(bad.m2,new.time)]       <- NA  # wipe out where m2 is flagged as bad
new.u2[time-min(new.time)+1]         <- u2
new.u2[match(bad.u2,new.time)]       <- NA  # wipe out where u2 is flagged as bad
new.sampfrac[time-min(new.time)+1]   <- sampfrac
new.logitP.cov[time-min(new.time)+1,]<- as.matrix(logitP.cov)
new.logitP.cov[ is.na(new.logitP.cov[,1]), 1] <- 1  # insert a 1 into first columns where not specified
new.logitP.cov[ is.na(new.logitP.cov)] <- 0         # other covariates are forced to zero not in column 1


# Check for and fix problems with the data
# If n1=m2=0, then set n1 to 1, and set m2<-NA
new.m2[new.n1==0] <- NA
new.n1[new.n1==0] <- 1

new.n1[match(bad.n1,new.time)] <- 1
new.m2[match(bad.n1,new.time)] <- NA 

# Adjust data when a stratum has less than 100% sampling fraction to "estimate" the number
# of unmarked fish that were captured. It is not necessary to adjust the n1 and m2 values 
# as these are used ONLY to estimate the capture efficiency. 
# In reality, there should be a slight adjustment
# to the precision to account for this change, but this is not done.
# Similarly, if the sampling fraction is more than 1, the adjustment is made back to a standard week.
new.u2 <- round(new.u2/new.sampfrac)

# Adjust for strata where sampling fraction=0. On these strata
# u2 is set to NA so that there is NO information on U2 for this stratum
new.u2[new.sampfrac<.001] <- NA

# Print out the revised data
jump.indicator <- rep('   ', max(time)-min(time)+1)
jump.indicator[jump.after-min(time)+1]<- '***'

cat("\n\n*** Revised data *** \n")
temp<- data.frame(time=new.time, n1=new.n1, m2=new.m2, u2=new.u2, 
       sampfrac=round(new.sampfrac,digits=2), new.logitP.cov=new.logitP.cov, 
       jump.indicator=jump.indicator)
print(temp) 
cat("\n\n")


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
cat("   Parameters for prior on beta.logitP[1] (intercept) (mean, 1/var):", round(mu_xiP,3), round(tau_xiP,5),
    " which corresponds to a median P of ", round(expit(mu_xiP),3), "\n")
cat("   Parameters for prior on tauP (residual variance of logit(P) after adjusting for covariates: ",tauP.alpha, tauP.beta, 
    " which corresponds to a mean/std dev of 1/var of:",
    round(tauP.alpha/tauP.beta,2),round(sqrt(tauP.alpha/tauP.beta^2),2),"\n")

cat("\n\n*** Initial seed for this run is: ",InitialSeed, "\n")

sink()

if (debug2) {
   cat("\nprior to formal call to TimeStratPetersenDiagError\n")
   browser()
}

if (debug) 
   {results <- TimeStratPetersenDiagError(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2=new.u2,
            jump.after=jump.after-min(time)+1,
            logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
            n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,  # set to low values for debugging purposes only
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, debug2=debug2, InitialSeed=InitialSeed, save.output.to.files=save.output.to.files)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenDiagError(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2=new.u2, 
            jump.after=jump.after-min(time)+1, 
            logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, debug2=debug2, InitialSeed=InitialSeed, save.output.to.files=save.output.to.files)
   }

# Now to create the various summary tables of the results

# A plot of the observered log(U) on the log scale, and the final mean log(U)
# Create the data frame needed for ggplot. 
# In the diagonal case, time, n1, m2, u2 are the same length

  plot.df   <- data.frame(time =new.time)
  plot.df$logUi <- log((new.u2+1)*(new.n1+2)/(new.m2+1)) # initial guess for logU

# extract the fitted U values
  results.row.names <- rownames(results$summary)
  etaU.row.index    <- grep("etaU", results.row.names)
  etaU<- results$summary[etaU.row.index,]
  plot.df$logU =etaU[,"mean"]
  plot.df$lcl =etaU[,"2.5%"]
  plot.df$ucl =etaU[,"97.5%"]

# extract the spline values
   logUne.row.index <- grep("logUne", results.row.names)
   logUne<- results$summary[logUne.row.index,"mean"]
   plot.df$spline <- results$summary[logUne.row.index,"mean"]

   #browser()
fit.plot <- ggplot(data=plot.df, aes_(x=~time))+
   ggtitle(title, subtitle="Fitted spline curve with 95% credible intervals for estimated log(U[i])")+
   geom_point(aes_(y=~logUi), color="red", shape=1)+  # open circle
   xlab("Time Index\nOpen/closed circles - initial and final estimates")+ylab("log(U[i]) + 95% credible interval")+
   geom_point(aes_(y=~logU), color="black", shape=19)+
   geom_line (aes_(y=~logU), color="black")+
   geom_errorbar(aes_(ymin=~lcl, ymax=~ucl), width=.1)+
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

# plot logit(P) over time
logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=new.m2, u2=new.u2, logitP.cov=new.logitP.cov, results=results)
if(save.output.to.files)ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
results$plots$logitP.plot <- logitP.plot

# Look at autocorrelation function for Utot
mcmc.sample <- data.frame(parm="Utot", sample=results$sims.matrix[,"Utot"], stringsAsFactors=FALSE)
acf.Utot.plot <- plot_acf(mcmc.sample)
if(save.output.to.files)ggsave(plot=acf.Utot.plot, filename=paste(prefix,"-Utot-acf.pdf",sep=""), height=4, width=6, units="in")
results$plots$acf.Utot.plot <- acf.Utot.plot


# Look at the shape of the posterior distribution
mcmc.sample1 <- data.frame(parm="Utot", sample=results$sims.matrix[,"Utot"], stringsAsFactors=FALSE)
mcmc.sample2 <- data.frame(parm="Ntot", sample=results$sims.matrix[,"Ntot"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2)
post.UNtot.plot <- plot_posterior(mcmc.sample)
post.UNtot.plot
if(save.output.to.files)ggsave(plot=post.UNtot.plot, filename=paste(prefix,"-UNtot-posterior.pdf",sep=""),
                               height=ifelse(length(unique(mcmc.sample$parm))==1,4,6), width=6, units="in")
results$plots$post.UNtot.plot <- post.UNtot.plot


#save the Bayesian predictive distribution (Bayesian p-value plots)
#browser()
discrep <-PredictivePosterior.TSPDE (new.n1, new.m2, new.u2,
             new.logitP.fixed,
				     expit(results$sims.list$logitP),
				     round(results$sims.list$U))
gof <- PredictivePosteriorPlot.TSPDE (discrep)  # get the bayesian p-values
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



#sink(results.filename, append=TRUE)
sink(report, append=TRUE)
# Global summary of results
cat("\n\n*** Summary of MCMC results *** \n\n")
  save.max.print <- getOption("max.print")
  options(max.print=.Machine$integer.max)
  
  print(results, digits.summary=3)#, max=.Machine$integer.max)
  
  options(max.print=save.max.print)

# Give an alternate computation of DIC based on the variance of the deviance
# Refer to http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/DIC-slides.pdf for derivation and why
# this alternate method may be superior to that automatically computed by OpenBugs

cat("\n\n*** Alternate DIC computation based on p_D = var(deviance)/2 \n")
results.row.names <- rownames(results$summary)
deviance.row.index<- grep("deviance", results.row.names)
deviance          <- results$summary[deviance.row.index,]
p.D <- deviance["sd"]^2/2
dic <- deviance["mean"]+p.D
cat("    D-bar: ", deviance["mean"],";  var(dev): ", deviance["sd"]^2,
    "; p.D: ", p.D, "; DIC: ", dic)

# Summary of population sizes. Extra code rounds results to integers except for Rhat
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
results$data <- list( time=time, n1=n1, m2=m2, u2=u2, sampfrac=sampfrac, 
                      jump.after=jump.after, 
                      bad.n1=bad.n1, bad.m2=bad.m2, bad.u2=bad.u2, 
                      logitP.cov=logitP.cov, version=version, date_run=date(),
                      title=title)

return(results)
} # end of function
