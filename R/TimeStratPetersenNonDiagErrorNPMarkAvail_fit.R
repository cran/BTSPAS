## 2020-12-15 CJS removed sampfrac from body of code
## 2020-12-15 CJS Fixed problem when specifyin u2==NA
## 2020-11-07 CJS Allowed user to specify prior for beta coefficient for logitP
## 2018-12-22 CJS add code to estimate mean movement vector (movep)
## 2018-12-19 CJS sampling fraction deprecated
## 2018-12-14 CJS bayesian p-value plots added
## 2018-12-06 CJS convert report to textConnections
## 2018-12-02 CJS convert trace plots to ggplot
## 2018-12-01 CJS converted acf, posterior plots to ggplot
## 2018-11-30 CJS Fixed problem of epsilon not being right length
## 2018-11-29 CJS Fixed problem of printing large results
## 2018-11-28 CJS remove reference to OpenBugs
## 2014-09-01 CJS conversion to jags
## 2012-08-30 CJS fixed problem in any() and all() in error checking with NAs
## 2011-02-21 CJS changed u2 to new.u2 in code for expanded.m2
## 2011-02-19 CJS First development

#' Wrapper (*_fit)  to call the function to fit a Time Stratified Petersen Estimator
#' with NON Diagonal Entries with an non-parametric travel time and fall back
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
#' @param u2 A numeric vector of the number of unmarked fish captured in each
#' stratum.  These will be expanded by the capture efficiency to estimate the
#' population size in each stratum. The length of u2 should be between the length of n1 and length n1 + number of columns in m2 -1
#' @template sampfrac
#' @param jump.after A numeric vector with elements belonging to \code{time}.
#' In some cases, the spline fitting the population numbers should be allowed
#' to jump.  For example, the population size could take a jump when hatchery
#' released fish suddenly arrive at the trap.  The jumps occur AFTER the strata
#' listed in this argument.
#' @template bad.n1
#' @template bad.m2
#' @template bad.u2
#' @template logitP.cov
#' @param logitP.fixed A numeric vector (could be null) of the time strata
#' where the logit(P) would be fixed. Typically, this is used when the capture
#' rates for some strata are 0 and logit(P) is set to -10 for these strata. The
#' fixed values are given in \code{logitP.fixed.values}
#' @param logitP.fixed.values A numerical vector (could be null) of the fixed
#' values for logit(P) at strata given by logitP.fixed. Typically this is used
#' when certain strata have a 0 capture rate and the fixed value is set to -10
#' which on the logit scale gives p[i] essentially 0. Don't specify values such
#' as -50 because numerical problems could occur in JAGS.
#' @param marked_available_n Information, usually from prior studies, on the
#' fraction of marks that will be available. The *_n and *_x are used to create
#' a "binomial" distribution for information on the marked availability. For
#' example, if *_n=66 and *_x=40, then you estimate that about 40/66=61\% of marks
#' are available and 39\% have dropped out or fallen back.
#' @param marked_available_x See marked_available_n
#' @template mcmc-parms
#' @param tauU.alpha One of the parameters along with \code{tauU.beta} for the
#' prior for the variance of the random noise for the smoothing spline.
#' @param tauU.beta One of the parameters along with \code{tauU.alpha} for the
#' prior for the variance of the random noise for the smoothing spline.
#' @param taueU.alpha One of the parameters along with \code{taueU.beta} for
#' the prior for the variance of noise around the spline.
#' @param taueU.beta One of the parameters along with \code{taueU.alpha} for
#' the prior for the variance of noise around the spline.
#' @param Delta.max Maximum transition time for marked fish, i.e. all fish
#' assumed to have moved by Delta.max unit of time
#' @param tauTT.alpha One of the parameters along with \code{tauTT.beta} for
#' the prior on 1/var of logit continuation ratio for travel times
#' @param tauTT.beta One of the parameters along with \code{tauTT.alpha} for
#' the prior on 1/var of logit continuation ratio for travel times
#' @param prior.beta.logitP.mean Mean of the prior normal distribution for
#' logit(catchability) across strata
#' @param prior.beta.logitP.sd   SD of the prior normal distribution for
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
#' @template InitialSeed
#' @template save.output.to.files
#' 
#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' @template author 
#' @template references
#' @keywords ~models ~smooth
#' @examples
#'  
#' ##---- See the vignettes  for examples of how to use this package
#' 
#' @export TimeStratPetersenNonDiagError_fit
#' @importFrom stats runif var sd


#' @export TimeStratPetersenNonDiagErrorNPMarkAvail_fit

TimeStratPetersenNonDiagErrorNPMarkAvail_fit<- function( title="TSPNDENP-avail", prefix="TSPNDENP-avail-", 
                         time, n1, m2, u2, sampfrac=rep(1,length(u2)), jump.after=NULL,
                         bad.n1=c(), bad.m2=c(), bad.u2=c(),
                         logitP.cov=rep(1,length(u2)),
                         logitP.fixed=NULL, logitP.fixed.values=NULL, 
                         marked_available_n, marked_available_x,
                         n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
                         tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
                         prior.beta.logitP.mean = c(logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),rep(0,  ncol(as.matrix(logitP.cov))-1)),
                         prior.beta.logitP.sd   = c(2,                                           rep(10, ncol(as.matrix(logitP.cov))-1)), 
                         tauP.alpha=.001, tauP.beta=.001,
                         Delta.max=NULL,tauTT.alpha=.1,tauTT.beta=.1,
                         run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
                         debug=FALSE, debug2=FALSE,
                         InitialSeed=ceiling(stats::runif(1,min=0, max=1000000)),
                         save.output.to.files=TRUE) {
  ## Fit a Time Stratified Petersen model with NON-diagonal entries and with smoothing on U allowing for random error
  ## and fall back after tagging. This is based on the Skeena River study, where only 40/66 (60%) acoustically tagged fish
  ## were observed above the canyon spot and hence 40% of tagged fish never migrated forward of their tagging release spot.
  ## This reduces the number of tagged fish available for recapture and so, if not accounted for, leads to 
  ## positive biases in the estimates of abundance.

  ## This is the classical stratified Petersen model where the recoveries can take place for this and multiple
  ## strata later. Transisions of marked fish are modelled non-parametrically.
  ##
  
  version <- '2021-11-01'
  options(width=200)
  
  ## Input parameters are
  ##    title  - title for the analysis (character string)
  ##    prefix - prefix used for files created with the analysis results
  ##             this should be in standard Window's format, eg. JC-2002-ST-TSPNDE
  ##             to which is appended various suffixes for plots etc (character string)
  ##    time   - vector of stratum numbers. For example, 9:38 would indicate that the
  ##             Trinity River system sampled weeks 9 to 38. 
  ##             These values are never used in the analysis and only serve as labels for the weeks and for plotting purposes.
  ##             They should be contiguous equally spaced values and be the same length as u2.
  ##    n1, m2, u2 - the input data consisting of fish marked and released, recapture, and unmarked captured
  ##             Note that m2 is a MATRIX. The first column are those fish recaptured in the stratum of release
  ##             and the subsequent columns are those recoveries in later strata.
  ##             This is expanded to the full matrix [i,j] for released released in stratum i and recovered in stratum j
  ##             The vector u2 should be long enough to account for any fish that are recaptured later on
  ##             from releases late in the season. The bottom right diagonal of m2 may be all zeros - that is ok
  ##             Notice that length(u2) can be longer than length(n1)+nrow(m2).
  ##    sampfrac - Deprecated. DO NOT USE.
  ##    jump.after - in some cases, a single spline is still not flexible enough to cope with rapid
  ##                 changes in the run curve. For example, in the Trinity River project, a larger
  ##                 hatchery release occurs around stratum 14. This is a vector indicating the
  ##                 strata AFTER which the spline curve is allowed to jump.
  ##                 null or vector of arbitrary length.
  ##    bad.n1  - list of stratum numbers where the value of n1 is suspect.
  ##    bad.m2  - list of stratum numbers where the value of m2 is suspect.
  ##              For example, the capture rate could be extremely low.
  ##              These are set to NA prior to the call to JAGS
  ##    bad.u2  - list of stratum numbers where the value of u2 is suspect.
  ##    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
  ##               for the first element of the covariance matrix and 0 for the rest of the covariates.
  ##               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
  ##               if you have any covariates. The default, if not specified, is a constant (the mean logit)
  ##    marked_available_n, marked_available_x  Information on the movement forward rate. Treat this a binomial data
  ##       which will be applied uniformly over all released. For example, use *_n=60 and *_x=40 to represent
  ##       data from the telemetry study that had 40/60 tagged fish move forward. You can vary the precision of the
  ##       estimate of the marked_availability_fraction by chaning the _n and _x values to create the illusion
  ##       of better and worse information on the availability value.
  ##
  ##    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
  ##    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline 
  #     prior.beta.logitP.mean, prior.beta.logitP.sd   - parameters for the prior on mean logit(P)'s [The intercept term]
  #                              The other covariates are assigned priors of a mean of 0 and a sd of 30
  ##    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
  ##    Delta.max - maximum transition time for marked fish
  ##    tauTT.alpha, tauTT.beta - parameters of the prior on 1/var of logit continuation ratio for travel times
  ##    run.prob  - percentiles of run timing wanted 
  ##    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data

# force input vectors to be vectors. Note that m2 is NOT a vector
time     <- as.vector(time)
n1       <- as.vector(n1)
u2       <- as.vector(u2)
sampfrac <- as.vector(sampfrac)
  
  ##  Do some basic error checking
  ##  1. Check that length of n1, m2, u2, sampfrac, time are consistent with each other.
  ##  In the non-diagonal case, they don't have to match
if(length(n1)!=nrow(m2))
    stop("***** ERROR ***** Length of n1 and number of rows of m2 must be equal. They are:",
        length(n1)," ",nrow(u2),"\n")
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

  if(stats::var(c(length(u2),length(sampfrac),length(time)))>0)
    stop("***** ERROR ***** Lengths of u2, sampfrac, time must all be equal. They are:",
         length(u2),' ',length(sampfrac),' ',length(time),"\n")

    if(length(logitP.cov) %% length(u2) != 0)
      stop("***** ERROR ***** Dimension of covariate vector doesn't match length of u2. They are:",
        length(u2),' ',length(logitP.cov),' ',dim(logitP.cov),"\n")

  ##  2. Check that rowsum of m2<= n1
  if(any(apply(m2,1,sum, na.rm=TRUE)>n1))
    stop("***** ERROR ***** m2[i,+] must be <= n1[i]. The arguments are \n n1:",paste(n1,collapse=","),
         "\n m2:",paste(m2,collapse=","),"\n")

  ##  3. Elements of bad.m2 and jump.after must belong to time
  if(!all(bad.n1 %in% time, na.rm=TRUE))
    stop("***** ERROR ***** bad.n1 must be elements of strata identifiers. You entered \n bad.n1:",
         paste(bad.n1,collapse=","),"\n Strata identifiers are \n time:",
         paste(time,collapse=","),"\n")

  if(!all(bad.m2 %in% time, na.rm=TRUE))
    stop("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",
         paste(bad.m2,collapse=","),"\n Strata identifiers are \n time:",
         paste(time,collapse=","),"\n")

  if(!all(bad.u2 %in% time, na.rm=TRUE))
    stop("***** ERROR ***** bad.u2 must be elements of strata identifiers. You entered \n bad.u2:",
         paste(bad.u2,collapse=","),"\n Strata identifiers are \n time:",
         paste(time,collapse=","), "\n")

  if(!all(jump.after %in% time, na.rm=TRUE))
    stop("***** ERROR ***** jump.after must be elements of strata identifiers. You entered \n jump.after:",
         paste(jump.after,collapse=","),"\n Strata identifiers are \n time:",
         paste(time,collapse=","), "\n")

  #  4. check that index of logitP.fixed belong to time
  if(!all(logitP.fixed %in% time, na.rm=TRUE)){
    cat("***** ERROR ***** logitP.fixed must be elements of strata identifiers. You entered \n logitP.fixed:",
        paste(logitP.fixed,collapse=","),"\n Strata identifiers are \n time:",
        paste(time,collapse=","), "\n")
    return()}
  if(length(logitP.fixed)!=length(logitP.fixed.values)){
    cat("***** ERROR ***** Lengths of logitP.fixed and logitP.fixed.values must all be equal. They are:",
        length(logitP.fixed),length(logitP.fixed.values),"\n")
    return()}
 
  # 5. Check that some basic information on marked availability is given
  if( is.na(marked_available_n) | is.na(marked_available_x) | marked_available_x > marked_available_n){
    cat("***** ERROR ***** Bad marked_availability values. You entered:",marked_available_n," ",marked_available_x,"\n")
    return()}

  #7 check that the length of u2 
  if(length(u2) < length(n1) | length(u2) > (length(n1)+ ncol(m2)-1)){
    cat("***** ERROR ***** Length(u2) must between length(n1) and length(n1)+ncol(m2) \n")
    return()
  }

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

  ## Define maximum travel time if not supplied by user
  if(is.null(Delta.max)) Delta.max <- ncol(m2)-1
 
  ## Define output filename
  results.filename <- paste(prefix,"-results.txt",sep="")   

  ## Open sink to output file
  stdout <- vector('character')
  report <- textConnection('stdout', 'wr', local = TRUE)
  sink(report)

  cat(paste("Time Stratified Petersen with Non-Diagonal recaptures, error in smoothed U, non-parametric modelling of travel times, and incorporating mark availability- ", date()))
  cat("\nVersion: ", version)
  
  cat("\n\n", title, "Results \n\n")
  
  ## m2(i,+) < n1(i)

  cat("*** Raw data *** (padded to match length of u2) \n")
  jump.indicator <- rep('   ', length(u2))
  jump.indicator[time %in% jump.after]<- '***'
  ex.n1 <- c(n1, rep(NA, length(u2)-length(n1)))
  ex.m2 <- rbind(m2,matrix(NA, nrow=length(u2)-length(n1), ncol=ncol(m2))) 
  temp<- data.frame(time=time, n1=ex.n1, m2=ex.m2, u2=u2, logitP.cov=logitP.cov, jump=jump.indicator)
  print(temp) 
  cat("\n\n")

  cat("*** Marked Availability prior information *** \n")
  cat("    Set marked available n=", marked_available_n," with x=",marked_available_x,"\n\n\n")

  ## Print information about jump points
  cat("Jump point are after strata: ", jump.after)
  if(length(jump.after)==0) cat("none - A single spline is fit")

  ## Print information about delta max
  cat("\nMaximum travel time (Delta.max): ",Delta.max)
  cat("\nFixed logitP indices are: ", logitP.fixed)
  if(length(logitP.fixed)==0) cat("none - NO fixed values")
  cat("\nFixed logitP values  are: ", logitP.fixed.values)
  if(length(logitP.fixed)==0) cat("none - NO fixed values")
  
  ## Obtain the Pooled Petersen estimator prior to fixup of bad.n1, bad.m2, and bad.u2 values
  cat("\n\n*** Pooled Petersen Estimate prior to fixing bad n1, m2, or u2 values  CHECK - CHECK - CHECK - CHECK ***\n\n")
  cat("    *** NOT ADJUSTED FOR MARK AVAILABILITY/Dropout/Fallback ***\n")
  temp.n1 <- n1
  temp.m2 <- m2
  temp.u2 <- u2
  
  cat("Total n1=", sum(temp.n1,na.rm=TRUE),";  m2=",sum(temp.m2,na.rm=TRUE),";  u2=",sum(temp.u2,na.rm=TRUE),"\n\n")
  pp <- SimplePetersen(sum(temp.n1,na.rm=TRUE), sum(temp.m2,na.rm=TRUE), sum(temp.u2,na.rm=TRUE))
  cat("Est U(total) not adjusted for fallback ", format(round(pp$U.est),big.mark=","),
      "  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
  cat("Est N(total) not adjusted for fallback ", format(round(pp$N.est),big.mark=","),
      "  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")
  
  # adjustment for dropout
  dr <- 1-marked_available_x/marked_available_n # dropout probability
  se_dr <- sqrt(dr*(1-dr)/marked_available_n)
  cat("\n\nAdjusting for fallback/dropout  \n")
  cat("Estimated dropout is", dr, "with se of ", se_dr, "\n")

  # adjust the petersen estimator for drop out including the uncertainty in the dropout probability
  pp.adj <- pp
  pp.adj$N.est <- pp.adj$N.est * (1-dr)
  pp.adj$N.se  <- sqrt(pp$N.se^2 * se_dr^2+
                     pp$N.se^2 * (1-dr)^2 +
                     pp$N.est^2 * se_dr^2)
  pp.adj$U.est <- pp.adj$U.est * (1-dr)
  pp.adj$U.se  <- sqrt(pp$U.se^2 * se_dr^2+
                     pp$U.se^2 * (1-dr)^2 +
                     pp$U.est^2 * se_dr^2)
  
  cat("Est U(total) adjusting for dropout is ", format(round(pp.adj$U.est),big.mark=","),
      "  (SE ", format(round(pp.adj$U.se), big.mark=","), ")\n")
  cat("Est N(total) adjusting for dropout is ", format(round(pp.adj$N.est),big.mark=","),
      "  (SE ", format(round(pp.adj$N.se), big.mark=","), ")\n\n\n")
  
  
  ## Test if pooling can be done
  cat("*** Test if pooled Petersen is allowable. [Check if fraction captured equal] ***\n\n")
  select <- (n1>0) & (!is.na(n1)) & (!is.na(apply(m2,1,sum)))
  temp.n1 <- n1[select]
  temp.m2 <- m2[select,]
  test <- TestIfPool( temp.n1, apply(temp.m2,1,sum))
  cat("(Large Sample) Chi-square test statistic ", test$chi$statistic," has p-value", test$chi$p.value,"\n\n")
  temp <- cbind(time[1:length(n1)][select],test$chi$observed, round(test$chi$expected,1), round(test$chi$residuals^2,1))
  colnames(temp) <- c('time','n1-m2*','m2*','E[n1-m2]','E[m2]','X2[n1-m2]','X2[m2]')
  print(temp)
  cat("\n Be cautious of using this test in cases of small expected values. \n\n")

  
  ## Adjust the data for the explicity bad values or other problems
  new.time <- time
  new.n1   <- n1
  new.m2   <- m2
  new.u2   <- u2
  new.logitP.cov <- logitP.cov
  

  ## Set the bad values to missing 
  new.n1[time[1:length(n1)] %in% c(bad.n1, bad.m2)] <- 0
  new.m2[time[1:length(n1)] %in% c(bad.m2, bad.n1)] <- 0
  new.u2[time %in% bad.u2]                          <- NA
  
  ## Print out the revised data
  cat("\n\n*** Revised data *** \n")
  jump.indicator <- rep('   ', length(u2))
  jump.indicator[time %in% jump.after]<- '***'
  ex.n1 <- c(new.n1, rep(NA, length(new.u2)-length(new.n1)))
  ex.m2 <- rbind(new.m2,matrix(NA, nrow=length(new.u2)-length(new.n1), ncol=ncol(new.m2))) 
  temp<- data.frame(time=new.time, n1=ex.n1, m2=ex.m2, u2=new.u2, logitP.cov=new.logitP.cov,
                    jump.after=jump.indicator)
  print(temp) 
  cat("\n\n")

  cat("*** Marked Availability prior information *** \n")
  cat("    Set marked available n=", marked_available_n," with x=",marked_available_x,"\n\n\n")
  
  ## The NP analysis does not need the expanded m2 array, but this is
  ## needed late on. So, we'd better compute it here. The last column
  ## of this matrix will be the number of individuals from each
  ## stratum that are not recaptured.
  ##
  
  expanded.m2 <- matrix(0, nrow=length(new.n1), ncol=length(new.u2)+1)
  for(i in 1:length(new.n1)){
    expanded.m2[i,1:length(new.u2)] <- c(rep(0,i-1),new.m2[i,],rep(0,length(new.u2)))[1:length(new.u2)]
    expanded.m2[i,length(new.u2)+1] <- new.n1[i] - sum(new.m2[i,])
  }
  
  cat("*** Expanded m2 array with column sum and u2 ***\n\n")
  save.max.print <- getOption("max.print")
  options(max.print=.Machine$integer.max)
 
  temp <- rbind(expanded.m2, apply(expanded.m2,2,sum, na.rm=TRUE))
  rownames(temp)[nrow(temp)] <- 'Column totals'
  temp <- rbind(temp, c(u2, rep(NA, ncol(expanded.m2)-length(u2)) ))
  rownames(temp)[nrow(temp)] <- "Untagged (u2)"
  temp <- rbind(temp, c(new.u2, rep(NA, ncol(expanded.m2)-length(new.u2)) ))
  rownames(temp)[nrow(temp)] <- "Untagged - after fixups"
  
  new.logitP.fixed <- rep(NA, length(new.u2))
  new.logitP.fixed[match(logitP.fixed, time)] <- logitP.fixed.values
  
  temp <- rbind(temp, c(new.logitP.fixed, rep(NA, ncol(expanded.m2)-length(new.u2)) ))
  rownames(temp)[nrow(temp)] <- "Logit P fixed"
  rownames(temp)[1:length(n1)] <- 1:length(n1)
  print(temp)
  options(max.print=save.max.print)
  sink()
  
  # some further checking on u2. Make sure that every columns where there are recoveries has a u2
  # browser()
  if( (length(u2)+1) <= (ncol(temp)-1)) {
     if(any( temp["Column totals", (length(u2)+1):(ncol(temp)-1)] >0)){
       cat("***** ERROR ***** Non-zero recoveries and u2 not available at end of experiment??? \n Check above matrix\n")
       return()
     }
  }

  sink(report, append=TRUE)
  # assign the logitP fixed values etc.
  new.logitP.fixed <- rep(NA, length(new.u2))
  new.logitP.fixed[match(logitP.fixed, time)] <- logitP.fixed.values

  ## We do need to add the column of not recaptured counts to the m2
  ## array.
  new.m2 <- cbind(new.m2,new.n1-apply(new.m2,1,sum))

  ## We construct a prior probability on the P(marks available) based on the information provided
  ## by assuming a beta prior that would give the binomial results
  ma.p.alpha <- marked_available_x
  ma.p.beta  <- marked_available_n - marked_available_x
  
  ## Print out information on the prior distributions used
  cat("\n\n*** Information on priors *** \n")
 
  ## 0) ma.p = p(marked_availability for subsequent recapture)
  cat("   P(marked fish available for subsequent recapture) has beta(",ma.p.alpha,ma.p.beta,") which corresponds \n",
      "   to a mean of ", round(ma.p.alpha/(ma.p.alpha+ma.p.beta),2),' and sd of ',
          round(sqrt(ma.p.alpha*ma.p.beta/(ma.p.alpha+ma.p.beta+1)/(ma.p.alpha+ma.p.beta)**2),3),"\n")

  ## 1) tauU = (variance of spline coefficients)^-1
  cat("   Parameters for prior on tauU (variance in spline coefficients): ", tauU.alpha, tauU.beta, 
      " which corresponds to a mean/std dev of 1/var of:",
      round(tauU.alpha/tauU.beta,2),round(sqrt(tauU.alpha/tauU.beta^2),2),"\n")

  ## 2) taueU = (variance of errors)^-1
  cat("   Parameters for prior on taueU (variance of log(U) about spline): ",taueU.alpha, taueU.beta, 
      " which corresponds to a mean/std dev of 1/var of:",
      round(taueU.alpha/taueU.beta,2),round(sqrt(taueU.alpha/taueU.beta^2),2),"\n")

  ## 3) prior.beta.logitP = priors for coefficients of covariates for logitP
  cat("   Parameters for prior on beta.logitP[1] (intercept) (mean, sd): \n", cbind(round(prior.beta.logitP.mean,3), round(prior.beta.logitP.sd,5)),"\n")

  ## 4) tauP = (variance of capture probabilites conditional on covariates)^-1
  cat("   Parameters for prior on tauP (residual variance of logit(P) after adjusting for covariates): ",tauP.alpha, tauP.beta, 
      " which corresponds to a mean/std dev of 1/var of:",
      round(tauP.alpha/tauP.beta,2),round(sqrt(tauP.alpha/tauP.beta^2),2),"\n")

  ## 5) tauTT = (variance of continuation ratios for theta)^-1
  cat("   Parameters for prior on tauTT (variance of continuation rations for travel times): ",tauTT.alpha, tauTT.beta, 
      " which corresponds to a mean/std dev of 1/var of:",
      round(tauTT.alpha/tauTT.beta,2),round(sqrt(tauTT.alpha/tauTT.beta^2),2),"\n")

  cat("\n\nInitial seed for this run is: ",InitialSeed, "\n")
  sink()

  if (debug2) {
    cat("\nprior to formal call to TimeStratPetersenNonDiagError\n")
    browser()
  }
 
  if (debug) 
   {results <- TimeStratPetersenNonDiagErrorNPMarkAvail(
         title=title, prefix=prefix, 
         time=new.time, n1=new.n1, m2=new.m2, u2=new.u2,
         jump.after=(1:length(u2))[time %in% jump.after],
         logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
         ma.p.alpha, ma.p.beta,
         n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,  # set to small values for debugging only
         prior.beta.logitP.mean=prior.beta.logitP.mean, 
         prior.beta.logitP.sd  =prior.beta.logitP.sd,
         tauU.alpha=tauU.alpha, tauU.beta=tauU.beta,
         taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
         Delta.max=Delta.max,tauTT.alpha=tauTT.alpha,tauTT.beta=tauTT.beta,
         debug=debug, debug2=debug2, InitialSeed=InitialSeed,
         save.output.to.files=save.output.to.files)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenNonDiagErrorNPMarkAvail(
        title=title, prefix=prefix, 
        time=new.time, n1=new.n1, m2=new.m2, u2=new.u2,
        jump.after=(1:length(u2))[time %in% jump.after],
        logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
        ma.p.alpha, ma.p.beta,
        n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims, 
        prior.beta.logitP.mean=prior.beta.logitP.mean, 
        prior.beta.logitP.sd  =prior.beta.logitP.sd,
        tauU.alpha=tauU.alpha, tauU.beta=tauU.beta,
        taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
        Delta.max=Delta.max,tauTT.alpha=tauTT.alpha,tauTT.beta=tauTT.beta,
        debug=debug, debug2=debug2, InitialSeed=InitialSeed,
        save.output.to.files=save.output.to.files)
   } 
  

  results$PP$using.all.data <-pp 
  results$PP$using.all.data.fallback <- pp.adj
  results$dr <- data.frame(est=dr, se=se_dr)
  
  ## Now to create the various summary tables of the results
  
  Nstrata.rel <- length(n1)
  Nstrata.cap <- ncol(expanded.m2) -1 ## don't forget that last column of m2 is number of fish never seen
  
  # A plot of the observered log(U) on the log scale, and the final mean log(U)
  plot.df   <- data.frame(time =new.time)
  plot.df$logUi <-log( c((new.u2[1:Nstrata.rel]+1)*(new.n1+2)/(apply(expanded.m2[,1:Nstrata.cap],1,sum)+1), rep(NA, length(u2)-Nstrata.rel)))

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

  fit.plot <- ggplot(data=plot.df, aes_(x=~time))+
    ggtitle(title, subtitle="Fitted spline curve with 95% credible intervals for estimated log(U[i])")+
    geom_point(aes_(y=~logUi), color="red", shape=1)+  # open circle
    xlab("Time Index\nOpen/closed circles - initial and final estimates")+
    ylab("log(U[i]) + 95% credible interval")+
    geom_point(aes_(y=~logU), color="black", shape=19)+
    geom_line (aes_(y=~logU), color="black")+
    geom_errorbar(aes_(ymin=~lcl, ymax=~ucl), width=.1)+
    geom_line(aes_(y=~spline),linetype="dashed")+
    scale_x_continuous(breaks=seq(min(plot.df$time,na.rm=TRUE),max(plot.df$time, na.rm=TRUE),2))+
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

  
  ## acf plot
  logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2, logitP.cov=new.logitP.cov, results=results) 
  if(save.output.to.files)ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
  results$plots$logitP.plot <- logitP.plot
  
  ## Look at autocorrelation function for Utot
  mcmc.sample <- data.frame(parm="Utot", sample=results$sims.matrix[,"Utot"], stringsAsFactors=FALSE)
  acf.Utot.plot <- plot_acf(mcmc.sample)
  if(save.output.to.files)ggsave(plot=acf.Utot.plot, filename=paste(prefix,"-Utot-acf.pdf",sep=""), height=4, width=6, units="in")
  results$plots$acf.Utot.plot <- acf.Utot.plot
  
  ## Look at the shape of the posterior distribution  browser()
  mcmc.sample1 <- data.frame(parm="Utot", sample=results$sims.matrix[,"Utot"], stringsAsFactors=FALSE)
  mcmc.sample2 <- data.frame(parm="Ntot", sample=results$sims.matrix[,"Ntot"], stringsAsFactors=FALSE)
  mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2)
  post.UNtot.plot <- plot_posterior(mcmc.sample)
  post.UNtot.plot
  if(save.output.to.files)ggsave(plot=post.UNtot.plot, filename=paste(prefix,"-UNtot-posterior.pdf",sep=""),
                               height=ifelse(length(unique(mcmc.sample$parm))<=2,4,6), width=6, units="in")
  results$plots$post.UNtot.plot <- post.UNtot.plot

  ## Bayesian P-values
  discrep <-PredictivePosterior.TSPNDENPMarkAvail(new.n1, expanded.m2, new.u2,
                                         new.logitP.fixed,
                                         expit(results$sims.list$logitP),
                                         round(results$sims.list$U),
                                         results$sims.list$Theta,
                                         results$sims.list$ma.p,
                                         Delta.max)
  gof <- PredictivePosteriorPlot.TSPNDE (discrep)
  if(save.output.to.files)ggsave(gof[[1]],filename=paste(prefix,"-GOF.pdf",sep=""),  height=8, width=8, units="in", dpi=300 )
  results$plots$gof <- gof

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
  
  ## Global summary of results
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
  
  ## Summary of population sizes
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

  
  ## add some of the raw data to the bugs object for simplicity in referencing it later
  results$data <- list( time=time, n1=n1, m2=m2, u2=u2, 
                       jump.after=jump.after, 
                       bad.n1=bad.n1, bad.m2=bad.m2, bad.u2=bad.u2, 
                       logitP.cov=logitP.cov,
                       version=version, date_run=date(),title=title)

  return(results)
} ## end of function
