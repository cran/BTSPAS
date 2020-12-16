# 2020-12-15 CJS Removed all references to sampfrac in the code
# 2020-11-07 CJS Allowed user to specify prior for beta coefficient for logitP
# 2018-12-19 CJS deprecation of sampling fraction
# 2018-12-06 CJS saved output to a textconnection that is saved
# 2018-12-02 CJS converted trace plots to ggplot
# 2018-12-01 CJS converted posterior plot to ggplot
# 2018-11-30 CJS converted acf plot to ggplot
# 2018-11-29 CJS fixed problem where print got cut off in large problems
# 2018-11-28 CJS removed reference of OpenBugs
# 2015-06-10 CJS converted gof plots to ggplot(). Bug fix.
# 2014-09-01 CJS converstion to JAGS
# 2012-08-30 CJS fixed problem with missing values in any() and all()
# 2011-06-13 CJS added p-values to results
# 2010-11-25 CJS pretty printing of final estimates of population sizes
# 2010-09-06 CJS forced input vectors to be vectors
# 2010-08-06 CJS added creation of traceplots
# 2010-08-03 CJS added version/date to final object
# 2010-03-29 CJS Inital version of code

#' @rdname TimeStratPetersenDiagErrorWHChinook_fit
#' @export TimeStratPetersenDiagErrorWHChinook2_fit

TimeStratPetersenDiagErrorWHChinook2_fit<- 
       function( title="TSPDE-WHChinook2", prefix="TSPDE-WHChinook2-", 
                 time, n1, m2, 
                 u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1,
                 clip.frac.H.YoY, clip.frac.H.1, sampfrac=rep(1,length(u2.A.YoY)), 
                 hatch.after.YoY=NULL, 
                 bad.m2=c(), bad.u2.A.YoY=c(), bad.u2.N.YoY=c(), bad.u2.A.1=c(), bad.u2.N.1=c(),
                 logitP.cov=as.matrix(rep(1,length(n1))),
                 n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
                 tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
                 prior.beta.logitP.mean = c(logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),rep(0,  ncol(as.matrix(logitP.cov))-1)),
                 prior.beta.logitP.sd   = c(sd(logit((m2+.5)/(n1+1)),na.rm=TRUE),        rep(10, ncol(as.matrix(logitP.cov))-1)), 
                 tauP.alpha=.001, tauP.beta=.001,
                 run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
                 debug=FALSE, debug2=FALSE, 
                 InitialSeed=ceiling(runif(1,min=0,1000000)),
                 save.output.to.files=TRUE) {
# Fit a Time Stratified Petersen model with diagonal entries and with smoothing on U allowing for random error,
# covariates for the the capture probabilities, and separating the YoY and Age1 wild vs hatchery fish
# The "diagonal entries" implies that no marked fish are recaptured outside the (time) stratum of release
#
   version <- '2021-01-01'
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
#    u2.A.YoY  - number of YoY unmarked fish with adipose fin clips
#    u2.N.YoY  - number of YoY unmarked fish with NO adipose fin clips
#               All YoY wild fish have NO adipose fin clips; however, hatchery fish are a mixture
#               of fish with adipose fin clips (a known percentage are marked) unmarked fish.
#               So u2.A.YoY MUST be hatchery fish.
#                  u2.N.YoY is a mixture of wild and hatchery fish.
#    u2.A.1  - number of Age1 unmarked fish with adipose fin clips
#    u2.N.1  - number of Age1 unmarked fish with NO adipose fin clips
#               All Age1 wild fish have NO adipose fin clips; however, hatchery fish are a mixture
#               of fish with adipose fin clips (a known percentage are marked) unmarked fish.
#               So u2.A.1 MUST be hatchery fish.
#                  u2.N.1 is a mixture of wild and hatchery fish.
#    clip.frac.H.YoY - what fraction of the YoY hatchery fish are clipped?
#    clip.frac.H.1   - what fraction of the Age1 hatchery fish are clipped (from last year's releases)?
#    sampfrac - Deprecated **** DO NOT USE ANYMORE **** sampling fraction to adjust for how many days of the week was the trap operating
#    hatch.after - julian week AFTER which hatchery fish are released 
#    bad.m2  - list of julian numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to JAGS
#    bad.u2.A.YoY - list of julian weeks where the value of u2.A.YoY is suspect. 
#               These are set to NA prior to the call to JAGS
#    bad.u2.N.YoY - list of julian weeks where the value of u2.N.YoY is suspect.
#               These are set to NA prior to the call to JAGS
#    bad.u2.A.1   - list of julian weeks where the value of u2.A.1 is suspect. 
#               These are set to NA prior to the call to JAGS
#    bad.u2.N.1   - list of julian weeks where the value of u2.N.1 is suspect.
#               These are set to NA prior to the call to JAGS
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

# force the input vectors to be vectors
time     <- as.vector(time)
n1       <- as.vector(n1)
m2       <- as.vector(m2)
u2.A.YoY <- as.vector(u2.A.YoY)
u2.N.YoY <- as.vector(u2.N.YoY)
u2.A.1   <- as.vector(u2.A.1)
u2.N.1   <- as.vector(u2.N.1)
sampfrac <- as.vector(sampfrac)

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time all match
if(var(c(length(n1),length(m2),length(u2.A.YoY),length(u2.N.YoY),length(u2.A.1),length(u2.N.1),
       length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1, sampfrac, time must all be equal. They are:",
        length(n1)," ",length(m2)," ",length(u2.A.YoY)," ",length(u2.N.YoY)," ",length(u2.A.1)," ",length(u2.N.1),
        length(sampfrac)," ",length(time),"\n")
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
if(any(m2>n1,na.rm=TRUE)){
   cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:",
       paste(n1,collapse=","),"\n m2:",
       paste(m2.collapse=","),"\n")
   return()}

#  3. Elements of bad.m2, bad.u2.A.YoY, bad.u2.A.1, bad.u2.N.YoY, bad.u2.N.1, and hatch.after.YoY must belong to time
if(!all(bad.m2 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",
       paste(bad.m2,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,  collapse=","), "\n")
   return()}
if(!all(bad.u2.A.YoY %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.A.YoY must be elements of strata identifiers. You entered \n bad.u2.A.YoY:",
       paste(bad.u2.A.YoY,collapse=","),"\n Strata identifiers are \n time:",
       paste(time        ,collapse=","), "\n")
   return()}
if(!all(bad.u2.A.1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.A.1 must be elements of strata identifiers. You entered \n bad.u2.A.1:",
       paste(bad.u2.A.1,collapse=","),"\n Strata identifiers are \n time:",
       paste(time      ,collapse=","), "\n")
   return()}
if(!all(bad.u2.N.YoY %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.N.YoY must be elements of strata identifiers. You entered \n bad.u2.N.YoY:",
       paste(bad.u2.N.YoY,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,        collapse=","), "\n")
   return()}
if(!all(bad.u2.N.1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.N.1 must be elements of strata identifiers. You entered \n bad.u2.N.1:",
       paste(bad.u2.N.1,collapse=","),"\n Strata identifiers are \n time:",
       paste(time,      collapse=","), "\n")
   return()}
if(!all(hatch.after.YoY %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** hatch.after.YoY must be elements of strata identifiers. You entered \n hatch.after.YoY:",
   paste(hatch.after.YoY,collapse=","),"\n Strata identifiers are \n time:",
   paste(time,           collapse=","), "\n")
   return()}

#  4. check that strata numbers are contiguous between smallest and largest value of the strata numbers
if( any(seq(min(time),max(time),1) != time,na.rm=TRUE)){
   cat("***** ERROR ***** Strata numbers must be contiguous. \n You entered :", paste(time,collapse=","), "\n")
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

# Deprication of sampling fraction.
if(any(sampfrac != 1)){
  cat("***** ERROR ***** Sampling fraction is depricated for any values other than 1. DO NOT USE ANYMORE. ")
  return()
}

results.filename <- paste(prefix,"-results.txt",sep="")   

stdout <- vector('character')
report <- textConnection('stdout', 'wr', local = TRUE)
sink(report)

cat(paste("Time Stratified Petersen with Diagonal recaptures, error in smoothed U, separating YoY and Age 1 wild and hatchery fish - ", date()))
cat("\nVersion:", version)

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** \n")
temp<- cbind(time, n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1, logitP.cov)
colnames(temp)<- c('time', 'n1','m2','u2.A.YoY', 'u2.N.YoY',"u2.A.1", "u2.N.1", 
                   paste("logitPcov[", 1:ncol(as.matrix(logitP.cov)),"]",sep="") )
print(temp) 
cat("\n\n")
cat("YoY Hatchery fish are released AFTER strata: ", hatch.after.YoY,"\n\n")
cat("YoY  Hatchery fish are clipped at a rate of :", clip.frac.H.YoY,"\n\n")
cat("Age1 Hatchery fish are clipped at a rate of :", clip.frac.H.1  ,"\n\n")
cat("The following strata had m2       set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.A.YoY set to missing: ", 
     if(length(bad.u2.A.YoY)>0){bad.u2.A.YoY} else {" NONE"}, "\n")
cat("The following strata had u2.N.YoY set to missing: ", 
     if(length(bad.u2.N.YoY)>0){bad.u2.N.YoY} else {" NONE"}, "\n")
cat("The following strata had u2.A.1   set to missing: ", 
     if(length(bad.u2.A.1)>0){bad.u2.A.1} else {" NONE"}, "\n")
cat("The following strata had u2.N.1   set to missing: ", 
     if(length(bad.u2.N.1)>0){bad.u2.N.1} else {" NONE"}, "\n")



# Pooled Petersen estimator over ALL of the data including when no releases take place, bad m2, bad.u2.A or bad.u2.N values.
cat("\n\n*** Pooled Petersen Estimate based on pooling over ALL strata***\n\n")
cat("Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",
     sum(u2.A.YoY, na.rm=TRUE)+sum(u2.N.YoY, na.rm=TRUE)+
     sum(u2.A.1  , na.rm=TRUE)+sum(u2.N.1,   na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), 
      sum(u2.A.YoY, na.rm=TRUE)+sum(u2.N.YoY, na.rm=TRUE)+
      sum(u2.A.1  , na.rm=TRUE)+sum(u2.N.1  , na.rm=TRUE))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

# estimate for YoY clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.A.YoY=",  sum(u2.A.YoY, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(n1,       na.rm=TRUE), 
     sum(m2,       na.rm=TRUE), 
     sum(u2.A.YoY, na.rm=TRUE))
cat("Est U.H.YoY(total) ", format(round(pp$U.est)/clip.frac.H.YoY,big.mark=","),
    "  (SE ",              format(round(pp$U.se) /clip.frac.H.YoY,big.mark=","), ")\n")
cat("Est N.H.YoY(total) ", format(round(pp$N.est)/clip.frac.H.YoY,big.mark=","),
    "  (SE ",              format(round(pp$N.se) /clip.frac.H.YoY,big.mark=","), ")\n\n\n")

# estimate for Age1 clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=",  sum(n1,     na.rm=TRUE),
    ";  m2=",     sum(m2,     na.rm=TRUE),
    ";  u2.A.1=", sum(u2.A.1, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(n1,     na.rm=TRUE), 
     sum(m2,     na.rm=TRUE), 
     sum(u2.A.1, na.rm=TRUE))
cat("Est U.H.1(total) ", format(round(pp$U.est)/clip.frac.H.1,big.mark=","),
      "  (SE ",          format(round(pp$U.se) /clip.frac.H.1,big.mark=","), ")\n")
cat("Est N.H.1(total) ", format(round(pp$N.est)/clip.frac.H.1,big.mark=","),
    "  (SE ",            format(round(pp$N.se) /clip.frac.H.1,big.mark=","), ")\n\n\n")

# estimate for YoY wild fish found by subtraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.W.YoY=",  sum((u2.N.YoY+u2.A.YoY-u2.A.YoY/clip.frac.H.YoY), na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum((u2.N.YoY+u2.A.YoY-u2.A.YoY/clip.frac.H.YoY), na.rm=TRUE))
cat("Est U.W.YoY(total) ", format(round(pp$U.est),big.mark=","),
        "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE\n")
cat("Est N.W.YoY(total) ", format(round(pp$N.est),big.mark=","),
        "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE\n\n\n")

# estimate for Age1 wild fish found by subtraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.W.1=",  sum((u2.N.1+u2.A.1-u2.A.1/clip.frac.H.1), na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum((u2.N.1+u2.A.1-u2.A.1/clip.frac.H.1), na.rm=TRUE))
cat("Est U.W.1(total) ", format(round(pp$U.est),big.mark=","),
      "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE\n")
cat("Est N.W.1(total) ", format(round(pp$N.est),big.mark=","),
      "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE\n\n\n")


# Obtain the Pooled Petersen estimator without excluding bad.m2, bad.u2.A.YoY, or bad.u2.N.YoY,
#        bad.u2.A.1, or bad.u2.N.1 values but after removing 0 or NA values
select <- (n1>0) & (!is.na(n1)) & (!is.na(m2)) & (!is.na(u2.A.YoY)) & (!is.na(u2.N.YoY)) &
                                                 (!is.na(u2.A.1))   & (!is.na(u2.N.1))
cat("\n\n*** Pooled Petersen Estimate AFTER excluding bad m2, u2.A.YoY, u2.A.1, u2.N.YoY, or u2.N.1 values  ***\n\n")
cat("The following strata are excluded because n1=0 or NA values in m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1 :", time[!select],"\n\n")

temp.n1       <- n1      [select]
temp.m2       <- m2      [select]
temp.u2.A.YoY <- u2.A.YoY[select]
temp.u2.N.YoY <- u2.N.YoY[select]
temp.u2.A.1   <- u2.A.1  [select]
temp.u2.N.1   <- u2.N.1  [select]

cat("Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2.YoY=",
     sum(temp.u2.A.YoY+temp.u2.N.YoY+
         temp.u2.A.1  +temp.u2.N.1, na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), 
     sum(temp.u2.A.YoY+temp.u2.N.YoY+
         temp.u2.A.1  +temp.u2.N.1  , na.rm=TRUE))
cat("Est U(total) ", format(round(pp$U.est),big.mark=","),"  (SE ", format(round(pp$U.se), big.mark=","), ")\n")
cat("Est N(total) ", format(round(pp$N.est),big.mark=","),"  (SE ", format(round(pp$N.se), big.mark=","), ")\n\n\n")

# estimate for YoY clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.A.YoY=",  sum(temp.u2.A.YoY, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1,       na.rm=TRUE), 
     sum(temp.m2,       na.rm=TRUE), 
     sum(temp.u2.A.YoY, na.rm=TRUE))
cat("Est U.H.YoY(total) ", format(round(pp$U.est)/clip.frac.H.YoY,big.mark=","),
        "  (SE ",          format(round(pp$U.se) /clip.frac.H.YoY,big.mark=","), ")\n")
cat("Est N.H.YoY(total) ", format(round(pp$N.est)/clip.frac.H.YoY,big.mark=","),
        "  (SE ",          format(round(pp$N.se) /clip.frac.H.YoY,big.mark=","), ")\n\n\n")

# estimate for YoY wild fish
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.W.YoY=",  sum((temp.u2.N.YoY+temp.u2.A.YoY-temp.u2.A.YoY/clip.frac.H.YoY), na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction YoY :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum((temp.u2.N.YoY+temp.u2.A.YoY-temp.u2.A.YoY/clip.frac.H.YoY), na.rm=TRUE))
cat("Est U.W.YoY(total) ", format(round(pp$U.est),big.mark=","),
        "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE \n")
cat("Est N.W.YoY(total) ", format(round(pp$N.est),big.mark=","),
        "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE \n\n\n")

# estimate for Age1 clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.A.1=",  sum(temp.u2.A.1, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum(temp.u2.A.1, na.rm=TRUE))
cat("Est U.H.1(total) ", format(round(pp$U.est)/clip.frac.H.1,big.mark=","),
      "  (SE ",          format(round(pp$U.se) /clip.frac.H.1,big.mark=","), ")\n\n\n")
cat("Est N.H.1(total) ", format(round(pp$N.est)/clip.frac.H.1,big.mark=","),
      "  (SE ",          format(round(pp$N.se) /clip.frac.H.1,big.mark=","), ")\n\n\n")

# estimate for Age1 wild fish
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.W.1=",  sum((temp.u2.N.1+temp.u2.A.1-temp.u2.A.1/clip.frac.H.1), na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction 1 :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum((temp.u2.N.1+temp.u2.A.1-temp.u2.A.1/clip.frac.H.1), na.rm=TRUE))
cat("Est U.W.1(total) ", format(round(pp$U.est),big.mark=","),
      "  (SE ",          format(round(pp$U.se) ,big.mark=","), ") APPROXIMATE \n")
cat("Est N.W.1(total) ", format(round(pp$N.est),big.mark=","),
      "  (SE ",          format(round(pp$N.se) ,big.mark=","), ") APPROXIMATE \n\n\n")


# Set the bad values to missing
temp.n1       <- n1
temp.m2       <- m2
temp.u2.A.YoY <- u2.A.YoY
temp.u2.N.YoY <- u2.A.YoY
temp.u2.A.1   <- u2.A.1
temp.u2.N.1   <- u2.A.1

temp.m2      [bad.m2      -min(time)+1] <- NA
temp.u2.A.YoY[bad.u2.A.YoY-min(time)+1] <- NA
temp.u2.N.YoY[bad.u2.N.YoY-min(time)+1] <- NA
temp.u2.A.1  [bad.u2.A.1  -min(time)+1] <- NA
temp.u2.N.1  [bad.u2.N.1  -min(time)+1] <- NA




# Obtain Stratified-Petersen estimator for each stratum after the removal of bad values
cat("*** Stratified Petersen Estimator for each stratum AFTER removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- (temp.u2.A.YoY + temp.u2.N.YoY)
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','(u2.A.YoY+u2.N.YoY)*adj', 'U.YoY[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("Est U.YoY(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
    "  (SE ",        format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum YoY Hatchery YoY AFTER removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- u2.A.YoY
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.A.YoY*adj', 'U.H.YoY[i]', 'SE(U[i])')
print(temp)
cat("** Estimates not adjusted for clip fraction above \n")
cat("Est U.H(total) ", format(round(sum(sp$U.est, na.rm=TRUE)/clip.frac.H.YoY),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))/clip.frac.H.YoY), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum YoY Wild YoY after removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- pmax(0,(u2.N.YoY+u2.A.YoY-u2.A.YoY/clip.frac.H.YoY))
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.YoY-est', 'U.W.YoY[i]', 'SE(U[i])')
print(temp)
cat("Est U.W.YoY(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
        "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))), big.mark=","), ") APPROXIMATE\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Age1 Hatchery AFTER removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- u2.A.1
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.A.1*adj', 'U.H.1[i]', 'SE(U[i])')
print(temp)
cat("** Estimates not adjusted for clip fraction above \n")
cat("Est U.H(total) ", format(round(sum(sp$U.est, na.rm=TRUE)/clip.frac.H.YoY),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$U.se^2, na.rm=TRUE))/clip.frac.H.1), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Age1 Wild after removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- pmax(0,(u2.N.1+u2.A.1-u2.A.1/clip.frac.H.1))
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$U.est), round(sp$U.se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1-est', 'U.W.1[i]', 'SE(U[i])')
print(temp)
cat("Est U.W.1(total) ", format(round(sum(sp$U.est, na.rm=TRUE)),big.mark=","),
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
new.u2.A.YoY   <- rep(0, max(time)-min(time)+1)
new.u2.N.YoY   <- rep(0, max(time)-min(time)+1)
new.u2.A.1     <- rep(0, max(time)-min(time)+1)
new.u2.N.1     <- rep(0, max(time)-min(time)+1)
new.logitP.cov <- matrix(NA, nrow=max(time)-min(time)+1, ncol=ncol(as.matrix(logitP.cov)))
new.time       <- min(time):max(time)


new.n1      [time-min(time)+1]         <- n1
new.m2      [time-min(time)+1]         <- m2
new.m2      [bad.m2-min(time)+1]       <- NA    # wipe out strata where m2 is known to be bad
new.u2.A.YoY[time-min(time)+1]         <- u2.A.YoY
new.u2.A.YoY[bad.u2.A.YoY-min(time)+1] <- NA    # wipe out strata where u2.A is known to be bad
new.u2.N.YoY[time-min(time)+1]         <- u2.N.YoY
new.u2.N.YoY[bad.u2.N.YoY-min(time)+1] <- NA    # wipe out strata where u2.N is known to be bad
new.u2.A.1  [time-min(time)+1]         <- u2.A.1
new.u2.A.1  [bad.u2.A.1  -min(time)+1] <- NA    # wipe out strata where u2.A is known to be bad
new.u2.N.1  [time-min(time)+1]         <- u2.N.1
new.u2.N.1  [bad.u2.N.1  -min(time)+1] <- NA    # wipe out strata where u2.N is known to be bad
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
new.u2.A.YoY <- round(new.u2.A.YoY)
new.u2.N.YoY <- round(new.u2.N.YoY)
new.u2.A.1   <- round(new.u2.A.1  )
new.u2.N.1   <- round(new.u2.N.1  )

# Print out the revised data
hatch.indicator <- rep('   ', max(time)-min(time)+1)
hatch.indicator[hatch.after.YoY-min(time)+1]<- '***'

cat("\n\n*** Revised data *** \n")
temp<- data.frame(time=new.time, n1=new.n1, m2=new.m2, 
       u2.A.YoY=new.u2.A.YoY, u2.N.YoY=new.u2.N.YoY, u2.A.1=new.u2.A.1, u2.N.1=new.u2.N.1,
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
   cat("\nprior to formal call to TimeStratPetersenDiagErrorWHChinook\n")
   browser()
}


if (debug) 
   {results <- TimeStratPetersenDiagErrorWHChinook2(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, 
            u2.A.YoY=new.u2.A.YoY, u2.N.YoY=new.u2.N.YoY, u2.A.1=new.u2.A.1, u2.N.1=new.u2.N.1, 
            hatch.after.YoY=hatch.after.YoY-min(time)+1, 
            clip.frac.H.YoY=clip.frac.H.YoY, clip.frac.H.1=clip.frac.H.1,
            logitP.cov=new.logitP.cov,
            n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,  # set to low values for debugging only
            prior.beta.logitP.mean=prior.beta.logitP.mean, 
            prior.beta.logitP.sd  =prior.beta.logitP.sd,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, InitialSeed=InitialSeed,
            save.output.to.files=save.output.to.files)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenDiagErrorWHChinook2(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, 
            u2.A.YoY=new.u2.A.YoY, u2.N.YoY=new.u2.N.YoY, u2.A.1=new.u2.A.1, u2.N.1=new.u2.N.1, 
            hatch.after.YoY=hatch.after.YoY-min(time)+1, 
            clip.frac.H.YoY=clip.frac.H.YoY, clip.frac.H.1=clip.frac.H.1,
            logitP.cov=new.logitP.cov,
            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims,
            prior.beta.logitP.mean=prior.beta.logitP.mean, 
            prior.beta.logitP.sd  =prior.beta.logitP.sd,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta, 
	          InitialSeed=InitialSeed,
            save.output.to.files=save.output.to.files)
   }

# Now to create the various summary tables of the results

# A plot of the initial and fitted values for YoY, Age 1, and Hatchery chinook
 # In the diagonal case, time, n1, m2, u2 are the same length
  Nstrata <- length(n1)

  plot.df   <- data.frame(time =new.time)
 
  # adjust the u2 for the clipping fractions
  plot.df$n1       <- new.n1
  plot.df$m2       <- new.m2
  plot.df$u2.H.YoY <- new.u2.A.YoY/clip.frac.H.YoY  # only a portion of the hatchery fish are clipped

  plot.df$u2.N.YoY <- new.u2.N.YoY
  plot.df$u2.H.1   <- new.u2.A.1  /clip.frac.H.1  # only a portion of the hatchery fish are clipped
  plot.df$u2.N.1   <- new.u2.N.1
  
  # estimate how many wild fish are present given that only a fraction of hatchery fish are marked
  plot.df$u2.W.1   <- pmax(plot.df$u2.N.1   - plot.df$u2.H.1  *(1-clip.frac.H.1)  ,0) # subtract the guestimated number of hatchery fish
  plot.df$u2.W.YoY <- pmax(plot.df$u2.N.YoY - plot.df$u2.H.YoY*(1-clip.frac.H.YoY),0) # subtract the questimated number of hatchery fish

  plot.df$u2.H.YoY[is.na(plot.df$u2.H.YoY)] <- 1  # in case of missing values
  plot.df$u2.W.YoY[is.na(plot.df$u2.W.YoY)] <- 1  # in case of missing values
  plot.df$u2.H.1  [is.na(plot.df$u2.H.1)  ] <- 1  # in case of missing values
  plot.df$u2.W.1  [is.na(plot.df$u2.W.1)]   <- 1  # in case of missing values

  get.est <- function(est.name, plot.df, hatch.after.YoY, results){
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

      if(est.name=="H.YoY"){
         est.df$logUguess[1:(hatch.after.YoY-min(plot.df$time)+1)]<- NA
         est.df$logU     [1:(hatch.after.YoY-min(plot.df$time)+1)]<- NA
         est.df$logUlcl  [1:(hatch.after.YoY-min(plot.df$time)+1)]<- NA
         est.df$logUucl  [1:(hatch.after.YoY-min(plot.df$time)+1)]<- NA
         est.df$spline   [1:(hatch.after.YoY-min(plot.df$time)+1)]<- NA
      }
      est.df
  }
  plot.data <-rbind( get.est("H.YoY",plot.df, hatch.after.YoY, results),
                     get.est("H.1"  ,plot.df, hatch.after.YoY, results),
                     get.est("W.YoY",plot.df, hatch.after.YoY, results),
                     get.est("W.1"  ,plot.df, hatch.after.YoY, results))
  fit.plot <- ggplot(data=plot.data, aes_(x=~time, color=~group))+
     ggtitle(title, subtitle="Fitted spline curve with 95% credible intervals")+
     geom_point(aes_(y=~logUguess), shape=16, position=position_dodge(width=.2))+  # guesses for population
     geom_point(aes_(y=~logU), shape=19, position=position_dodge(width=.2))+
     geom_line (aes_(y=~logU), position=position_dodge(width=.2))+
     geom_errorbar(aes_(ymin=~logUlcl, ymax=~logUucl), width=.1, position=position_dodge(width=.2))+
     geom_line(aes_(y=~spline),linetype="dashed", position=position_dodge(width=.2)) + 
     xlab("Time Index\nFitted/Smoothed/Raw values plotted for W(black) and H(blue)")+ylab("log(U[i]) + 95% credible interval")+
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


# plot logit P vs time
logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=new.m2, 
             	    u2=u2.A.YoY+u2.N.YoY+u2.A.1+u2.N.1,   logitP.cov=new.logitP.cov, results=results)
if(save.output.to.files)ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
results$plots$logitP.plot <- logitP.plot


# Look at autocorrelation function for Utot.W.YoY, Utot.H.YoY, Utota.W.1, Utot.H.1
mcmc.sample1<- data.frame(parm="Utot.W.YoY", sample=results$sims.matrix[,"Utot.W.YoY"], stringsAsFactors=FALSE)
mcmc.sample2<- data.frame(parm="Utot.H.YoY", sample=results$sims.matrix[,"Utot.H.YoY"], stringsAsFactors=FALSE)
mcmc.sample3<- data.frame(parm="Utot.W.1",   sample=results$sims.matrix[,"Utot.W.1"], stringsAsFactors=FALSE)
mcmc.sample4<- data.frame(parm="Utot.H.1",   sample=results$sims.matrix[,"Utot.H.1"], stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2, mcmc.sample3, mcmc.sample4)
acf.Utot.plot <- plot_acf(mcmc.sample)
if(save.output.to.files)ggsave(plot=acf.Utot.plot, filename=paste(prefix,"-Utot-acf.pdf",sep=""), height=4, width=6, units="in")
results$plots$acf.Utot.plot <- acf.Utot.plot


# Look at the shape of the posterior distribution
mcmc.sample1<- data.frame(parm="Utot.W.YoY", sample=results$sims.matrix[,"Utot.W.YoY"], stringsAsFactors=FALSE)
mcmc.sample2<- data.frame(parm="Utot.H.YoY", sample=results$sims.matrix[,"Utot.H.YoY"], stringsAsFactors=FALSE)
mcmc.sample3<- data.frame(parm="Utot.W.1",   sample=results$sims.matrix[,"Utot.W.1"],   stringsAsFactors=FALSE)
mcmc.sample4<- data.frame(parm="Utot.H.1",   sample=results$sims.matrix[,"Utot.H.1"],   stringsAsFactors=FALSE)
mcmc.sample <- rbind(mcmc.sample1, mcmc.sample2, mcmc.sample3, mcmc.sample4)
post.Utot.plot <- plot_posterior(mcmc.sample, ncol=2)
post.Utot.plot
if(save.output.to.files)ggsave(plot=post.Utot.plot, filename=paste(prefix,"-Utot-posterior.pdf",sep=""),
                               height=ifelse(length(unique(mcmc.sample$parm))<=2,4,6), width=6, units="in")
results$plots$post.Utot.plot <- post.Utot.plot



# make the Bayesian predictive distribution (Bayesian p-value plots)
#browser()
discrep <-PredictivePosterior.TSPDE.WHCH2 (time, new.n1, new.m2,   # get the discrepancy measures
          new.u2.A.YoY, new.u2.N.YoY, new.u2.A.1, new.u2.N.1, 
          clip.frac.H.YoY, clip.frac.H.1, 
          expit(results$sims.list$logitP), 
          round(results$sims.list$U.W.YoY), 
          round(pmax(results$sims.list$U.H.YoY,0)), 
          round(results$sims.list$U.W.1), 
          round(results$sims.list$U.H.1), 
          hatch.after.YoY) #don't forget that hatchery fish is 0 until hatch.after
#browser()
gof <- PredictivePosteriorPlot.TSPDE.WHCH2 (discrep)
if(save.output.to.files){
  pdf(file=paste(prefix,"-GOF.pdf",sep=""))
    plyr::l_ply(gof, function(x){plot(x)})
  dev.off()
}
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

# Summary of population sizes. Add pretty printing for the final results
cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
cat("Wild YoY \n")
temp<- results$summary[ grep("Utot.W.YoY", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nWild Age 1 \n")
temp<-results$summary[ grep("Utot.W.1"  , rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nHatchery YoY\n")
temp<- results$summary[ grep("Utot.H.YoY", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nHatchery Age 1\n")
temp<-results$summary[ grep("Utot.H.1"   , rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nGrand Total\n")
temp<- results$summary[ rownames(results$summary) == "Utot",]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)


#browser()
time.H <- time>hatch.after.YoY
cat("\n\n\n\n*** Summary of Quantiles of Run Timing.Wild *** \n")
cat(    "    This is based on the sample weeks provided and the U.W.YoY[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W.YoY, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))
cat(    "\n    This is based on the sample weeks provided and the U.W.1[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W.1, prob=run.prob)

temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))



cat("\n\n*** Summary of Quantiles of Run Timing.Hatchery *** \n")
cat(    "    This is based on the sample weeks provided and the U.H.YoY[i] values \n") 
q <- RunTime(time=time[time.H], U=results$sims.list$U.H.YoY[,time.H], prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))
cat(    "\n    This is based on the sample weeks provided and the U.H.1[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.H.1, prob=run.prob)
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
results$data <- list( time=time, n1=n1, m2=m2, 
                      u2.A.YoY=u2.A.YoY, u2.N.YoY=u2.N.YoY, u2.A.1=u2.A.1, u2.N.1=u2.N.1, 
                      clip.frac.H.YoY=clip.frac.H.YoY, clip.frac.H.1=clip.frac.H.1,
                      hatch.after.YoY=hatch.after.YoY, 
                      bad.m2=bad.m2,
                      bad.u2.A.YoY=bad.u2.A.YoY, bad.u2.N.YoY=bad.u2.N.YoY, 
                      bad.u2.A.1=bad.u2.A.1, bad.u2.N.1=bad.u2.N.1, 
                      logitP.cov=logitP.cov,
                      version=version, date_run=date(), title=title)

return(results)
} # end of function
