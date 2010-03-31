# 2010-03-29 CJS Inital version of code

TimeStratPetersenDiagErrorWHChinook2_fit<- 
       function( title="TSPDE-WHChinook2", prefix="TSPDE-WHChinook2-", 
                 time, n1, m2, 
                 u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1,
                 clip.frac.H.YoY, clip.frac.H.1, sampfrac, 
                 hatch.after.YoY=NULL, 
                 bad.m2=c(), bad.u2.A.YoY=c(), bad.u2.N.YoY=c(), bad.u2.A.1=c(), bad.u2.N.1=c(),
                 logitP.cov=rep(1,length(n1)),
                 n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
                 tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
                 mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)), 
                 tau_xiP=1/var(logit((m2+.5)/(n1+1)), na.rm=TRUE), 
                 tauP.alpha=.001, tauP.beta=.001,
                 run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
                 debug=FALSE, debug2=FALSE, openbugs=TRUE,
                 OPENBUGS.directory=file.path("c:","Program Files","OpenBugs"),
                 WINBUGS.directory=file.path("c:","Program Files","WinBUGS14"),
                 InitialSeed=trunc(runif(1,min=1, max=1000000000))) {
# Fit a Time Stratified Petersen model with diagonal entries and with smoothing on U allowing for random error,
# covariates for the the capture probabilities, and separating the YoY and Age1 wild vs hatchery fish
# The "diagonal entries" implies that no marked fish are recaptured outside the (time) stratum of release
#
   version <- '2010-03-28'

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
#    sampfrac - sampling fraction to adjust for how many days of the week was the trap operating
#              This is expressed as fraction i.e. 3 days out of 7 is expressed as 3/7=.42 etc.
#              If the trap was operating ALL days, then the SampFrac = 1. It is possible for the sampling
#              fraction to be > 1 (e.g. a mark is used for 8 days instead of 7. The data are adjusted
#              back to a 7 day week as well.
#    hatch.after - julian week AFTER which hatchery fish are released 
#    bad.m2  - list of julian numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to WinBugs/OpenBugs
#    bad.u2.A.YoY - list of julian weeks where the value of u2.A.YoY is suspect. 
#               These are set to NA prior to the call to WinBugs/OpenBugs
#    bad.u2.N.YoY - list of julian weeks where the value of u2.N.YoY is suspect.
#               These are set to NA prior to the call to WinBugs/OpenBugs
#    bad.u2.A.1   - list of julian weeks where the value of u2.A.1 is suspect. 
#               These are set to NA prior to the call to WinBugs/OpenBugs
#    bad.u2.N.1   - list of julian weeks where the value of u2.N.1 is suspect.
#               These are set to NA prior to the call to WinBugs/OpenBugs
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
#             and WINBUGS will run and stop waiting for your to exit and complete
#    openbugs - if TRUE, then openbugs is called; else WInbugs is called

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time all match
if(var(c(length(n1),length(m2),length(u2.A.YoY),length(u2.N.YoY),length(u2.A.1),length(u2.N.1),
       length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1, sampfrac, time must all be equal. They are:",
        length(n1),length(m2),length(u2.A.YoY),length(u2.N.YoY),length(u2.A.1),length(u2.N.1),
        length(sampfrac),length(time),"\n")
   return()}
if(length(logitP.cov) %% length(n1) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of n1 etc They are:",
        length(n1),length(logitP.cov),dim(logitP.cov),"\n")
   return()}
#  2. Check that m2<= n1
if(any(m2>n1)){
   cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:",n1,"\n m2:",m2,"\n")
   return()}
#  3. Elements of bad.m2, bad.u2.A.YoY, bad.u2.A.1, bad.u2.N.YoY, bad.u2.N.1, and hatch.after.YoY must belong to time
if(!all(bad.m2 %in% time)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",bad.m2,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.A.YoY %in% time)){
   cat("***** ERROR ***** bad.u2.A.YoY must be elements of strata identifiers. You entered \n bad.u2.A.YoY:",
       bad.u2.A.YoY,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.A.1 %in% time)){
   cat("***** ERROR ***** bad.u2.A.1 must be elements of strata identifiers. You entered \n bad.u2.A.1:",
       bad.u2.A.1,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.N.YoY %in% time)){
   cat("***** ERROR ***** bad.u2.N.YoY must be elements of strata identifiers. You entered \n bad.u2.N.YoY:",
       bad.u2.N.YoY,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.N.1 %in% time)){
   cat("***** ERROR ***** bad.u2.N.1 must be elements of strata identifiers. You entered \n bad.u2.N.1:",
       bad.u2.N.1,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(hatch.after.YoY %in% time)){
   cat("***** ERROR ***** hatch.after.YoY must be elements of strata identifiers. You entered \n hatch.after.YoY:",
   hatch.after.YoY,"\n Strata identifiers are \n time:",time, "\n")
   return()}
#  4. Check for existence of the openbugs/winbugs directories. This is not a guarantee that openbugs/wingbugs exists, but it is a start 
if(openbugs){ # user wants to use OpenBUGS
   if(file.access(OPENBUGS.directory) != 0){
      cat("***** ERROR ***** Can't find OPENBUGS directory. You entered", OPENBUGS.directory, "\n")
      return() }
} else {  # user wants to use WinBUGS
  if(file.access(WINBUGS.directory) != 0){
      cat("***** ERROR ***** Can't find WINBUGS directory. You entered", WINBUGS.directory, "\n")
      return() }
}
#  5. check that strata numbers are contiguous between smallest and largest value of the strata numbers
if( any(seq(min(time),max(time),1) != time)){
   cat("***** ERROR ***** Strata numbers must be contiguous. \n You entered :", time, "\n")
   return()
}


results.filename <- paste(prefix,"-results.txt",sep="")   


sink(results.filename)
cat(paste("Time Stratified Petersen with Diagonal recaptures, error in smoothed U, separating YoY and Age 1 wild and hatchery fish - ", date()))
cat("\nVersion:", version)

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** \n")
temp<- cbind(time, n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1, round(sampfrac,digits=2), logitP.cov)
colnames(temp)<- c('time', 'n1','m2','u2.A.YoY', 'u2.N.YoY',"u2.A.1", "u2.N.1", 
                   'SampFrac', paste("logitPcov[", 1:ncol(as.matrix(logitP.cov)),"]",sep="") )
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
cat("\n\n*** Pooled Petersen Estimate based on pooling over ALL strata adjusting for sampling fraction***\n\n")
cat("Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",
     sum(u2.A.YoY/sampfrac, na.rm=TRUE)+sum(u2.N.YoY/sampfrac, na.rm=TRUE)+
     sum(u2.A.1  /sampfrac, na.rm=TRUE)+sum(u2.N.1/sampfrac,   na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), 
      sum(u2.A.YoY/sampfrac, na.rm=TRUE)+sum(u2.N.YoY/sampfrac, na.rm=TRUE)+
      sum(u2.A.1  /sampfrac, na.rm=TRUE)+sum(u2.N.1  /sampfrac, na.rm=TRUE))
cat("Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")
# estimate for YoY clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.A.YoY=",  sum(u2.A.YoY/sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum(u2.A.YoY/sampfrac, na.rm=TRUE))
cat("Est U.H.YoY(total) ", format(round(pp$est)/clip.frac.H.YoY,big.mark=","),
    "  (SE ",          format(round(pp$se) /clip.frac.H.YoY,big.mark=","), ")\n\n\n")
# estimate for Age1 clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.A.1=",  sum(u2.A.1/sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum(u2.A.1/sampfrac, na.rm=TRUE))
cat("Est U.H.1(total) ", format(round(pp$est)/clip.frac.H.1,big.mark=","),
      "  (SE ",          format(round(pp$se) /clip.frac.H.1,big.mark=","), ")\n\n\n")

# estimate for YoY wild fish found by subtraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.W.YoY=",  sum((u2.N.YoY+u2.A.YoY-u2.A.YoY/clip.frac.H.YoY)/sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum((u2.N.YoY+u2.A.YoY-u2.A.YoY/clip.frac.H.YoY)/sampfrac, na.rm=TRUE))
cat("Est U.W.YoY(total) ", format(round(pp$est),big.mark=","),
        "  (SE ",          format(round(pp$se) ,big.mark=","), ") APPROXIMATE\n\n\n")
# estimate for Age1 wild fish found by subtraction
cat("Total n1=", sum(n1, na.rm=TRUE),
    ";  m2=",    sum(m2, na.rm=TRUE),
    ";  u2.W.1=",  sum((u2.N.1+u2.A.1-u2.A.1/clip.frac.H.1)/sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(n1, na.rm=TRUE), 
     sum(m2, na.rm=TRUE), 
     sum((u2.N.1+u2.A.1-u2.A.1/clip.frac.H.1)/sampfrac, na.rm=TRUE))
cat("Est U.W.1(total) ", format(round(pp$est),big.mark=","),
      "  (SE ",          format(round(pp$se) ,big.mark=","), ") APPROXIMATE\n\n\n")


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
temp.sampfrac <- sampfrac[select]

cat("Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2.YoY=",
     sum(temp.u2.A.YoY/temp.sampfrac+temp.u2.N.YoY/temp.sampfrac+
         temp.u2.A.1  /temp.sampfrac+temp.u2.N.1  /temp.sampfrac, na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), 
     sum(temp.u2.A.YoY/temp.sampfrac+temp.u2.N.YoY/temp.sampfrac+
         temp.u2.A.1  /temp.sampfrac+temp.u2.N.1  /temp.sampfrac, na.rm=TRUE))
cat("Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

# estimate for YoY clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.A.YoY=",  sum(temp.u2.A.YoY/temp.sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum(temp.u2.A.YoY/temp.sampfrac, na.rm=TRUE))
cat("Est U.H.YoY(total) ", format(round(pp$est)/clip.frac.H.YoY,big.mark=","),
    "  (SE ",          format(round(pp$se) /clip.frac.H.YoY,big.mark=","), ")\n\n\n")
# estimate for YoY wild fish
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.W.YoY=",  sum((temp.u2.N.YoY+temp.u2.A.YoY-temp.u2.A.YoY/clip.frac.H.YoY)/temp.sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction YoY :", clip.frac.H.YoY, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum((temp.u2.N.YoY+temp.u2.A.YoY-temp.u2.A.YoY/clip.frac.H.YoY)/temp.sampfrac, na.rm=TRUE))
cat("Est U.W.YoY(total) ", format(round(pp$est),big.mark=","),
    "  (SE ",          format(round(pp$se) ,big.mark=","), ") APPROXIMATE \n\n\n")

# estimate for Age1 clipped fish (hatchery) and expand by the clip fraction
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.A.1=",  sum(temp.u2.A.1/temp.sampfrac, na.rm=TRUE),"\n")
cat("Clip fraction :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum(temp.u2.A.1/temp.sampfrac, na.rm=TRUE))
cat("Est U.H.1(total) ", format(round(pp$est)/clip.frac.H.1,big.mark=","),
    "  (SE ",          format(round(pp$se) /clip.frac.H.1,big.mark=","), ")\n\n\n")
# estimate for Age1 wild fish
cat("Total n1=", sum(temp.n1, na.rm=TRUE),
    ";  m2=",    sum(temp.m2, na.rm=TRUE),
    ";  u2.W.1=",  sum((temp.u2.N.1+temp.u2.A.1-temp.u2.A.1/clip.frac.H.1)/temp.sampfrac, na.rm=TRUE),
    "[Formed by interpolation based on clip rate]\n")
cat("Clip fraction 1 :", clip.frac.H.1, "\n\n")
pp <- SimplePetersen(
     sum(temp.n1, na.rm=TRUE), 
     sum(temp.m2, na.rm=TRUE), 
     sum((temp.u2.N.1+temp.u2.A.1-temp.u2.A.1/clip.frac.H.1)/temp.sampfrac, na.rm=TRUE))
cat("Est U.W.1(total) ", format(round(pp$est),big.mark=","),
    "  (SE ",          format(round(pp$se) ,big.mark=","), ") APPROXIMATE \n\n\n")


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
temp.u2 <- (temp.u2.A.YoY + temp.u2.N.YoY)/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','(u2.A.YoY+u2.N.YoY)*adj', 'U.YoY[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("Est U.YoY(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ",        format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum YoY Hatchery YoY AFTER removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- u2.A.YoY/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.A.YoY*adj', 'U.H.YoY[i]', 'SE(U[i])')
print(temp)
cat("** Estimates not adjusted for clip fraction above \n")
cat("Est U.H(total) ", format(round(sum(sp$est, na.rm=TRUE)/clip.frac.H.YoY),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$se^2, na.rm=TRUE))/clip.frac.H.YoY), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum YoY Wild YoY after removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- pmax(0,(u2.N.YoY+u2.A.YoY-u2.A.YoY/clip.frac.H.YoY)/sampfrac)
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.YoY-est', 'U.W.YoY[i]', 'SE(U[i])')
print(temp)
cat("Est U.W.YoY(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ") APPROXIMATE\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Age1 Hatchery AFTER removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- u2.A.1/sampfrac
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.A.1*adj', 'U.H.1[i]', 'SE(U[i])')
print(temp)
cat("** Estimates not adjusted for clip fraction above \n")
cat("Est U.H(total) ", format(round(sum(sp$est, na.rm=TRUE)/clip.frac.H.YoY),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$se^2, na.rm=TRUE))/clip.frac.H.1), big.mark=","), ")\n\n\n")

cat("*** Stratified Petersen Estimator for each stratum Age1 Wild after removing bad m2 values after adjusting for sampling fration ***\n\n")
temp.u2 <- pmax(0,(u2.N.1+u2.A.1-u2.A.1/clip.frac.H.1)/sampfrac)
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1-est', 'U.W.1[i]', 'SE(U[i])')
print(temp)
cat("Est U.W.1(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ",          format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ") APPROXIMATE\n\n\n")



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
new.sampfrac   <- rep(0, max(time)-min(time)+1)
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
new.u2.A.YoY <- round(new.u2.A.YoY/new.sampfrac)
new.u2.N.YoY <- round(new.u2.N.YoY/new.sampfrac)
new.u2.A.1   <- round(new.u2.A.1  /new.sampfrac)
new.u2.N.1   <- round(new.u2.N.1  /new.sampfrac)

# Adjust for strata where sampling fraction=0. On these strata
# u2.A and u2.N is set to NA so that there is NO information on U2 for this stratum
new.u2.A.YoY[new.sampfrac<.001] <- NA
new.u2.N.YoY[new.sampfrac<.001] <- NA
new.u2.A.1  [new.sampfrac<.001] <- NA
new.u2.N.1  [new.sampfrac<.001] <- NA

# Print out the revised data
hatch.indicator <- rep('   ', max(time)-min(time)+1)
hatch.indicator[hatch.after.YoY-min(time)+1]<- '***'

cat("\n\n*** Revised data *** \n")
temp<- data.frame(time=new.time, n1=new.n1, m2=new.m2, 
       u2.A.YoY=new.u2.A.YoY, u2.N.YoY=new.u2.N.YoY, u2.A.1=new.u2.A.1, u2.N.1=new.u2.N.1,
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
            n.chains=3, n.iter=2000, n.burnin=300, n.sims=300,  # set to low values for debugging only
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug,  openbugs=openbugs, InitialSeed=InitialSeed ,
            OPENBUGS.directory=OPENBUGS.directory, WINBUGS.directory=WINBUGS.directory)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenDiagErrorWHChinook2(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, 
            u2.A.YoY=new.u2.A.YoY, u2.N.YoY=new.u2.N.YoY, u2.A.1=new.u2.A.1, u2.N.1=new.u2.N.1, 
            hatch.after.YoY=hatch.after.YoY-min(time)+1, 
            clip.frac.H.YoY=clip.frac.H.YoY, clip.frac.H.1=clip.frac.H.1,
            logitP.cov=new.logitP.cov,
            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims,
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            openbugs=openbugs, InitialSeed=InitialSeed,
            OPENBUGS.directory=OPENBUGS.directory, WINBUGS.directory=WINBUGS.directory)
   }

# Now to create the various summary tables of the results

#browser()
# A plot of the observered log(U) on the log scale, and the final mean log(U)
plot_logU <- function(title, time, n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1, clip.frac.H.YoY, clip.frac.H.1,
             hatch.after.YoY, results){
#  Plot the observed and fitted logU values along with posterior limits
#  n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1 are the raw data
#  results is the summary table from WinBugs

   Nstrata <- length(n1)

# get rough guesses as to the number of hatchery and wildfish in the sample.
   u2.H.YoY <- u2.A.YoY/clip.frac.H.YoY  # only a portion of the YoY hatchery fish are clipped
   u2.W.YoY <- pmax(u2.N.YoY - u2.H.YoY*(1-clip.frac.H.YoY),0) # subtract the questimated number of hatchery fish
   u2.H.YoY[is.na(u2.H.YoY)] <- 1  # in case of missing values
   u2.W.YoY[is.na(u2.W.YoY)] <- 1  # in case of missing values
   # Estimate number of Age1 wild and hatchery fish based on clip rate
   u2.H.1   <- u2.A.1/clip.frac.H.1  # only a portion of the AGE1 hatchery fish are clipped
   u2.W.1   <- pmax(u2.N.1 - u2.H.1*(1-clip.frac.H.1),0) # subtract the questimated number of hatchery fish
   u2.H.1[is.na(u2.H.1)] <- 1  # in case of missing values
   u2.W.1[is.na(u2.W.1)] <- 1  # in case of missing values

   avg.P <- sum(m2,na.rm=TRUE)/sum(n1, na.rM=TRUE)
   Uguess.W.YoY <- pmax((u2.W.YoY+1)*(n1+2)/(m2+1), u2.W.YoY/avg.P, na.rm=TRUE)  # try and keep Uguess larger than observed values
   Uguess.H.YoY <- pmax((u2.H.YoY+1)*(n1+2)/(m2+1), u2.H.YoY/avg.P, na.rm=TRUE)
   Uguess.H.YoY[1:hatch.after.YoY] <- 0   # no YoY hatchery fish prior to release from hatchery
   Uguess.W.1   <- pmax((u2.W.1+1)*(n1+2)/(m2+1), u2.W.1/avg.P, na.rm=TRUE)  # try and keep Uguess larger than observed values
   Uguess.H.1   <- pmax((u2.H.1+1)*(n1+2)/(m2+1), u2.H.1/avg.P, na.rm=TRUE)


   min_logU <- min( log(c(Uguess.W.YoY,Uguess.H.YoY,Uguess.W.1,Uguess.H.1)+1), na.rm=TRUE)
   max_logU <- max( log(c(Uguess.W.YoY,Uguess.H.YoY,Uguess.W.1,Uguess.H.1)+1), na.rm=TRUE)

   # which rows contain the etaU.W[xx] and etaU.H[xx] ?
   results.row.names <- rownames(results$summary)
   etaU.row.index.W.YoY <- grep("etaU.W.YoY", results.row.names)
   etaU.row.index.H.YoY <- grep("etaU.H.YoY", results.row.names)
   etaU.row.index.W.1   <- grep("etaU.W.1"  , results.row.names)
   etaU.row.index.H.1   <- grep("etaU.H.1"  , results.row.names)
   etaU.W.YoY  <- results$summary[etaU.row.index.W.YoY,]
   etaU.H.YoY  <- results$summary[etaU.row.index.H.YoY,]
   etaU.W.1    <- results$summary[etaU.row.index.W.1,]
   etaU.H.1    <- results$summary[etaU.row.index.H.1,]
 
   min_logU <- min( c(min_logU, etaU.W.YoY[,"mean"], etaU.H.YoY[,"mean"]), na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.1  [,"mean"], etaU.H.1  [,"mean"]), na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.YoY[,"mean"], etaU.H.YoY[,"mean"]), na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.1  [,"mean"], etaU.H.1  [,"mean"]), na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.YoY[,"2.5%"], etaU.H.YoY[,"2.5%"]), na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.1  [,"2.5%"], etaU.H.1  [,"2.5%"]), na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.YoY[,"2.5%"], etaU.H.YoY[,"2.5%"]), na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.1  [,"2.5%"], etaU.H.1  [,"2.5%"]), na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.YoY[,"97.5%"],etaU.H.YoY[,"97.5%"]),na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.1  [,"97.5%"],etaU.H.1  [,"97.5%"]),na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.YoY[,"97.5%"],etaU.H.YoY[,"97.5%"]),na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.1  [,"97.5%"],etaU.H.1  [,"97.5%"]),na.rm=TRUE)
   min_logU <- max(min_logU, -1)  # no point in plotting values less than 1 for abundance

   #browser()
   # the wild fish are estimated for ALL strata
   # the hatchery fish are estimated ONLY for those weeks after hatch.after
   time.W <- time
   time.H <- time[time>hatch.after.YoY]
   # plot the raw log(U) values
   plot(time.W, log(Uguess.W.YoY), type="p", pch="w", 
        main=paste(title,"\nFitted spline curve to raw U.W[i] U.H[i] with 95% credible intervals"),
        sub='Fitted/Smoothed/Raw values plotted for W(solid) and H(dash)',
        ylab='log(U[i])',
        xlab='Time Index', ylim=c(min_logU,max_logU))  # initial points on log scale.
   points(time.H, log(Uguess.H.YoY[time>hatch.after.YoY]), pch="h")
   points(time.W, log(Uguess.W.1)                    , pch="W")
   points(time.W, log(Uguess.H.1)                    , pch="H")

   # plot the mean of the etaU
   points(time.W,    etaU.W.YoY[,"mean"], type="p", pch=19)  # fitted values for wild fish
   points(time.H+.1, etaU.H.YoY[,"mean"], type="p", pch=22)  # fitted values for hatchery fish
   points(time.W-.1, etaU.H.1  [,"mean"], type="p", pch=19)  # fitted values for hatchery fish
   points(time.W-.2, etaU.W.1  [,"mean"], type="p", pch=22)  # fitted values for hatchery fish
   lines(time.W,    etaU.W.YoY[,"mean"])  # add smoothed spline through points
   lines(time.H+.1, etaU.H.YoY[,"mean"], lty=2)
   lines(time.W-.1, etaU.H.1  [,"mean"])  # add smoothed spline through points
   lines(time.W-.2, etaU.W.1  [,"mean"], lty=2)
   # plot the 2.5 -> 97.5 posterior values
   segments(time.W,    etaU.W.YoY[,"2.5%"], time.W,    etaU.W.YoY[,"97.5%"])
   segments(time.H+.1, etaU.H.YoY[,"2.5%"], time.H+.1, etaU.H.YoY[,"97.5%"], lty=2)
   segments(time.W-.1, etaU.H.1  [,"2.5%"], time.W-.1, etaU.W.1  [,"97.5%"])
   segments(time.W-.2, etaU.W.1  [,"2.5%"], time.W-.2, etaU.H.1  [,"97.5%"], lty=2)

   # plot the spline curve before the error term is added.
   # extract the bU coefficients
   logUne.row.index.W.YoY <- grep("logUne.W.YoY", results.row.names)
   logUne.row.index.H.YoY <- grep("logUne.H.YoY", results.row.names)
   logUne.row.index.W.1   <- grep("logUne.W.1"  , results.row.names)
   logUne.row.index.H.1   <- grep("logUne.H.1"  , results.row.names)
   logUne.W.YoY<- results$summary[logUne.row.index.W.YoY,"mean"]
   logUne.H.YoY<- results$summary[logUne.row.index.H.YoY,"mean"]
   logUne.W.1  <- results$summary[logUne.row.index.W.1  ,"mean"]
   logUne.H.1  <- results$summary[logUne.row.index.H.1  ,"mean"]
   lines(time.W, logUne.W.YoY, lty=1)  # plot the curve
   lines(time.H, logUne.H.YoY, lty=2)  # plot the curve # too busy to plot the spline points as well
   lines(time.W, logUne.W.1  , lty=1)  # plot the curve
   lines(time.W, logUne.H.1  , lty=2)  # plot the curve # too busy to plot the spline points as well

}

plot_logitP <- function(title, time, n1, m2, logitP.cov, results){
#  Plot the observed and fitted logit(p) values along with posterior limits
#  n1, m2, u2.A, u2.N are the raw data (u2 has been adjusted upward for sampling fraction < 1 prior to call)
#  logitP.cov is the covariate matrix for modelling the logit(P)'s
#  results is the summary table from WinBugs
#
#  Get the minimum and maximum values for the axis on the logit scale

   min_logitP <- 100
   max_logitP <- -100

   if(debug2){
      cat("plot_logitP\n")
      browser()
   }
   Nstrata <- length(n1)
   raw_logitP <- logit((m2+1)/(n1+2))        # based on raw data
   min_logitP <- min( c(min_logitP, raw_logitP), na.rm=TRUE)
   max_logitP <- max( c(max_logitP, raw_logitP), na.rm=TRUE)

   # which rows contain the logitP[xx] ?
   results.row.names <- rownames(results$summary)
   logitP.row.index    <- grep("^logitP", results.row.names)
   est_logitP<- results$summary[logitP.row.index,]
   min_logitP <- min( c(min_logitP, est_logitP[,"mean"]), na.rm=TRUE)
   max_logitP <- max( c(max_logitP, est_logitP[,"mean"]), na.rm=TRUE)
   min_logitP <- min( c(min_logitP, est_logitP[,"2.5%"]), na.rm=TRUE)
   max_logitP <- max( c(max_logitP, est_logitP[,"2.5%"]), na.rm=TRUE)
   min_logitP <- min( c(min_logitP, est_logitP[,"97.5%"]),na.rm=TRUE)
   max_logitP <- max( c(max_logitP, est_logitP[,"97.5%"]),na.rm=TRUE)



   main.title <- paste(title,"\nPlot of logit(p[i]) with 95% credible intervals")
   if(ncol(logitP.cov)>1){main.title<- title}
       sub.title <- paste("Horizontal line is estimated beta.logitP[1]",
                 "\nInner fence is c.i. on beta.logitP[1]",
                 "\nOuter fence is 95% range on logit(p)")
   if(ncol(logitP.cov)>1){sub.title <- "Dashed line is second covariate"}
   plot(time, raw_logitP, 
       main=main.title,
       sub=sub.title,
       ylab='logit(p[i])', xlab='Stratum', ylim=c(min_logitP,max_logitP))  # initial points on log scale.

 
   # plot the posterior mean of the logitP if there is only one column for a covariate
   points(time, est_logitP[,"mean"],type="p", pch=19) # the final estimates
   lines(time, est_logitP[,"mean"])  # join the mean of the fitted logitP
 
   # plot the 2.5 -> 97.5 posterior values
   segments(time, est_logitP[,"2.5%"], time, est_logitP[,"97.5%"])

   if(ncol(logitP.cov)==1){  # if only 1 column for covariate vector, usually an intercept
      # plot the posterior mean of the beta.logitP[1] term which is usually
      #      the intercept in most models with covariates along with 95% credible interval
      intercept.row.index    <- grep("beta.logitP[1]", results.row.names, fixed=TRUE)
      intercept <- results$summary[intercept.row.index,]
      segments(time[1], intercept["mean"],time[Nstrata], intercept["mean"])
      segments(time[1], intercept["2.5%"],time[Nstrata], intercept["2.5%"], lty=2)
      segments(time[1], intercept["97.5%"],time[Nstrata],intercept["97.5%"], lty=2)

      # plot the posterior "95% range" for the logit(P)'s based on N(xip, sigmaP^2)
      sigmaP.row.index <- grep("sigmaP", results.row.names)
      sigmaP <- results$summary[sigmaP.row.index,]
      segments(time[1], intercept["mean"]-2*sigmaP["mean"], time[Nstrata], intercept["mean"]-2*sigmaP["mean"], lty=3)
      segments(time[1], intercept["mean"]+2*sigmaP["mean"], time[Nstrata], intercept["mean"]+2*sigmaP["mean"], lty=3)
   }
   if(ncol(logitP.cov)>1){  # if exactly 2 covariates, plot the second covarite over time as well
      par(new=TRUE)   # reuse the same plot
      plot(time, logitP.cov[,2], type="l", lty=2, axes=FALSE, xlab="", ylab="")  # plot the covariate
   }

   # plot residuals of the logit(P)'s against the various covariates
   # to be done in my next life
}


pdf(file=paste(prefix,"-logU.pdf",sep=""))
plot_logU(title=title, time=new.time, n1=new.n1, m2=new.m2, 
          u2.A.YoY=new.u2.A.YoY, u2.N.YoY=new.u2.N.YoY, u2.A.1=new.u2.A.1, u2.N.1=new.u2.N.1, 
          clip.frac.H.YoY, clip.frac.H.1, hatch.after.YoY=hatch.after.YoY, results=results)
dev.off()

pdf(file=paste(prefix,"-logitP.pdf",sep=""))
plot_logitP(title=title, time=new.time, n1=new.n1, m2=new.m2, logitP.cov=new.logitP.cov, results=results) 
dev.off()

# Look at autocorrelation function for Utot.W.YoY, Utot.H.YoY, Utota.W.1, Utot.H.1
# This is plotted as a split screen with 4 cells
pdf(file=paste(prefix,"-Utot-acf.pdf",sep=""))
split.screen(c(2,2))
screen(1)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
acf(results$sims.matrix[,"Utot.W.YoY"], main=paste(title,"\nACF U total.W.YoY"))
screen(2)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
acf(results$sims.matrix[,"Utot.H.YoY"], main=paste(title,"\nACF U total.H.YoY"))
screen(3)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
acf(results$sims.matrix[,"Utot.W.1"], main=paste(title,"\nACF U total.W.1"))
screen(4)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
acf(results$sims.matrix[,"Utot.H.1"], main=paste(title,"\nACF function for U total.H.1"))
close.screen(all.screen=TRUE)
dev.off()

# Look at the shape of the posterior distribution
# This is plotted as a split screen with 4 cells
pdf(file=paste(prefix,"-Utot-posterior.pdf",sep=""))
split.screen(c(2,2))
screen(1)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
plot( x=density(as.vector(results$sims.array[,,"Utot.W.YoY"])), 
    main=paste(title,'\nPosterior of U-total.W.YoY'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.W.YoY",c("2.5%","97.5%")])  # add vertical reference lines
screen(2)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
plot( x=density(as.vector(results$sims.array[,,"Utot.H.YoY"])), 
    main=paste(title,'\nPosterior of U-total.H.YoY'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.H.YoY",c("2.5%","97.5%")])  # add vertical reference lines
screen(3)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
plot( x=density(as.vector(results$sims.array[,,"Utot.W.1"])), 
    main=paste(title,'\nPosterior of U-total.W.1'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.W.1",c("2.5%","97.5%")])  # add vertical reference lines
screen(4)
par(cex=.5)
par(mai=c(.40,.40,.40,.40))
plot( x=density(as.vector(results$sims.array[,,"Utot.H.1"])), 
    main=paste(title,'\nPosterior of U-total.H.1'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.H.1",c("2.5%","97.5%")])  # add vertical reference lines
close.screen(all.screen=TRUE)
dev.off()


# make the Bayesian predictive distribution (Bayesian p-value plots)
pdf(file=paste(prefix,"-GOF.pdf",sep=""))
#browser()
discrep <-PredictivePosterior.TSPDE.WHCH2 (time, new.n1, new.m2,   # get the discrepancy measures
          new.u2.A.YoY, new.u2.N.YoY, new.u2.A.1, new.u2.N.1, 
          clip.frac.H.YoY, clip.frac.H.1, 
          expit(results$sims.list$logitP), 
          round(results$sims.list$U.W.YoY), 
          cbind(matrix(0,nrow=nrow(results$sims.list$U.H.YoY),ncol=(hatch.after.YoY-min(time)+1)),
                round(results$sims.list$U.H.YoY)), 
          round(results$sims.list$U.W.1), 
          round(results$sims.list$U.H.1), 
          hatch.after.YoY) #don't forget that hatchery fish is 0 until hatch.after
PredictivePosteriorPlot.TSPDE.WHCH2 (discrep)
dev.off()


sink(results.filename, append=TRUE)
# What was the initial seed
cat("\n\n*** Initial Seed for this run ***: ", results$Seed.initial,"\n")

# Global summary of results
cat("\n\n*** Summary of MCMC results *** \n\n")
print(results, digits.summary=3)

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

# Summary of population sizes
cat("\n\n*** Summary of Unmarked Population Size ***\n")
cat("Wild YoY \n")
print(round(results$summary[ grep("Utot.W.YoY", rownames(results$summary)),]))
cat("\n\nWild Age 1 \n")
print(round(results$summary[ grep("Utot.W.1"  , rownames(results$summary)),]))
cat("\n\nHatchery YoY\n")
print(round(results$summary[ grep("Utot.H.YoY", rownames(results$summary)),]))
cat("\n\nHatchery Age 1\n")
print(round(results$summary[ grep("Utot.H.1"   , rownames(results$summary)),]))
cat("\n\nGrand Total\n")
print(round(results$summary[ rownames(results$summary) == "Utot",]))

#browser()
cat("\n\n*** Summary of Quantiles of Run Timing.Wild *** \n")
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
q <- RunTime(time=time[time>hatch.after.YoY], U=results$sims.list$U.H.YoY, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))
cat(    "\n    This is based on the sample weeks provided and the U.H.1[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.H.1, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n")
cat(paste("*** end of fit *** ", date()))

sink()


# add some of the raw data to the bugs object for simplicity in referencing it later
results$data <- list( time=time, n1=n1, m2=m2, 
                      u2.A.YoY=u2.A.YoY, u2.N.YoY=u2.N.YoY, u2.A.1=u2.A.1, u2.N.1=u2.N.1, 
                      clip.frac.H.YoY=clip.frac.H.YoY, clip.frac.H.1=clip.frac.H.1,
                      sampfrac=sampfrac, 
                      hatch.after.YoY=hatch.after.YoY, 
                      bad.m2=bad.m2,
                      bad.u2.A.YoY=bad.u2.A.YoY, bad.u2.N.YoY=bad.u2.N.YoY, 
                      bad.u2.A.1=bad.u2.A.1, bad.u2.N.1=bad.u2.N.1, 
                      logitP.cov=logitP.cov)

return(results)
} # end of function
