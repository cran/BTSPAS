% 2013-12-30 CJS added engine=jags argument to call
% 2013-05-13 CJS fixed all \usage lines more than 90 characters long to avoid R CMD check warnings
% 2012-01-13 CJS commented out the \details, \note, \seealso sections to avoid R CMD check warnings
% 2011-05-03 CJS fixed documentation to match latest release (drop openbugs etc)
%                openbugs/winbugs directory
% 2010-11-25 CJS added bad.u2 argument
% 2010-03-12 CJS Updated argument list on _fit
% 2010-02-06 CJS Added alias to function
% 2009-12-01 CJS Documentation added.
\name{TimeStratPetersenDiagError_fit}
\alias{TimeStratPetersenDiagError_fit}
\alias{TimeStratPetersenDiagError}
\title{Wrapper (*_fit) and function to call the Time Statified Petersen Estimator with Diagonal Entries function.}
\description{Takes the number of marked fish released, the number of recaptures, 
and the number of unmarked fish and 
uses Bayesian methods to fit a fit a spline through the population numbers and a hierarchical model for the
trap efficiencies over time.  The output is written to files and an MCMC object is also created with samples
from the posterior.

Normally, the wrapper (*_fit) function is called which then calls the fitting routine.

Use the \code{\link{TimeStratPetersenNonDiagError_fit}} function for cases where
recaptures take place outside the stratum of release.}

\usage{
TimeStratPetersenDiagError_fit(title="TSDPE", prefix="TSPDE-", 
    time, n1, m2, u2, sampfrac, 
    jump.after = NULL, bad.n1=c(), bad.m2 = c(), bad.u2=c(),
    logitP.cov = rep(1, length(n1)), 
    n.chains = 3, n.iter = 2e+05, n.burnin = 1e+05, n.sims = 2000, 
    tauU.alpha = 1, tauU.beta = 0.05, taueU.alpha = 1, taueU.beta = 0.05, 
    mu_xiP = logit(sum(m2, na.rm = TRUE)/sum(n1, na.rm = TRUE)),  
    tau_xiP = 1/var(logit((m2 + 0.5)/(n1 + 1)), na.rm = TRUE),  
    tauP.alpha = 0.001, tauP.beta = 0.001, run.prob = seq(0, 1, 0.1),  
    debug = FALSE, debug2 = FALSE,
    engine=c("jags","openbugs")[1],
    InitialSeed=ceiling(runif(1,min=0, max=if(engine=="jags"){1000000}else{14})))

TimeStratPetersenDiagError(title, prefix, time, n1, m2, u2, jump.after = NULL,
     logitP.cov, 
     n.chains = 3, n.iter = 2e+05, n.burnin = 1e+05, n.sims = 2000, 
     tauU.alpha = 1, tauU.beta = 0.05, taueU.alpha = 1, taueU.beta = 0.05, 
     mu_xiP = logit(sum(m2, na.rm = TRUE)/sum(n1, na.rm = TRUE)), 
     tau_xiP = 1/var(logit((m2 + 0.5)/(n1 + 1)), na.rm = TRUE), 
     tauP.alpha = 0.001, tauP.beta = 0.001, debug = FALSE, debug2 = FALSE, 
     engine=c("jags","openbugs")[1],
     InitialSeed )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{title}{A character string used for a title on reports and graphs}
  \item{prefix}{A character string used as the prefix for created files. All created graph files are of the form prefix-xxxxx.pdf. }
  \item{time}{A numeric vector of time used to label the strata. For example, this could be julian week for data stratified at a weekly level.  }
  \item{n1}{A numeric vector of the number of marked fish released in each time stratum. }
  \item{m2}{A numeric vector of the number of marked fish from n1 that are recaptured in each time stratum. All recaptures
            take place within the stratum of release.  }
  \item{u2}{A numeric vector of the number of unmarked fish captured in each stratum. These will be expanded by the capture efficiency
            to estimate the population size in each stratum.}
  \item{sampfrac}{A numeric vector with entries between 0 and 1 indicating what fraction of the stratum was sampled. For example, if strata
       are calendar weeks, and sampling occurred only on 3 of the 7 days, then the value of \code{sampfrac} for that stratum would be 3/7.}
  \item{jump.after}{A numeric vector with elements belonging to \code{time}.
      In some cases, the spline fitting the population numbers should be allowed to jump. For example, the population size
      could take a jump when hatchery released fish suddenly arrive at the trap. The jumps occur AFTER the strata listed in this argument.}
  \item{bad.n1}{A numeric vector with elements belonging to \code{time}. In some cases, something goes wrong in the stratum, and the
      number of marked fish released should be ignored.  The values of \code{m2} for this stratum will also be set to NA for these strata.}
  \item{bad.m2}{A numeric vector with elements belonging to \code{time}. In some cases, something goes wrong in the stratum, and the
      number of recovered fish should be ignored. For example, poor handling is suspected to induce handling induced mortality in the
      marked fish and so only very few are recovered. The values of \code{m2} will be set to NA for these strata.}
  \item{bad.u2}{A numeric vector with elements belonging to \code{time}. In some cases, something goes wrong in the stratum, and the
      number of unmarked fish should be ignored. For example, the trap didn't work properly in this stratum.
      The values of \code{u2} will be set to NA for these strata.}
  \item{logitP.cov}{A numeric matrix for covariates to fit the logit(catchability). Default is a single intercept, i.e. all strata
      have the same mean logit(catchability). } \cr \cr

  \item{n.chains}{Number of parallel MCMC chains to fit.}
  \item{n.iter}{Total number of MCMC iterations in each chain.}
  \item{n.burnin}{Number of burn-in iterations.}
  \item{n.sims}{Number of simulated values to keeps for posterior distribution.} \cr


  \item{tauU.alpha}{One of the parameters along with \code{tauU.beta} for the prior for the variance of the random noise for the smoothing spline.}
  \item{tauU.beta}{One of the parameters along with \code{tauU.alpha} for the prior for the variance of the random noise for the smoothing spline.}
  \item{taueU.alpha}{One of the parameters along with \code{taueU.beta} for the prior for the variance of noise around the spline.}
  \item{taueU.beta}{One of the parameters along with \code{taueU.alpha} for the prior for the variance of noise around the spline.}
  \item{mu_xiP}{One of the parameters for the prior for the mean of the logit(catchability) across strata}
  \item{tau_xiP}{One of the parameter for the prior for the mean of the logit(catchability) across strata}
  \item{tauP.alpha}{One of the parameters for the prior for the variance in logit(catchability) among strata}
  \item{tauP.beta}{One of the parameters for the prior for the variance in logit(catchability) among strata}
  \item{run.prob}{Numeric vector indicating percentiles of run timing should be computed.}
  \item{debug}{Logical flag indicating if a debugging run should be made. In the debugging run, the number of samples in the posterior
       is reduced considerably for a quick turn around. }
  \item{debug2}{Logical flag indicated if additional debugging information is produced. Normally the functions will halt at \code{browser()}
       calls to allow the user to peek into the internal variables. Not useful except to package developers.}
  \item{engine}{Which MCMC sampler should be used. Default is "JAGS". Also available is "OpenBugs" or "Winbugs". Case not relevant.}
  \item{InitialSeed}{Numeric value used to initialize the random numbers used in the MCMC iterations.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{ An MCMC object with samples from the posterior distribution. A series of graphs and text file are also created in the
working directory.}
\references{ Refer to the Trinity River Restoration Project report by Schwarz, C.J. et al. (2009)
available at \url{http://www.stat.sfu.ca/~cschwarz/Consulting/Trinity/Phase2}. Please
contact \email{cschwarz@stat.sfu.ca} for more details.
%% ~put references to the literature/web site here ~
}
\author{Bonner, S.J. \email{s.bonner@stat.ubc.ca} and Schwarz, C. J. 
\email{cschwarz@stat.sfu.ca}}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{ 
##---- See the demo files for examples of how to use this package
##
##     demo("demo-TSPDE",     package='BTSPAS', ask=FALSE)  # the simplest usage
##     demo("demo-TSPDE-cov", package='BTSPAS', ask=FALSE)  # including a covariate for logit(P)
##
} % end of examples section

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~models}
\keyword{~smooth}
