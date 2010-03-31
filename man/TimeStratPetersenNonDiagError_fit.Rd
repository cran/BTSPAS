% 2010-03-29 CJs Updated with *NP functions and updated argument list
% 2010-03-12 CJS Updated argument list on _fit
% 2010-03-03 CJS Add ability to fix some logitP to certain values.
% 2010-02-06 CJS Added documentation on called function as well.
% 2009-12-08 CJS Documentation added.
\name{TimeStratPetersenNonDiagError_fit}
\alias{TimeStratPetersenNonDiagError_fit}
\alias{TimeStratPetersenNonDiagError}
\alias{TimeStratPetersenNonDiagErrorNP_fit}
\alias{TimeStratPetersenNonDiagErrorNP}
\title{Wrapper (*_fit) and function to call the Time Statified Petersen Estimator with NON Diagonal Entries function.}
\description{Takes the number of marked fish released, the number of recaptures, 
and the number of unmarked fish and 
uses Bayesian methods to fit a fit a spline through the population numbers and a hierarchical model for the
trap efficiencies over time.  The output is written to files and an MCMC object is also created with samples
from the posterior.

Normally the user makes a call to the *_fit function which then calls the fitting function.

Use the \code{\link{TimeStratPetersenDiagError_fit}} function for cases where
recaptures take place ONLY in the stratum of release, i.e. the diagonal case.

The non-diagonal case fits a log-normal distribution for the travel time. The *NP functions fit
a non-parametric distribution for the travel times.
}
\usage{
TimeStratPetersenNonDiagError_fit(title="TSPNDE", prefix="TSPNDE-", 
    time, n1, m2, u2, sampfrac, jump.after = NULL, 
    bad.n1=c(), bad.m2 = c(), bad.u2=c(),
    logitP.cov = rep(1, length(u2)), 
    logitP.fixed=NULL, logitP.fixed.values=NULL,
    n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
    tauU.alpha = 1, tauU.beta = 0.05, taueU.alpha = 1, taueU.beta = 0.05, 
    mu_xiP = logit(sum(m2, na.rm = TRUE)/sum(n1, na.rm = TRUE)),  
    tau_xiP = .6666,
    tauP.alpha = 0.001, tauP.beta = 0.001, run.prob = seq(0, 1, 0.1),  
    debug = FALSE, debug2 = FALSE, openbugs = TRUE,  
    OPENBUGS.directory = file.path("c:", "Program Files", "OpenBugs"), 
    WINBUGS.directory = file.path("c:", "Program Files", "WinBUGS14"),
    InitialSeed = trunc(runif(1, min = 1, max = 1e+09)))

TimeStratPetersenNonDiagError (
    title, prefix, time, n1, m2, u2,
    jump.after=NULL,
    logitP.cov, logitP.fixed=rep(NA,length(u2)), 
    n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
    tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
    mu_xiP=-2, tau_xiP=.6666, 
    tauP.alpha=.001, tauP.beta=.001, 
    debug=FALSE, debug2=FALSE, openbugs=TRUE, InitialSeed,
    OPENBUGS.directory, WINBUGS.directory )

TimeStratPetersenNonDiagErrorNP_fit( 
    title="TSPNDENP", prefix="TSPNDENP-", 
    time, n1, m2, u2, sampfrac, jump.after=NULL,
    bad.n1=c(), bad.m2=c(), bad.u2=c(),
    logitP.cov=rep(1,length(u2)),
    logitP.fixed=NULL, logitP.fixed.values=NULL, 
    n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
    tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
    mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),
    tau_xiP=.6666,   # need a better initial value for variation in catchability
    tauP.alpha=.001, tauP.beta=.001,
    Delta.max=NULL,tauTT.alpha=.1,tauTT.beta=.1,
    run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
    debug=FALSE, debug2=FALSE, openbugs=TRUE,
    OPENBUGS.directory=file.path("c:","Program Files","OpenBugs"),
    WINBUGS.directory=file.path("c:","Program Files","WinBUGS14"),
    InitialSeed=trunc(runif(1,min=1, max=1000000000))) 

TimeStratPetersenNonDiagErrorNP (
     title, prefix, time, n1, m2, u2,
     jump.after=NULL,
     logitP.cov=rep(1,length(u2)), logitP.fixed=rep(NA,length(u2)),
     n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
     tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
     Delta.max,tauTT.alpha=.1, tauTT.beta=.1,
     mu_xiP=-2, tau_xiP=.6666, 
     tauP.alpha=.001, tauP.beta=.001, 
     debug=FALSE, debug2=FALSE, openbugs=TRUE, InitialSeed,
     OPENBUGS.directory, WINBUGS.directory )
#
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{title}{A character string used for a title on reports and graphs}
  \item{prefix}{A character string used as the prefix for created files. All created graph files are of the form prefix-xxxxx.pdf. }
  \item{time}{A numeric vector of time used to label the strata. For example, this could be julian week for data stratified at a weekly level.  }
  \item{n1}{A numeric vector of the number of marked fish released in each time stratum. }
  \item{m2}{A numeric matrix of the number of fish released in stratum [i] and recovered in [j-1] strata later.
            For example m2[3,5] is the number of marked fish released in stratum 3 and recovered 4 strata later in stratum 7.
            The first column is the number of marked fish recovered in the stratum of release, i.e. 0 strata later.
            Use the \code{\link{TimeStratPetersenDiagError_fit}} function for cases where
            recaptures take place ONLY in the stratum of release, i.e. the diagonal case.}
  \item{u2}{A numeric vector of the number of unmarked fish captured in each stratum. 
            These will be expanded by the capture efficiency
            to estimate the population size in each stratum.}
  \item{sampfrac}{A numeric vector with entries between 0 and 1 indicating 
       what fraction of the stratum was sampled. For example, if strata
       are calendar weeks, and sampling occurred only on 3 of the 7 days, 
        then the value of \code{sampfrac} for that stratum would be 3/7.}
  \item{jump.after}{A numeric vector with elements belonging to \code{time}.
      In some cases, the spline fitting the population numbers should be allowed to jump. 
      For example, the population size
      could take a jump when hatchery released fish suddenly arrive at the trap. 
      The jumps occur AFTER the strata listed in this argument.}
  \item{bad.n1}{A numeric vector with elements belonging to \code{time}. 
      In some cases, something goes wrong in the stratum, and the
      number of marked fish releases should be discarded.
      The values of \code{n1} will be set to NA for these strata.}
  \item{bad.m2}{A numeric vector with elements belonging to \code{time}. 
      In some cases, something goes wrong in the stratum, and the
      number of recovered marked fish should be ignored. For example, poor handling is 
      suspected to induce handling induced mortality in the
      marked fish and so only very few are recovered. 
      The values of \code{m2} in the entire row will be set to NA for these strata.}
  \item{bad.u2}{A numeric vector with elements belonging to \code{time}. 
      In some cases, something goes wrong in the stratum, and the
      number of unmarked fish captred should be ignored.
      The values of \code{u2} in the entire row will be set to NA for these strata.}
  \item{logitP.cov}{A numeric matrix for covariates to fit the logit(catchability). 
      Default is a single intercept, i.e. all strata
      have the same mean logit(catchability). }
  \item{logitP.fixed}{A numeric vector (could be null) of the time strata where the
      logit(P) whould be fixed. Typically, this is used when the capture rates for
      some strata are 0 and logit(P) is set to -10 for these strata. The fixed
      values are given in \code{logitP.fixed.values}}
  \item{logitP.fixed.values}{A numerical vector (could be null) of the fixed values
      for logit(P) at strata given by logitP.fixed. Typically this is used when 
      certain strata have a 0 capture rate and the fixed value is set to -10  which
      on the logit scale gives p[i] essentially 0. Don't specify values such as -50
      because numerical problems could occur in WinBugs/OpenBugs.} \cr \cr


   \bold{The following parameters are usually not set by the user.} \cr

  \item{n.chains}{Number of parallel MCMC chains to fit.}
  \item{n.iter}{Total number of MCMC iterations in each chain.}
  \item{n.burnin}{Number of burn-in iterations.}
  \item{n.sims}{Number of simulated values to keeps for posterior distribution.} \cr

  \item{tauU.alpha}{One of the parameters along with \code{tauU.beta} for the prior for the variance of the random noise for the smoothing spline.}
  \item{tauU.beta}{One of the parameters along with \code{tauU.alpha} for the prior for the variance of the random noise for the smoothing spline.}
  \item{taueU.alpha}{One of the parameters along with \code{taueU.beta} for the prior for the variance of noise around the spline.}
  \item{taueU.beta}{One of the parameters along with \code{taueU.alpha} for the prior for the variance of noise around the spline.}
  \item{Delta.max}{Maximum transition time for marked fish, i.e. all fish assumed to have moved by Delta.max unit of time}
  \item{tauTT.alpha}{One of the parameters along with \code{tauTT.beta} for the prior on 1/var of logit continuation
                    ratio for travel times}
  \item{tauTT.beta}{One of the parameters along with \code{tauTT.alpha} for the prior on 1/var of logit continuation
                    ratio for travel times}
  \item{mu_xiP}{One of the parameters for the prior for the mean of the logit(catchability) across strata}
  \item{tau_xiP}{One of the parameter for the prior for the mean of the logit(catchability) across strata}
  \item{tauP.alpha}{One of the parameters for the prior for the variance in logit(catchability) among strata}
  \item{tauP.beta}{One of the parameters for the prior for the variance in logit(catchability) among strata}
  \item{run.prob}{Numeric vector indicating percentiles of run timing should be computed.}
  \item{debug}{Logical flag indicating if a debugging run should be made. In the debugging run, the number of samples in the posterior
       is reduced considerably for a quick turn around. }
  \item{debug2}{Logical flag indicated if additional debugging information is produced. Normally the functions will halt at \code{browser()}
       calls to allow the user to peek into the internal variables. Not useful except to package developers.}
  \item{openbugs}{Logical flag indicating if OpenBugs (TRUE) or WinBUGS (FALSE) should be used for the MCMC portion.}
  \item{OPENBUGS.directory}{File path indicating where OpenBugs is installed.}
  \item{WINBUGS.directory}{File path indicating where WinBugs is installed.}
  \item{InitialSeed}{Numeric value used to initialize the random numbers used in the MCMC iterations.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ An MCMC object with samples from the posterior distribution. A series of graphs and text file are also created in the
working directory.}
\references{ Refer to the Trinity River Restoration Project report by Schwarz, C.J. et al. (2009)
available at \url{http://www.stat.sfu.ca/~cschwarz/Consulting/Trinity/Phase2}. Please
contact \email{cschwarz@stat.sfu.ca} for more details.
%% ~put references to the literature/web site here ~
}
\author{Bonner, S.J. \email{s.bonner@stat.ubc.ca} and Schwarz, C. J. 
\email{cschwarz@stat.sfu.ca}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{ 
##---- See the demo files for examples of how to use this package
##
##     demo("demo-TSPNDE",     package='BTSPAS', ask=FALSE)  # the simplest usage
##
} % end of examples section

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~models}
\keyword{~smooth}
