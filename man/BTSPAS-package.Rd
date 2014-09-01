\name{BTSPAS-package}
\alias{BTSPAS-package}
\alias{BTSPAS}
\docType{package}
\title{
Bayesian Time-Stratified Population Analysis
}
\description{
BTSPAS provides advanced Bayesian methods to estimate
abundance and run-timing from temporally-stratified
Petersen mark-recapture experiments. Methods include
hierarchical modelling of the capture probabilities
and spline smoothing of the daily run size. \cr
\cr
Full details are available at
\url{http://www.stat.sfu.ca/~cschwarz/Consulting/Trinity/Phase2}

}
\details{
\tabular{ll}{
Package: \tab BTSPAS\cr
Version: \tab 2014.09\cr
Date: \tab 2014-09-01\cr
Depends: \tab R2WinBUGS, coda, splines, actuar, plyr, ggplot2\cr
License: \tab GPL (>=2)\cr
URL: \tab http://www.r-project.org\cr
Repository: \tab R-Forge\cr
Repository/R-Forge/Project: \tab btspas\cr
Repository/R-Forge/Revision: \tab 8\cr
Date/Publication: \tab 2014-09-02 02:49:34\cr
Packaged: \tab 2009-12-02 21:03:35 UTC; rforge\cr
Built: \tab R 2.10.0; ; 2009-12-03 23:01:26 UTC; windows\cr
}

Index: \cr

\link{TimeStratPetersenDiagError_fit}	Wrapper to call the Time Statified Petersen Estimator with Diagonal Entries function.\cr
\link{TimeStratPetersenDiagErrorWHChinook_fit}	Wrapper to call the Time Statified Petersen Estimator with Diagonal Entries to separate wild and hatchery Chinook function.\cr
\link{TimeStratPetersenDiagErrorWHChinook2_fit}	Wrapper to call the Time Statified Petersen Estimator with Diagonal Entries to separate wild and hatchery Chinook function when age 1 fish residualize overwinter.\cr
\link{TimeStratPetersenDiagErrorWHSteel_fit}	Wrapper to call the Time Statified Petersen Estimator with Diagonal Entries to separate wild and hatchery Steelhead function. \cr
\link{TimeStratPetersenNonDiagError_fit}	Wrapper to call the Time Statified Petersen Estimator with NON Diagonal Entries. \cr
}
\author{
Carl J Schwarz \email{cschwarz@stat.sfu.ca} and
Simon J Bonner \email{s.bonner@stat.ubc.ca}

Maintainer: Simon J Bonner \email{s.bonner@stat.ubc.ca}
}
\keyword{ package }
