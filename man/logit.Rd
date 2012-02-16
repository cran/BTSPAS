% 2012-01-10 CJS commented out the \details to avoid warnings on build
% 2009-12-01 CJS
\name{logit}
\alias{logit}
\alias{expit}
\title{Logit and anti-logit function.}
\description{ Compute the logit or anti-logit.}
\usage{ 
logit(p) 
expit(theta)
}
\arguments{
  \item{p}{probability between 0 and 1.}
  \item{theta}{logit between -infinity and +infinity}
}
%\details{ }
 \value{Computed logit or anti-logit}
\author{C.J.Schwarz \email{cschwarz@stat.sfu.ca}}
\examples{
##---- compute the logit and its inverse
logitp <- logit(.3)
p <- expit(-.84)
} % end of examples

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~misc}
