% 2012-10-01 CJS Comment out the \details, \note, \references, \seealso sections to avoid
%                warnings when a build is done.
% 2011-02-16 CJS Add documentation about computing TimeToTargetRunSize
\name{TimeToTargetRunSize}
\alias{TimeToTargetRunSize}
\title{Computes and plots posterior distribution of time to get target run size. For example, the time
 to reach a cumulative run of 10,000 fish.}
\description{Takes a sim.list object from the MCMC runs, computes the posterior distribution
of the time to the target runsize, plots the posterior}

\usage{
TimeToTargetRunSize( U, time, targetU, file_prefix, ci_prob=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{Elements of sim.list from MCMC object for U - the estimate runsize in each stratum}
  \item{time}{Vector of stratum time indices}
  \item{targetU}{The targeted cumulative run size. E.g. 10,000}
  \item{file_prefix}{Character string giving prefix for plot. A plot will be produced 
        of the posterior in the filename paste(file_prefix,"-target.pdf",sep="")).}
  \item{ci_prob}{What size of credible interval should be computed?}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A list with a sample of the posterior (index), quantiles (quantiles), mean (mean), median(median), and 
       standard deviation (sd), and target value (targetU)}
%\references{ 
%}
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
\dontrun{
# Compute the posterior of time to reach 10,000 fish. Results contains the MCMC object
# 
results$TimeToTargetRunSize <- TimeToTargetRunSize( 
        U=results$sims.list$U,
        time=results$data$time,
        targetU=10000,
        file_prefix = 'Time10000')

} % end of dontrun
} % end of examples section

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~models}
\keyword{~plots}
