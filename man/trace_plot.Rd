% 2012-01-10 CJS commented out \details, \note, \references, \seealso to avoid
%                warnings during build
% 2010-08-06 CJS Add documentation about trace_plot
\name{trace_plot}
\alias{trace_plot}
\title{Creates trace plots of specified parameters showing the multiple chains and the value of Rhat}
\description{Takes the MCMC object returned from a split and produces trace_plots for the
listed parameters. It shows a separate line on the plot for each chain and also shows the
value of Rhat}

\usage{
trace_plot(title=" ", results=NULL, parms_to_plot=NULL, panels=c(1,1),
   mai=if(prod(panels)>1){c(.4,.4,.4,.4)} else {c(1.02,0.82,0.82,0.42)},
   cex=if(prod(panels)>1){.5} else {1}
   ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{title}{A character string used for a title on reports and graphs}
  \item{results}{The MCMC object containing the results from the call to WinBugs/OpenBugs}
  \item{parms_to_plot}{A character vector of names of parameters to plot. These must match 
        exactly to the parameter names used in the simulation.}
  \item{panels}{How many plots to put on a page. It used split.screen to format the plots.}
  \item{mai}{Margins (inches) from the 4 borders. See par().}
  \item{cex}{Character expansion factor. See par() for more details.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{ Nothing returned. Creates the plots}
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
# Create trace plots of the logitP parameters
# 
# extract the names of the variables from the MCMC object
varnames <- names(results$sims.array[1,1,]) 
# get the parms that start with logitP
parm.names <- varnames[grep("^logitP", varnames)] 
# create a pdf file of the plots
pdf(file="trace-logitP.pdf",sep=""))
trace_plot(title=title, results=results, 
    parms_to_plot=parm.names, panels=c(3,2))
dev.off()

# Or if you want an interactive display
par(ask=TRUE)
trace_plot(title=title, results=results, 
    parms_to_plot=parm.names, panels=c(2,1))
} % end of dontrun
} % end of examples section

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~models}
\keyword{~trace plot}
