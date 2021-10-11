#' Compute percentiles of the run timing distribution.
#' 
#' Take the posterior sample of U[1,...nstrata] and compute the percentiles of the run
#' timing. 
#' This uses the quantile() function from the "actuar" package which is designed to compute
#' quantiles of grouped data.
#' It is assumed that there are no fish in the system prior to the first point
#            in time, and after the last point in time

#' @template time 
#' @param U matrix of posterior samples. Each row is a sample from the posterior.
#          Columns correspond to U[1]...U[nstrata]
#' @param prob Quantiles of the run timing to estimate. 

#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the working directory.
#' This information is now added to the fit object as well and so it is unlikely
#' that you will use this function.
#' @template author 
#' @template references 
#' @export RunTime
#' @import plyr
#' @importFrom actuar grouped.data
#' @importFrom stats quantile

# 2018-12-14 CJS converted from a for() loop to adply()

RunTime <- function(time, U, prob=seq(0,1,.1)) {
  timing <- c(min(time):(1+max(time)))
  q.U <- plyr::adply(U, 1, function(U.sample, timing){
       quant <- stats::quantile(actuar::grouped.data(Group=timing, Frequency=U.sample), prob=prob)
       quant
  }, timing=timing, .id=NULL)
  q.U
}
