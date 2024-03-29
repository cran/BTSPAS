#' Generate Predictive Posterior Plots (Bayesian p-values) for number of models.
#' 
#' This is an internal function, not normally of use to users.


#' @aliases PredictivePosterior_TSPDE_WHChinook 
#' @aliases PredictivePosterior.TSPDE.WHCH2 
#' @aliases PredictivePosterior.TSPDE.WHSteel
#' @aliases PredictivePosterior.TSPNDE 
#' @aliases PredictivePosterior.TSPNDENP
#' @aliases PredictivePosteriorPlot.TSPDE 
#' @aliases PredictivePosteriorPlot.TSPDE.WHCH 
#' @aliases PredictivePosteriorPlot.TSPDE.WHCH2 
#' @aliases PredictivePosteriorPlot.TSPDE.WHSteel 
#' @aliases PredictivePosteriorPlot.TSPNDE 
#' @importFrom stats sd dbinom rbinom
#' @import plyr

#' @keywords internal

# 2018-12-15 CJS added in logitP fixed entries.

PredictivePosterior.TSPDE <- function (n1, m2, u2, 
                                       logitP.fixed, p, 
                                       U) {
#  Generate Predictive Posterior Plot (Bayesian p-value) given the data
#  for a TimeStratified Petersen with Diagonal Elements and error
#    n1, m2, u2  = vectors of input data
#    p, U        = matrix of values (rows=number of posterior samples, columns=strata)
#                  These are returned from the call to JAGS
#
#cat("Call to PredictivePosterior\n")
#browser()
discrep <- matrix(0, nrow=0, ncol=12)
select.m2 <- !is.na(m2)
select.u2 <- !is.na(u2)

for(i in 1:nrow(p)){
   # generate sample data
   gen.m2 <- stats::rbinom(ncol(p), n1, p[i,])
   gen.u2 <- stats::rbinom(ncol(p), U[i,], p[i,])
   # compute a discrepancy measure
   # Observed vs expected values for recaptures of marked fish
     temp <- sqrt(m2) - sqrt(n1*p[i,])
     d1.m2.o <- sum( temp[select.m2]^2, na.rm=TRUE)
     temp <- sqrt(gen.m2) - sqrt(n1*p[i,])
     d1.m2.s <- sum( temp[select.m2]^2, na.rm=TRUE)
     
   # Observed vs expected values for captures of unmarked fish
     temp <- sqrt(u2) - sqrt(U[i,]*p[i,])
     d1.u2.o <- sum( temp[select.u2]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2) - sqrt(U[i,]*p[i,])
     d1.u2.s <- sum( temp[select.u2]^2, na.rm=TRUE)
     
   # Deviance (-2*log-likelihood )
     temp <- stats::dbinom(m2,     n1, p[i,], log=TRUE)
     d2.m2.o <- -2*sum(temp[select.m2])
     temp <- stats::dbinom(gen.m2, n1, p[i,], log=TRUE)
     d2.m2.s <- -2*sum(temp[select.m2])
     temp <- stats::dbinom(u2,     U[i,], p[i,], log=TRUE)
     d2.u2.o <- -2*sum(temp[select.u2])
     temp <- stats::dbinom(gen.u2, U[i,], p[i,], log=TRUE)
     d2.u2.s <- -2*sum(temp[select.u2])
     
   # combined discrepancy measures
     d1.o <- d1.m2.o + d1.u2.o
     d1.s <- d1.m2.s + d1.u2.s
     d2.o <- d2.m2.o + d2.u2.o
     d2.s <- d2.m2.s + d2.u2.s
   # update the array
     discrep <- rbind(discrep, 
              c(d1.m2.o, d1.m2.s, d2.m2.o, d2.m2.s, 
                d1.u2.o, d1.u2.s, d2.u2.o, d2.u2.s, 
                d1.o   , d1.s,    d2.o,    d2.s)) 
}
#browser()
discrep
}
