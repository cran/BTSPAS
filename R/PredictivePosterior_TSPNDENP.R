#' @rdname PredictivePosterior.TSPDE
#' @importFrom stats sd dbinom dmultinom rbinom rmultinom 

# 2021-10-05 CJS Fixed computation of discrepancy measure when some m2.expanded are missing
# 2020-12-15 CJS Fixed computation of discrepancy measures when u2 is missing
# 2018-11-30 CJS Changed defintion of m2.expanded propogates down here
# 2018-11-27 CJS Removed openbugs stuff
# 2014-09-01 CJS Fixed bug when logitP.fixed is fixed in first position
# 2014-09-01
#    There was also a subtle bug in dealing with the multinomial distribution where the length of p 
#    (that had to be padded to deal with an OpenBugs problem
#    had to have the indicies explicitly stated.


PredictivePosterior.TSPNDENP <- function (n1,
                                          m2.expanded,
                                          u2,
                                          logitP.fixed,
                                          p,
                                          U,
                                          Theta,
                                          Delta.max) {
#  Generate Predictive Posterior Plot (Bayesian p-value) given the data
#  for a TimeStratified Petersen with Non-Diagonal entries and a non-parametric movement model
#    n1          - vector of number of releases
#    m2.expanded - matrix of number of marks recovered from each value of n1
#    u2          - vector of recovieres of unmarked fish
#    logitP.fixed - which values of p(i) are fixed
#    Delta.max    - maximum strata involved in movement from diagonal
#    U, Theta  = matrix of values (rows=number of posterior samples, columns=strata)
#                  These are returned from the call to  JAGS
#


s <- length(n1)
t <- length(u2)

select.u2 <- !is.na(u2) # which terms involving u2 to use?

## Transform saved iterations for theta from vectors to full movement matrices
Theta.bkp <- Theta

Theta <- lapply(1:nrow(Theta.bkp),function(k){
  M <- Theta.bkp[k,,]
  tmp <- matrix(0,nrow=s,ncol=t)
  for(i in 1:length(n1)){
    tmp[i,i:min(t,i+Delta.max)] <- M[i,1:min(t-i+1,Delta.max+1)]
  }
  tmp
})

## Simulate data for each iteration
#browser()
simData <- lapply(1:nrow(p),function(k) simTSPNDE(n1,U[k,],p[k,],Theta[[k]]))

#browser()
## Compute discrepancy measures
#browser()
discrep <- t(sapply(1:nrow(p),function(k){

  ## 1) Observed vs expected values for recaptures of marked fish
  ## a) Observed data
  temp1.o <- sqrt(m2.expanded[,1:t]) - sqrt(n1 * t(t(Theta[[k]]) * p[k,1:t]))
  d1.m2.o <- sum(temp1.o^2,na.rm=TRUE)

  ## b) Simulate data
  temp1.s <- sqrt(simData[[k]]$m2[,1:t]) - sqrt(n1 * t(t(Theta[[k]]) * p[k,1:t]))
  d1.m2.s <- sum(temp1.s^2,na.rm=TRUE)

  ## 2) Observed vs expected values for captures of unmarked fish
  ## a) Observed data
  temp2.o <- sqrt(u2) - sqrt(U[k,] * p[k,1:t])
  d1.u2.o <- sum(temp2.o[select.u2]^2,na.rm=TRUE)

  ## b) Simulate data
  temp2.s <- sqrt(simData[[k]]$u2) - sqrt(U[k,] * p[k,1:t])
  d1.u2.s <- sum(temp2.s[select.u2]^2,na.rm=TRUE)

  ## 3) Deviance (-2*log-likelihood)
  ## a) Observed data
  #browser()
  d2.m2.o <- -2 * sum(sapply(1:length(n1),function(i){
    cellProbs <- Theta[[k]][i,] * p[k,1:t]  # 2014-09-01 need to ignore extra p's at end which were needed for OPENbugs quirk
    cellProbs <- c(cellProbs,1-sum(cellProbs))
    res <- 0
    if(!is.na(sum(m2.expanded[i,1:t]))){  # data is present, so compute the deviance; otherwise nothing
       res<- stats::dmultinom(c(m2.expanded[i,1:t], n1[i]-sum(m2.expanded[i,1:t])),n1[i],cellProbs,log=TRUE)
    }
    res
  }))

  d2.u2.o <- -2 * sum(stats::dbinom(u2[select.u2],U[k,select.u2],p[k,select.u2],log=TRUE), na.rm=TRUE)

  d2.o <- d2.m2.o + d2.u2.o

  ## b) Simulated data
  d2.m2.s <- -2 * sum(sapply(1:length(n1),function(i){
    cellProbs <- Theta[[k]][i,] * p[k,1:t]  # 2014-09-01 ditto to previous fix
    cellProbs <- c(cellProbs,1-sum(cellProbs))
    res <- 0
    if(!is.na(sum(m2.expanded[i,1:t]))){  # data is present, so compute the deviance on sim data; otherwise nothing
       res<- stats::dmultinom(simData[[k]]$m2[i,],n1[i],cellProbs,log=TRUE)
    }
    res
  }))

  d2.u2.s <- -2 * sum(stats::dbinom(simData[[k]]$u2[select.u2],U[k,select.u2],p[k,select.u2],log=TRUE), na.rm=TRUE)

  d2.s <- d2.m2.s + d2.u2.s

  c(d1.m2.o, d1.m2.s, d2.m2.o, d2.m2.s,
    d1.u2.o, d1.u2.s, d2.u2.o, d2.u2.s,
    d1.m2.o+d1.u2.o, d1.m2.s+d1.u2.s, d2.o, d2.s)
}))

discrep
}

simTSPNDE <- function(n1,U,p,Theta){
  ## Simulate data from the TSPNDE model conditional on values of n and U.

  s <- length(n1)
  t <- length(U)

  ## 1) Simulate matrix of recoveries
  m2 <- t(sapply(1:length(n1),function(i){
    cellProbs <- Theta[i,] * p[1:t]
    cellProbs <- c(cellProbs,1-sum(cellProbs))
    if( any(cellProbs < 0)){browser()}
    stats::rmultinom(1,n1[i],cellProbs)[1:t]
  }))

  ## 2) Add number of individuals not recovered to last column of m2
  m2 <- cbind(m2,n1-apply(m2,1,sum))

  ## 3) Simulate captures of unmarked fish
  u2 <- stats::rbinom(t,U,p)

  return(list(m2=m2,u2=u2))
}
