#' Logit and anti-logit function.
#' 
#' Compute the logit or anti-logit.
#' 
#' 
#' @aliases logit expit
#' @param p probability between 0 and 1.
#' @param theta logit between -infinity and +infinity
#' @return Computed logit or anti-logit
#' @author C.J.Schwarz \email{cschwarz@@stat.sfu.ca}
#' @keywords ~misc
#' @examples
#' 
#' ##---- compute the logit and its inverse
#' logitp <- logit(.3)
#' p <- expit(-.84)
#' 
#' @export logit expit
#' 
logit <- function(p){
#   logit of p
    log(p/(1-p))
}



#' @rdname logit 
expit <- function(theta){
# anti logit function
   1/(1+exp(-theta))
}
