
#' Create an acf plot of a parameter
#'
#' @param mcmc.sample Data frame with 2 columns. parm, and sample. A separate ACF plot is generated for each parameter using facet_wrap.
#' @param ncol Number of columns in the plot. 
#'
#' @return acf plot(s) as an ggplot2 object 
#' @keywords internal
#' @import plyr ggplot2
#'
#'
plot_acf <- function(mcmc.sample, ncol=2){
  acf.parm <- plyr::ddply(mcmc.sample, "parm", function(x){
     acf.list <- acf(x$sample, plot=FALSE)
     data.frame(lag=acf.list$lag, acf=acf.list$acf, stringsAsFactors=FALSE)
  })
  acfplot <- ggplot(data=acf.parm, aes(x = lag, y = acf)) +
     ggtitle("Autocorrelation")+
     geom_hline(aes(yintercept = 0)) +
     geom_segment(aes(xend = lag, yend = 0))+
     facet_wrap(~parm, ncol=ncol)
  acfplot
}