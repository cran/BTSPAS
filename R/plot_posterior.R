#' Create posterior plot of a parameter with credible interval limits shown as vertical lines
#'
#' @param mcmc.sample Data frame with 2 columns. parm, and sample. A separate posterior plot is generated for each parameter.
#' @param alpha Used to determine credible interval.
#' @param ncol Number of columns in the plot (default=1).
#' @return Posterior plot(s) as an ggplot2 object 
#' @keywords internal
#' @import plyr ggplot2 scales
#'
#'

# R CMD check gets upset with ggplot because it thinks there are no visible bindings for quant etc.
# see https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# for how to get around this
plot_posterior <- function(mcmc.sample, alpha=0.05, ncol=1){
  qparm <- plyr::ddply(mcmc.sample, "parm", function(x){
     quants<- quantile(x$sample, probs=c(alpha/2, 1-alpha/2))
     data.frame(quants=quants, stringsAsFactors=FALSE)
  })
  post_stat <- plyr::ddply(mcmc.sample, "parm", plyr::summarize,
                    post.mean=signif(mean(sample),5),
                    post.sd  =signif(sd  (sample),5))
  
  postplot <- ggplot(data=mcmc.sample, aes_(x = ~sample, y =~ ..density..)) +
     ggtitle(paste("Posterior plots with ", formatC(100*(1-alpha), format="f", digits=0),"% credible intervals",sep=""))+
     geom_vline(data=qparm, aes_(xintercept =~ quants)) +
     geom_density()+
     facet_wrap(~parm, ncol=ncol,scales="free")+xlab("Value of parameter")+
     geom_text(data=post_stat, aes_(x=Inf, y=Inf, label=~paste("Posterior mean: ",post.mean,"\nPosterior sd  : ",post.sd,sep=""),
                                   vjust=1, hjust=1))+
     scale_x_continuous(labels=scales::comma)
  postplot
}

