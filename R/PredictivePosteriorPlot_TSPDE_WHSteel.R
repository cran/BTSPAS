#' @rdname PredictivePosterior.TSPDE


# 2015-06-10 CJS converted to ggplot()
# 2014-09-01 CJS Inf in discrep set to NA
# 2012-01-22 CJS made X/Y axis limits the same so p-value prints properly
# 2011-06-13 CJS returned bayesian p-values

PredictivePosteriorPlot.TSPDE.WHSteel <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the bp has 12 columns
#     (1-2)  o,s Freeman-Tukey measures for m2
#     (3-4)  o,s Freeman-Tukey measures for u2.W.YoY
#     (5-6)  o,s Freeman-Tukey measures for u2.W.1
#     (7-8)  o,s Freeman-Tukey measures for U2.H.1
#     (9-10) o,s combined Freeman-Tukey 

# Change any Inf to NA
temp <- discrep == Inf | discrep == -Inf
if(sum(temp)>0){cat(sum(temp), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

#browser()
titles <- c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.W.YoY",
            "Freeman-Tukey for u2.W.1",
            "Freeman-Tukey for u2.H.1", 
            "Total Freeman-Tukey")

saved_p_values <- rep(NA, length(titles))

discrep.df <- data.frame(discrep)
plot.list<- llply(1:5, function(i){
  p.value <- sum(discrep[,2*i-1]<discrep[,2*i],na.rm=TRUE)/nrow(discrep)
  saved_p_values[i] <<- p.value
  
  bp.plot <- ggplot(data=discrep.df, aes_string(x=colnames(discrep.df)[2*i], y=colnames(discrep.df)[2*i-1]))+
    ggtitle(titles[i])+
    geom_point()+
    xlab("Simulated")+ylab("Observed")+
    geom_abline(intercept=0, slope=1)+
    annotate("text", x=Inf,y=-Inf, hjust=1, vjust=0,
               label=paste("Bayesian GOF P:",formatC(p.value, digits=2, format="f")))
  bp.plot 
})

bigplot <- do.call(arrangeGrob, c(plot.list, list(ncol=2)))

gof <- list(bp.plot=bigplot,  bp.values=data.frame(test.names=titles, p.value=saved_p_values, stringsAsFactors=FALSE))
gof
}
