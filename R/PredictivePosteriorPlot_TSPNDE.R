#' @rdname PredictivePosterior.TSPDE

# 2019-02-13 CJS convert from arrangeGrob to facet_wrap because arrangeGrob dosn't work in paralle.
# 2015-06-10 CJS convert to ggplot()
# 2014-09-01 CJS dealing with Inf and -Inf in the discrepancy measures
# 2012-01-22 CJS made X/Y axis limits the same so p-value prints properly
# 2011-06-13 CJS returned bayesian p-values

PredictivePosteriorPlot.TSPNDE <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the bp has 12 columns
#     (1-2)  o,s Freeman-Tukey measures for m2
#     (3-4)  o,s Deviance for m2
#     (5-6)  o,s Freeman-Tukey measures for u2
#     (7-8)  o,s Deviance for u2
#      9-10  o,s Freeman-Tukey for m2+u2
#     11-12  o,s Deviance for m2+u2

# Change any Inf to NA
temp <- discrep == Inf | discrep == -Inf
if(sum(temp)>0){cat(sum(temp), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

discrep.long <- data.table::melt( data.table::as.data.table(discrep), 
                                    measure.vars=list(seq(1,ncol(discrep),2), seq(2,ncol(discrep),2)),
                                    value.name=c("Observed","Simulated"),
                                    variable.name="Statistic",
                                    variable.factor=FALSE)

titles <- data.frame(Statistic=as.character(1:6), Title=c( 
            "Freeman-Tukey for m2", 
            "Deviance for m2", 
            "Freeman-Tukey for u2", 
            "Deviance for u2", 
            "Total Freeman-Tukey",
            "Total Deviance"), stringsAsFactors=FALSE)
discrep.long <- merge(discrep.long, titles)

# compute the bayesian p-values
p_values <-plyr::ddply(discrep.long, c("Statistic","Title"), function(x){
       p.value=mean(x$Observed < x$Simulated)
       data.frame(p.value=p.value)
})
p_values$label = paste("Bayesian GOF P:",formatC(p_values$p.value, digits=2, format="f"))
  
gof.plot <-ggplot(data=discrep.long, aes_(x=~Simulated, y=~Observed))+
       geom_point()+
       geom_abline(intercept=0, slope=1)+
       geom_text(data=p_values, x=Inf,y=-Inf, hjust=1, vjust=0, label=p_values$label)+
       facet_wrap(~Title, ncol=2, nrow=3, scales="free")
 
gof <- list(bp.plot=gof.plot,  bp.values=data.frame(test.names=titles, p.value=p_values, stringsAsFactors=FALSE))
#browser()
gof
}
