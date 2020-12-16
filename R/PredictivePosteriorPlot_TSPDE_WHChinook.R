# 2020-12-15 CJS Fixed problem with u2==NA in the plots
# 2015-06-10 CJS Updated to produce ggplot2 plots
# 2014-09-21 CJS change all Inf to NA
# 2012-01-22 CJS made X/Y axis limits the same so p-value prints properly
# 2011-06-13 CJS returned bayesian p-values

PredictivePosteriorPlot.TSPDE.WHCH <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the discrep matrix has the following columns
#     (1-2)  o,s Freeman-Tukey measures for m2
#     (3-4)  o,s Freeman-Tukey measures for u2.A
#     (5-6)  o,s Freeman-Tukey measures for u2.N
#     (7-8)  o,s Freeman-Tukey for m2+u2.A+u2.N

# Change any Inf to NA
temp <- is.infinite(discrep) & !is.na(discrep)
if(sum(temp, na.rm=TRUE)>0){cat(sum(temp, na.rm=TRUE), " infinite discrepancy measures set to NA\n")}
discrep[ temp ] <- NA

#browser()
discrep.long <- data.table::melt( data.table::as.data.table(discrep), 
                                    measure.vars=list(seq(1,ncol(discrep),2), seq(2,ncol(discrep),2)),
                                    value.name=c("Observed","Simulated"),
                                    variable.name="Statistic",
                                    variable.factor=FALSE)

titles <- data.frame(Statistic=as.character(1:(ncol(discrep/2))), Title=c( 
            "Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A", 
            "Freeman-Tukey for u2.N", 
            "Total Freeman-Tukey"), stringsAsFactors=FALSE)

discrep.long <- merge(discrep.long, titles)

# compute the bayesian p-values
p_values <-plyr::ddply(discrep.long, c("Statistic","Title"), function(x){
       p.value=mean(x$Observed < x$Simulated, na.rm=TRUE)
       data.frame(p.value=p.value)
})
p_values$label = paste("Bayesian GOF P:",formatC(p_values$p.value, digits=2, format="f"))
  
gof.plot <-ggplot(data=discrep.long, aes_(x=~Simulated, y=~Observed))+
       geom_point()+
       geom_abline(intercept=0, slope=1)+
       geom_text(data=p_values, x=Inf,y=-Inf, hjust=1.05, vjust=-0.2, label=p_values$label)+
       facet_wrap(~Title, ncol=2, nrow=3, scales="free")
 
gof <- list(bp.plot=gof.plot,  bp.values=data.frame(test.names=titles, p.value=p_values, stringsAsFactors=FALSE))
gof
}
