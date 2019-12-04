#' @rdname PredictivePosterior.TSPDE
#' @import ggplot2 ggforce plyr
#' @importFrom data.table setDT

# We need to use the importFrom here to get the data.table package in the namespace but to avoid warning messages
# from R CMD check about melt() and dcast() being overwritten from the reshape2 package.

# 2018-12-03 CJS convert to using facets rather that marrangeGrob as the latter had problems
# 2015-06-10 CJS Converted to ggplot()
# 2014-09-01 CJS Change Inf to NA 
# 2012-01-22 CJS Made X/Y axis limits the same so that Bayesian p-value prints properly
# 2011-06-13 CJS returned p-values
# 2010-03-29 CJS First creation of routine

PredictivePosteriorPlot.TSPDE.WHCH2 <- function( discrep, ncol=2, nrow=2 ) {
  #   Given the discrepancy measures, creates a set of panel plots.
  #   It is assumed that the discrepancy measure has 16 columns for the Bayesian p-value plot
  #     ( 1- 2)  o,s Freeman-Tukey measures for m2
  #     ( 3- 4)  o,s Freeman-Tukey measures for u2.A.YoY
  #     ( 5- 6)  o,s Freeman-Tukey measures for u2.N.YoY
  #     ( 7- 8)  o,s Freeman-Tukey measures for u2.A.1
  #     ( 9-10)  o,s Freeman-Tukey measures for u2.N.1
  #     (11-12)  o,s Freeman-Tukey for u2.A.YoY+u2.N.YoY
  #     (13-14)  o,s Freeman-Tukey for u2.A.1  +u2.N.1
  #     (15-16)  o,s Freeman-Tukey for all data (m2, YoY and Age 1)`

  # Change any Inf to NA
  temp <- discrep == Inf | discrep == -Inf
  if(sum(temp)>0){cat(sum(temp), " infinite discrepancy measures set to NA\n")}
  discrep[ temp ] <- NA
  

  # Convert from wide to long format and add labels
  discrep.long <- data.table::melt( data.table::as.data.table(discrep), 
                                    measure.vars=list(seq(1,ncol(discrep),2), seq(2,ncol(discrep),2)),
                                    value.name=c("Observed","Simulated"),
                                    variable.name="Statistic",
                                    variable.factor=FALSE)

  titles <- data.frame(Statistic=as.character(1:8), Title=c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A.YoY", 
            "Freeman-Tukey for u2.N.YoY", 
            "Freeman-Tukey for u2.A.1", 
            "Freeman-Tukey for u2.N.1", 
            "Freeman-Tukey for YoY",
            "Freeman-Tukey for Age 1",
            "Total Freeman-Tukey"), stringsAsFactors=FALSE)
  discrep.long <- merge(discrep.long, titles)
 
  # compute the bayesian p-values         
  p_values <-plyr::ddply(discrep.long, c("Statistic","Title"), function(x){
       p.value=mean(x$Observed < x$Simulated)
       data.frame(p.value)
  })
  p_values$label = paste("Bayesian GOF P:",formatC(p_values$p.value, digits=2, format="f"))
  #browser()
  gof <- plyr::llply(1:2, function (page){
    ggplot(data=discrep.long, aes_(x=~Simulated, y=~Observed))+
       geom_point()+
       geom_abline(intercept=0, slope=1)+
       geom_text(data=p_values, x=Inf,y=-Inf, hjust=1, vjust=0, label=p_values$label)+
       facet_wrap_paginate(~Title, ncol=ncol, nrow=nrow, page=page, scales="free")
  })
  gof
}  # end of function
