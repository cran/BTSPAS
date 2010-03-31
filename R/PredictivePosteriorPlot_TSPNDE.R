PredictivePosteriorPlot.TSPNDE <- function( discrep  ) {
#   Given the discrepancy measures, creates a set of panel plots.
#   It is assumed that the bp has 12 columns
#     (1-2)  o,s Freeman-Tukey measures for m2
#     (3-4)  o,s Deviance for m2
#     (5-6)  o,s Freeman-Tukey measures for u2
#     (7-8)  o,s Deviance for u2
#      9-10  o,s Freeman-Tukey for m2+u2
#     11-12  o,s Deviance for m2+uw

split.screen(figs=c(3,2))  # 3 rows and 2 columns
titles <- c("Freeman-Tukey for m2", 
            "Deviance for m2",
            "Freeman-Tukey for u2", 
            "Deviance for u2",
            "Total Freeman-Tukey",
            "Total deviance")

for(i in 1:6){
  screen(i)
  par(cex=.5)
  par(mai=c(.40,.40,.40,.40)) # margins of plot relative to plotting position
  plot(discrep[,(2*i-1):(2*i)], ylab="Observed", xlab="Simulated", main=titles[i], cex.main=1.5)
  abline(a=0, b=1)
  p.value <- sum(discrep[,2*i-1]<discrep[,2*i])/nrow(discrep)
  # locate where to plot
  x.loc <- min(discrep[,2*i-1])+.10*(max(discrep[,2*i-1])-min(discrep[,2*i-1]))
  y.loc <- max(discrep[,2*i  ])-.10*(max(discrep[,2*i  ])-min(discrep[,2*i  ]))
  text(x.loc, y.loc, labels=paste("Bayesian GOF P:",formatC(p.value, digits=2, format="f")), cex=1.5, adj=c(0,0))  
}
close.screen(all=TRUE)     # exit from plots  
}