# 2010-03-29 CJS First creation of routine

PredictivePosteriorPlot.TSPDE.WHCH2 <- function( discrep  ) {
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

#browser()
titles <- c("Freeman-Tukey for m2", 
            "Freeman-Tukey for u2.A.YoY", 
            "Freeman-Tukey for u2.N.YoY", 
            "Freeman-Tukey for u2.A.1", 
            "Freeman-Tukey for u2.N.1", 
            "Freeman-Tukey for YoY",
            "Freeman-Tukey for Age 1",
            "Total Freeman-Tukey")

for(page in 1:2){
   split.screen(figs=c(2,2))  # 4 rows and 2 columns
   for(i in 1:4){
     screen(i)
     par(cex=.5)
     par(mai=c(.40,.40,.40,.40)) # margins of plot relative to plotting position
     plot(discrep[,((page-1)*4+(2*i-1)):((page-1)*4+(2*i))], ylab="Observed", xlab="Simulated", 
         main=titles[(page-1)*4+i], cex.main=1.5)
     abline(a=0, b=1)
     p.value <- sum(discrep[,(page-1)*4+2*i-1]<discrep[,(page-1)*4+2*i])/nrow(discrep)
     # locate where to plot
     x.loc <- min(discrep[,(page-1)*4+2*i-1])+.10*(max(discrep[,(page-1)*4+2*i-1])-min(discrep[,(page-1)*4+2*i-1]))
     y.loc <- max(discrep[,(page-1)*4+2*i  ])-.10*(max(discrep[,(page-1)*4+2*i  ])-min(discrep[,(page-1)*4+2*i  ]))
     text(x.loc, y.loc, labels=paste("Bayesian GOF P:",formatC(p.value, digits=2, format="f")), cex=1.5, adj=c(0,0))  
   }
   close.screen(all=TRUE)     # exit from plots for this page
   }
}