# 2020-08-07  CJS ggforce::facet_wrap_paginate() has a bug where it fails if the last page has only a few plots
#                 we added some dummy parameters to the last page (zzz_dummy1, etc)

#' Creates trace plots of specified parameters showing the multiple chains and
#' the value of Rhat
#' 
#' Takes the MCMC object returned from a split and produces trace_plots for the
#' listed parameters. It shows a separate line on the plot for each chain and
#' also shows the value of Rhat
#' 
#' 
#' @param title A character string used for a title on reports and graphs
#' @param results The MCMC object containing the results from the call to JAGS
#' @param parms_to_plot A character vector of names of parameters to plot.
#' These must match exactly to the parameter names used in the simulation.
#' @param ncol,nrow How many plots to put on a page (number of rows and columns)
#' @return List of ggplot2 objects using facet_wrap_paginate (...., page=...) with each element of the list
#' corresponding to one page of the plot. 
#' @template author
#' @import ggforce reshape2 ggplot2 plyr
#' @keywords internal

#' @examples
#' \dontrun{
#' # Create trace plots of the logitP parameters
#' # 
#' # Trace plots of logitP
#' varnames <- names(results$sims.array[1,1,])  # extract the names of the variables 
#' trace.plot <- plot_trace(title=title, 
#'                          results=results, 
#'                          parms_to_plot=varnames[grep("^logitP", varnames)])
#' if(save.output.to.files){
#'   pdf(file=paste(prefix,"-trace-logitP.pdf",sep=""))
#'   plyr::l_ply(trace.plot, function(x){plot(x)})
#'   dev.off()
#'}
#' } % end of dontrun
#' 

plot_trace <- function(title=" ", results=NULL, parms_to_plot=NULL, nrow=2, ncol=2){
#
# Takes the MCMC object from the fit (could be TPSDE etc), a list of parameters and produces
# the traceplots.
#   
# title - title of the plot
# results - the MCMC object containing the necessary information
# parms_to_plot - character vector containing the names of the parms to plot
#                e.g. c("logitP[1]", "logitP[2]")
#               - this should be an exact match#

  varnames <- colnames(results$sims.matrix)
  index <- match(parms_to_plot, varnames) # find where these parms exist in the array

  trace.df <- reshape2::melt(results$sims.array[,,varnames[index]],
                      varnames=c("Simulation","Chain","Parameter"),
                      value.name="Value")
  npages <- ceiling(length(index)/ncol/nrow)
  
  # we need to add "extra" panels to deal with facet_wrap_paginate() error when 
  # last page has only a single plot
  #browser()
  if(npages*nrow*ncol > length(unique(trace.df$Parameter))){
     trace.df <- plyr::rbind.fill(trace.df,
                                  data.frame(Parameter=paste0("zzz_dummy",1:(npages*nrow*ncol - length(unique(trace.df$Parameter)))),
                                             Simulation=1,
                                             Value=0,
                                             Chain=1, stringsAsFactors=TRUE # be sure that ordering is kept of parameter names
                                             )
     )
  }

  allplots <-plyr::llply(1:npages, function(page){
     ggplot(data=trace.df, aes_(x=~Simulation, y=~Value, color=~as.factor(Chain)))+
              ggtitle(title)+
              geom_line()+
              scale_color_discrete(name="Chain")+
              ggforce::facet_wrap_paginate(~Parameter, ncol=ncol, nrow=nrow, page=page, scales="free_y")
  })
  allplots
} # end of function

