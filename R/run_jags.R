## 2024-10-21 CJS added code to deal with bug in R2jags::jags() where pD is not placed in the BUGSoutput object if pD==0
## 2018-12-21 CJS converted to simple call to jags using R2jags package. This fixed problem that arrays not stored properly in output
## 2014-09-01 CJS added code to dump out mcmc.list to coda files to match functionality of OpenBugs
##                Set the seed here.
## 2013-12-31 CJS added code to dump out data and initial values to files as in run.openbugs
## 2013-12-30 CJS changed program argument in as.bugs.array to JAGS
## 2013-09-22 sjb Created file. Copied from run_openbugs.R

#' @import R2jags utils
#' @importFrom coda as.mcmc.list
#' @importFrom stats runif
#' @keywords internal
#' 


run.jags <-
    function(modelFile,
             dataFile,
             dataList,
             initFiles,
             initVals,
             parameters,
             nChains,
             nIter,
             nBurnin,
             nSims,
             overRelax,
             initialSeed,
             working.directory,
             debug=FALSE){
  cat("\n\n*** Start of call to JAGS \n")
  cat("Working directory: ",working.directory,"\n")

  nIterPostBurnin <- nIter - nBurnin   # Number of iterations to run after burnin
        
  nThin <- round(nIterPostBurnin/nSims) # Thinning to obtain desired number of samples
	#cat("In run_jags.R\n")
  #browser()

	## Set seed. We need to set a separate seed for each chain.
	## We start by setting the seed in R, and then generating nchain values between 1 and 10000000.
  cat("Initial seed for JAGS set to:", initialSeed, "\n")
  set.seed(initialSeed)
  initVals <- llply(initVals, function(x){
    # add to this list
    x$.RNG.seed <- round(stats::runif(1, min=1, max=1000000))
    x$.RNG.name <- "base::Wichmann-Hill"
    cat("Random number seed for chain ", x$.RNG.seed, "\n")
    x
  })

  # Dump out the data list (useful for debugging)
  #file.remove(dataFile)
  with(dataList, dump(names(dataList), file = dataFile))
	 
  #browser()

	# Dump out the initial values (useful for debugging)
  for (i in 1:nChains) {
    #file.remove(initFiles[i])
    initial.values <- initVals[[i]]
    with(initial.values, dump(names(initial.values), file = initFiles[i]))
  }
  #browser()     
  parametersToSave <- unique(c(parameters))

  results1 <- R2jags::jags( 
  #results1 <- my.jags( 
      data       =dataList,   # list of data variables
      inits      =initVals,   # list/function for initial values
      parameters =parametersToSave,# list of parameters to monitor
      model.file =modelFile, # where the model is saved by the cat above
      n.chains   =nChains,
      n.iter     =nIter,           # total iterations INCLUDING burn in
      n.burnin   =nBurnin,          # number of burning iterations
      n.thin     =nThin,                # how much to thin
      DIC=TRUE,                    # is DIC to be computed in both ways
      pD=TRUE, 
      jags.seed  = initialSeed,
      working.dir=working.directory      # store results in current working directory
      )

	## save the sample of the posteriors  as a coda file (for debugging)
	## taken from http://stackoverflow.com/questions/12078152/how-can-i-convert-an-mcmc-list-to-a-bugs-object
  s2 <- as.array(coda::as.mcmc.list(results1$BUGSoutput))
  lapply(seq_len(dim(s2)[3]), function(i) {
      write.table(cbind(rep(seq_len(nrow(s2[,,i])), ncol(s2)), c(s2[,,i])), 
              paste0(working.directory, '/CODAchain', i, '.txt'),
              row.names=FALSE, col.names=FALSE)
  })
  cat(paste(colnames(s2), 1+(seq_len(ncol(s2))-1) * nrow(s2), nrow(s2)*seq_len(ncol(s2))), 
               sep='\n', 
               file=file.path(working.directory, 'codaIndex.txt'))
  #browser()

  
  results <- results1$BUGSoutput
  
  # there is bug in the R2jags::jags code which does not include pD if pD is zero (despite the pD=TRUE being set)
  # see 
  # we need to add back the pD 
  if(is.null(results$pD)){
    results$pD <- 0
    results$DIC2 <- NA
  }
  
  results$model <- results1$model
  results$parameters.to.save <- results1$parameters.to.save
  ## Return results
  cat("\n\n*** Finished JAGS ***\n\n")

  return(results)
}
