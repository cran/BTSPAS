# 2018-11-27 CJS - removed call to OpenBugs. Left function here in case we want to change samplers in the future.
# 2013-12-30 CJS - function to switch between the three samplers as needed. 
#                  This way we wont't have to modify much code (hopefully)

#' @keywords internal
run.MCMC <-
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
             overRelax=FALSE,
             initialSeed,
             working.directory,
             debug=FALSE){

   results <- run.jags(modelFile=modelFile,
                           dataFile=dataFile,
                           dataList=dataList,
                           initFiles=initFiles,
                           initVals=initVals,
                           parameters=parameters,
                           nChains=nChains,
                           nIter=nIter,
                           nBurnin=nBurnin,
                           nSims=nSims,
                           overRelax=overRelax,
                           initialSeed=initialSeed,
                           working.directory=working.directory,
                           debug=debug)
return(results)
} # end of function

