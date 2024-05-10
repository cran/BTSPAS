# 2015-06-11 CJS Updated for fixed to Bayesian p-value plots
# 2013-12-30 CJS Updated for JAGS
# 2013-01-25 CJS Changed messages from .onLoad to .onAttach in accordance with R policies
# 2012-01-10 CJS Change name of .First.lib to  .onLoad with changes for R2.14+
# 2011-06-01 SB  remove usage of winbugs; change initial seed; speed up mixing
# 2009-12-06 CJS test if the standard winbugs/openbugs directory exists and warn the user - now removed.
# 2009-12-01 CJS added openbugs/winbugs to arguments of all functions. No need for global variables
#

#' Message to display when package is loaded
#' 
#' @keywords internal
#'

.onAttach <- function(libname,pkgname){

  packageStartupMessage("***** BTSPAS: Bayesian Time Stratified Petersen Analysis System - Version 2024-05-09 (2024.5.9) ***** \n\n",
        "      Help available with  help(package='BTSPAS') \n",
        '      Several vignettes are available. See browseVignettes(package="BTSPAS") \n\n')
}







