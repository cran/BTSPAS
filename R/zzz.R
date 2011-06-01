# 2011-06-01 SB  remove usage of winbugs; change initial seed; speed up mixing
# 2009-12-06 CJS test if the standard winbugs/openbugs directory exists and warn the user.
# 2009-12-01 CJS added openbugs/winbugs to arguments of all functions. No need for global variables
#

.First.lib <- function(libname,pkgname){

  cat("***** BTSPAS: Bayesian Time Stratified Petersen Analysis System - Version 2011.06 (2011-06-01) ***** \n\n",
        "      Help available with  help(package='BTSPAS')  or  \n",
        "                           help(BTSPAS)                \n\n",
        '      List of demonstrations available with demo(package="BTSPAS") \n',
        '      A demo is run as (e.g.): demo("demo-TSPDE", package="BTSPAS")\n\n')
  ## Turn off the ask options when demonstrations are run
  options(demo.ask=FALSE)
}







