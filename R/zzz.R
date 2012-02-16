# 2012-01-10 CJS Change name of .First.lib to  .onLoad with changes for R2.14+
# 2011-06-01 SB  remove usage of winbugs; change initial seed; speed up mixing
# 2009-12-06 CJS test if the standard winbugs/openbugs directory exists and warn the user.
# 2009-12-01 CJS added openbugs/winbugs to arguments of all functions. No need for global variables
#

#.First.lib <- function(libname,pkgname){
.onLoad <- function(libname,pkgname){

  packageStartupMessage("***** BTSPAS: Bayesian Time Stratified Petersen Analysis System - Version 2012.0215 (2012-02-15) ***** \n\n",
        "      Help available with  help(package='BTSPAS')  or  \n",
        "                           help(BTSPAS)                \n\n",
        '      List of demonstrations available with demo(package="BTSPAS") \n',
        '      A demo is run as (e.g.): demo("demo-TSPDE", package="BTSPAS")\n\n')
  ## Turn off the ask options when demonstrations are run
  options(demo.ask=FALSE)
}







