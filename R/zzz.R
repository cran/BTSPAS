# 2009-12-06 CJS test if the standard winbugs/openbugs directory exists and warn the user.
# 2009-12-01 CJS added openbugs/winbugs to arguments of all functions. No need for global variables
#

.First.lib <- function(libname,pkgname){

  .First.windows <- function(){
    cat("***** BTSPAS: Bayesian Time Stratified Petersen Analysis System (2010-04-26) ***** \n\n",
        "      Help available with  help(package='BTSPAS')  or  \n",
        "                           help(BTSPAS)                \n\n",
        '      List of demonstrations available with demo(package="BTSPAS") \n',
        '      A demo is run as (e.g.): demo("demo-TSPDE", package="BTSPAS")\n\n')
    # Turn off the ask options when demonstrations are run
    options(demo.ask=FALSE)
    # check for existence of the standard install locations for open/win bugs
    path <-file.path("c:","Program Files","OpenBugs")
    cat("OpenBUGS assumed to be at:",path,"\n")
    if( 0 != file.access(path)){
       cat("***** WARNING ***** Could not find the OpenBugs filepath: \n   ",path,
           "\n You may need to specify where OpenBugs is found in the call to the *_fit programs if you are using OpenBugs \n")
    }
    path <-file.path("c:","Program Files","WinBUGS14")
    cat("WinBUGS assumed to be at:",path,"\n")
    if( 0 != file.access(path)){
       cat("***** WARNING ***** Could not find the WinBUGS filepath: \n   ",path,
           "\n You may need to specify where WinBUGS is found in the call to the *_fit program if you are using WinBUGS \n")
    }

  }

  .First.other <- function(){
    cat("It appears that you are running R on ",.Platform$OS.type,".\n",
        "The functions in BTSPAS make use of the BRugs package.\n",
        "which is only available for the Windows operating\n",
        "system. Unfortunately, this means that BTSPAS can only be run\n",
        "on Windows at this time. We hope to be able make a *nix version\n",
        "of the software soon.\n\n",sep="")
  }
    
  if(.Platform$OS.type=="windows")
    .First.windows()

  else
    .First.other()
}



  

  
