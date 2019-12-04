# BTSPAS

Bayesian Time Stratified Petersen Analysis System

## Versions and installation

  * **Github** To install the latest development version from Github, 
    install the newest version of the **devtools** package; then run
```
devtools::install_github("cschwarz-stat-sfu-ca/BTSPAS", dependencies = TRUE,
                        build_vignettes = TRUE)
```

## Features
Provides advanced Bayesian methods to estimate
abundance and run-timing from temporally-stratified
Petersen mark-recapture experiments as described in Bonner and Schwarz(2011). 
Methods include
hierarchical modelling of the capture probabilities
and spline smoothing of the daily run size. This version 
uses JAGS to sample from the posterior distribution.

This is a more specific application of a stratified Petersen (Darroch, 1961; Plante et al 1988; Schwarz and Taylor, 1998) where, for example,
fish are tagged a a fish wheel over a series of weeks, migrate, and then are recaptured. The stratification is temporal rather than
geographical. This temporal stratification enable estimation of the run profile and various other parameters of interest such as abundance.

More information is available at
http://www.stat.sfu.ca/~cschwarz/Consulting/Trinity/Phase2



## References
Bonner Simon, J., & Schwarz Carl, J. (2011). 
Smoothing Population Size Estimates for Time-Stratified MarkRecapture Experiments Using Bayesian P-Splines. 
Biometrics, 67, 1498–1507.
https://doi.org/10.1111/j.1541-0420.2011.01599.x 

Darroch, J. N. (1961). The two-sample capture-recapture census when tagging and sampling are stratified. Biometrika, 48, 241–260.
https://www.jstor.org/stable/2332748

Plante, N., L.-P Rivest, and G. Tremblay. (1988). Stratified Capture-Recapture Estimation of the Size of a Closed Population. Biometrics 54, 47-60.
https://www.jstor.org/stable/2533994

Schwarz, C. J., & Dempson, J. B. (1994).
Mark-recapture estimation of a salmon smolt population. 
Biometrics, 50, 98–108.

Schwarz, C. J., & Taylor, C. G. (1998). The use of the stratified-Petersen estimator in fisheries management with an illustration of estimating the number of pink salmon (Oncorhynchus gorbuscha) that return to spawn in the Fraser River. Canadian Journal of Fisheries and Aquatic Sciences, 55, 281–296.
https://doi.org/10.1139/f97-238


