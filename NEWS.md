
# BTSPAS 2024=04-01

* CRAN fixups to correct error in example. No change in functionality.

# BTSPAS 2023-03-31

* Documentation updates. No change in functionality.

# BTSPAS 2011.11.2

* Added trunc.logitP to truncate logit(P) when plotting to avoid problems with extreme cases (issue #30)
* Editorial changes

# BTSPAS 2011.11.1

* Fixed importing sd() from stat package that conflicts with revised sd() in actuar package
* Minor bug fixes to deal with extreme cases
* Editorial changes

# BTSPAS 2021.1.1

* User now able to specify prior for beta parameters of logitP vs covariates.
Refer to the vignette on covariates.
* All values of n1, m2 must be non-missing and positive. If you set the bad.n1
and bad.m2 values, then n1 and m2 are all set to 0.
* Fixed slight error in computation of gof statistics for the non-diagonal cases. We now exclude 
strata that are flagged as being bad from the computations.

# BTSPAS 2020.9.1

* Fixed a bug in trace plots (from ggforce::facet_wrap_paginate) 
when the last page doesn't have a full set of plots. Dummy plots added.
* Secondary axis added on right side of logit(P) fit and log(U) fit
 
# BTSPAS  2020.2.1

* Added a vignette about forcing the run to 0 at the start and end of the study.

# BTSPAS 2020.1.1

* Runtimings now added to output object.
* More testing of fall back options.

# BTSPAS 2019.0101
 
* sampfrac argument deprecated. Only values of 1 now accepted. This never did work properly and will be
removed in future releases.
* OpenBugs no longer supported.
* Diagonal model now allows you to fix some P's -- typically used when second trap is not running.
* All plots converted to ggplot format and stored in the plots object in the final result. An argument is available in the fitting routines to stop saving output to files. 
* Bayesian p.value plot created for MarkAvail routines (issue #17) 
* Several vignettes created to replace and enhance demo files.


