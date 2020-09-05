\name{NEWS}
\title{BTSPAS News}
\encoding{UTF-8}

%\section{CHANGES IN VERSION 2012.xxxx}{
%  \subsection{USER-VISIBLE CHANGES}{
%    \itemize{
%       \item None.
%    }
%  }
%   \subsection{BUG FIX}{
%     \itemize{
%       \item Add na.rm=TRUE to any() and all() function calls in data checking.
%     }
%   }
%}

\section{CHANGES IN VERSION 2020.9.1}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Fixed a bug in trace plots (from ggforce::facet_wrap_paginate) 
       when the last page doesn't have a full set of plots. Dummy plots added.
      \item Secondary axis added on right side of logit(P) fit and log(U) fit
     }
   }
}
\section{CHANGES IN VERSION 2020.2.1}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added a vignette about forcing the run to 0 at the start and end of the study.
     }
   }
}
\section{CHANGES IN VERSION 2020.1.1}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Runtimings now added to output object.
      \item More testing of fall back options.
     }
   }
}


\section{CHANGES IN VERSION 2019.0101}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item sampfrac argument deprecated. Only values of 1 now accepted. This never did work properly and will be
            removed in future releases.
      \item OpenBugs no longer supported.
      \item Diagonal model now allows you to fix some P's -- typically used when second trap is not running.
      \item All plots converted to ggplot format and stored in the plots object in the final result. An argument is available in the 
      fitting routines to stop saving output to files. 
      \item Bayesian p.value plot created for MarkAvail routines (issue #17) 
      \item Several vignettes created to replace and enhance demo files.
    }
  }
   \subsection{BUG FIX}{
     \itemize{
       \item  Corrected error in creating u2copy to improve mixing. Issue #18..
       \item  Corrected error where output of Theta in non-diagonal cases was stored incorrectly in sims.matrix and sim.array objects. Code used sim.list where Theta were stored correctly, so no effect on fitting.
     }
   }
}


\section{CHANGES IN VERSION 2015.0702}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Bayesian P-values plots are now created using ggplot() and saved in results.
    }
  }
   \subsection{BUG FIX}{
     \itemize{
       \item  Bayesian p-value plots were incorrect for hatchery/wild mixtures.
       \item  Updated Imports/Suggests to conform to latest CRAN standards
     }
   }
}

\section{CHANGES IN VERSION 2014.0901}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
        \item Now uses JAGS rather than OpenBugs so BTSPAS should run on all platforms.
              OpenBUGS can still be run using the 'engine' argument on the *_fit() call.
      \item Minor changes to package file structures to meet R package conventions. No change in 
             functionality.
      \item Minor bug fixes discovered in conversion to JAGS.
 
    }
  }
   \subsection{BUG FIX}{
     \itemize{
       \item Initialized missing values in n1, u2, m2 to sensible values to avoid OpenBugs from 
             wandering too far off and generating nonsense values.
       \item Corrected initialization of epsilon when logitP.fixed is set 
     }
   }
}


\section{CHANGES IN VERSION 2012.0217}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
       \item Minor changes to package file structures to meet R package conventions. No change in 
             functionality.
    }
  }
   \subsection{BUG FIX}{
     \itemize{
       \item None.
     }
   }
}


\section{CHANGES IN VERSION 2012.0215 and 2012.0216}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
       \item Plots in non-diagonal case use all time periods rather than 
             just those from releases.
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item Cleanup of some function calls; Fix arguments for xlim in some
            plots; don't try and print Rhat if number of chains is 1.
    }
  }
}

\section{CHANGES IN VERSION 2012.0202}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
       \item NONE
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item Fixed and issue where if you had missing values in your data, 
            or zero in n1 and m2 (which
            the program automatically changes to NA), the scaling factor
            to compute initial values was returned as NA which caused lm()
            to fail when initial values are generated. 
    }
  }
}

\section{CHANGES IN VERSION 2012.02}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
       \item Posterior Predictive Plots modified to make axes symmetric.
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item Predictive posterior plots (Bayesian P-value plot) did not have symmetric X/Y axis so the p-value may not print.
            The axes are now forced to be symmetric.
      \item Misc changes to correct typos, documentation, etc. 
    }
  }
}


\section{CHANGES IN VERSION 2012.01}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
       \item BTSPAS now uses OpenBugs exclusively so it can run on Windoze, Mac, and Unix. Users may need to install
             OpenBugs in place of WinBugs. No other change in functionality.
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item NONE.
    }
  }
}

