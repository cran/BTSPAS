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

