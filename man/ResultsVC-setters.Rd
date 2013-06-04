\name{ResultsVC-setters}
\alias{resetResultsVC,ResultsVC-method}
\alias{resetResultsVC}
\alias{setResultsVC,ResultsVC-method}
\alias{setResultsVC}


\title{
Setter functions for class \code{ResultsVC}
}
\description{
These functions serve as constructors for \code{ResultsVC} objects, i.e. these should be used to assign values to object slots.
\code{setResultsVC()} is primarily used internally and is called by \code{doFmriVC()}.
\code{resetResultsVC()} can be used to modify slots which are of interest for the user.
}

\usage{
\S4method{setResultsVC}{ResultsVC}(object,configVCobj,method,methodResults)

\S4method{resetResultsVC}{ResultsVC}(object, VC.plot=seq(min(configVCobj@VC),max(configVCobj@VC),
	length.out=1000), t.plot=seq(0,30,length.out=1000), alpha=0.05)
}

\arguments{
  \item{object}{An object of class \code{ResultsVC}.}
  \item{configVCobj}{An object of class \code{ConfigVC}.}
  \item{method}{The method which has been applied for estimation and inference. Either "classical" (penalized least squares) or "bayes". }
  \item{methodResults}{A list with estimation results.}
  \item{VC.plot}{Vector of VC variable values for which betaHat estimates should be calculated.} 
  \item{t.plot}{Vector of (peristimulus) time values for which hrf estimates should be calculated. }
  \item{alpha}{Significance resp. credibility level. Default: 0.05.}
}

\value{
\item{object}{An object of class \code{ResultsVC}.}
}

\references{
Bothmann, L. (2012). Statistische Modellierung von EEG-abhaengigen Stimuluseffekten in der fMRT-Analyse. Diploma thesis. LMU Munich: Germany. 
}
\author{Ludwig Bothmann and Stefanie Kalus}

\seealso{
	\code{\linkS4class{ResultsVC}}, \link[=RfmriVC-package]{RfmriVC}
}
%\examples{
%\dontrun{
%See package RfmriVC documentation: package?RfmriVC
%}
%}

