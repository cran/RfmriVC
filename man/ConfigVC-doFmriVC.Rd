\name{ConfigVC-doFmriVC}
\alias{doFmriVC}
%%\alias{doFmriVC-method}
\alias{doFmriVC,ConfigVC-method}

\alias{doFmriVC.Bayes}
%%\alias{doFmriVC.Bayes-method}
\alias{doFmriVC.Bayes,ConfigVC-method}

\alias{doFmriVC.Classical}
%%\alias{doFmriVC.Classical-method}
\alias{doFmriVC.Classical,ConfigVC-method}




\title{
Function to fit a varying stimulus coefficient within an fMRI regression model.
}
\description{
Depending on the configuration settings in a \code{ConfigVC} object, these functions start the estimation process: 
\code{doFmriVC(object)} either starts the bayesian or classical procedure depending on the value of slot \code{method}.
Alternatively, the classical procedure can be started directly by calling \code{doFmriVC.Classical(object)} 
and the bayesian procedure can be strated directly by calling \code{doFmriVC.Bayes(object)}.  
}
\usage{
\S4method{doFmriVC}{ConfigVC}(object)

\S4method{doFmriVC.Classical}{ConfigVC}(object)

\S4method{doFmriVC.Bayes}{ConfigVC}(object)
}

\arguments{
  \item{object}{Object of class \code{ConfigVC}.}
}

\value{
resultsObj
  \item{\code{resultsObj}:}{Object of class \code{ResultsVC}.}
}
\references{
Bothmann, L. (2012). Statistische Modellierung von EEG-abhaengigen Stimuluseffekten in der fMRT-Analyse. Diploma thesis. LMU Munich: Germany. 
}
\author{Ludwig Bothmann and Stefanie Kalus}

\seealso{
	\code{\linkS4class{ConfigVC}}, \link[=RfmriVC-package]{RfmriVC}
}
\examples{
\dontrun{
See package RfmriVC documentation: package?RfmriVC
}
}