\name{RfmriVC-save}
\alias{saveConfig}
\alias{saveConfig,ConfigVC-method}

\alias{saveResults}
\alias{saveResults,ResultsVC-method}

\alias{loadConfig}
\alias{loadResults}

\title{Functions to save and load \code{ConfigVC} and \code{RfmriVC} objects from the RfmriVC package}
\description{
\code{saveConfig} writes an external representation of a \code{ConfigVC} object to the specified file.
\code{loadConfig} reloads datasets written with the function \code{saveConfig}.
\code{saveResults} writes an external representation of a \code{ResultsVC} object to the specified file.
\code{loadResults} reloads datasets written with the function \code{saveResults}.
}
\usage{
\S4method{saveConfig}{ConfigVC}(object,path)

loadConfig(path)

\S4method{saveResults}{ResultsVC}(object,path)

loadResults(path)
}

\arguments{
  \item{object}{Object of class \code{ConfigVC} or \code{ResultsVC}.}
  \item{path}{Path to folder, where objects should be saved resp. was saved.}
}

\value{
  \item{\code{object}:}{Object of either class \code{ConfigVC} or \code{ResultsVC}.}
}

\details{Extra functions are necessary because otherwise the nifti-members are not written to file.}

\references{
Bothmann, L. (2012). Statistische Modellierung von EEG-abhaengigen Stimuluseffekten in der fMRT-Analyse. Diploma thesis. LMU Munich: Germany. 
}
\author{Ludwig Bothmann and Stefanie Kalus}

\seealso{
	\code{\link{save}}, \code{\link{load}}, \code{\linkS4class{ConfigVC}}, \code{\linkS4class{ResultsVC}}, \link[=RfmriVC-package]{RfmriVC}
}
\examples{
\dontrun{
See package RfmriVC documentation: package?RfmriVC
}
}