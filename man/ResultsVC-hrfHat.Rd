\name{ResultsVC-hrfHat}
\alias{hrfHat}
\alias{hrfHatPoint}
\alias{hrfHat,ResultsVC-method}
\alias{hrfHatPoint,ResultsVC-method}


\title{
Calculation of HRF estimates
}
\description{
Postprocessing routine: Functions to calculate the HRF estimates with corresponding significance / credibility intervals.

\code{hrfHatPoint()} calculates the hrf estimate for one pair of t and VC value.
\code{hrfHat()} calculates the hrf estimate for scalar or vector valued t and VC settings.
}
\usage{
\S4method{hrfHatPoint}{ResultsVC}(object,VCVal=median(object@VC),tVal,voxelIndex3D,alpha=object@alpha)

\S4method{hrfHat}{ResultsVC}(object,VC=median(object@VC),t=object@t.plot,voxelIndex3D,alpha=object@alpha)
}

\arguments{
  \item{object}{An object of class \code{ResultsVC}.}
  \item{VCVal}{Varying coefficient variable value for which hrf(VC,t) should be calculated.}
  \item{VC}{Either a value or a vector of the VC variable for which hrf(VC,t) should be calculated.}
  \item{tVal}{Peristimulus time value for which hrf(VC,t) should be calculated.}
  \item{t}{Either a value or a vector of (peristimulus) time for which hrf(VC,t) should be calculated.}
  \item{voxelIndex3D}{A 3D vector with the voxel index for which betaHat estimates should be calculated.}
  \item{alpha}{Significance resp. credibility level. Default: 0.05.}
}

\value{
\code{hrfHatPoint} returns
\item{hrfHatPoint}{Vector with hrfHat-estimate (1st element), lower and upper confidence interval limits (2nd + 3rd elements) for the requested t and VC value for a selected voxel.}


\code{hrfHat} returns
\item{hrfHat}{Vector, matrix or array with hrfHat-estimate (1st element/column/layer), lower and upper confidence interval limits (2nd + 3rd elements/columns/layers) for a selected voxel.}
}

\references{
Bothmann, L. (2012). Statistische Modellierung von EEG-abhaengigen Stimuluseffekten in der fMRT-Analyse. Diploma thesis. LMU Munich: Germany. 
}
\author{Ludwig Bothmann and Stefanie Kalus}

\seealso{
	\code{\linkS4class{ResultsVC}}, \link[=RfmriVC-package]{RfmriVC}
}
\examples{
\dontrun{
See package RfmriVC documentation: package?RfmriVC
}
}

