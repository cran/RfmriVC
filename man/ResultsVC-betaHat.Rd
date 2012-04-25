\name{ResultsVC-betaHat}
\alias{betaHat}
\alias{betaHat,ResultsVC-method}
%%\alias{betaHat-methods}

\alias{betaHatNifti}
\alias{betaHatNifti,ResultsVC-method}
%%\alias{betaHatNifti-method}

\title{Calculation of the stimulus varying coefficient courses}

\description{
Postprocessing routine: Functions to calculate the course of the stimulus varying coefficient of a selected HRF basis function with corresponding significance/credibility intervals
}
\usage{
\S4method{betaHat}{ResultsVC}(object,voxelIndex3D,VC.plot=object@VC.plot,which.beta=1,alpha=object@alpha)

\S4method{betaHatNifti}{ResultsVC}(object,VC.plot=object@VC.plot,which.beta=1)
}


\arguments{
  \item{object}{An object of class \code{ResultsVC}.}
  \item{voxelIndex3D}{A 3D vector with the voxel index for which betaHat estimates should be calculated.}
  \item{VC.plot}{Vector of VC variable values for which betaHat estimates should be calculated.}
  \item{which.beta}{Number of requested basis function  (1 = canonical HRF, 2 = time derivative, 3 = dispersion derivative).}
  \item{alpha}{Significance resp. credibility level. Default: 0.05.}
}

\value{
\code{betaHat} returns
\item{betaHat}{Matrix of size length(VC.plot) x 3 with betaHat-estimate (1st column), lower and upper confidence interval limits (2nd + 3rd column) for a selected voxel.}

\code{betaHatNifti} returns
\item{betaHatTotal}{4D nifti image with betaHat vector estimate for each voxel.}
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

