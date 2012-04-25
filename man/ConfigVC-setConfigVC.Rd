\name{ConfigVC-setters}
\alias{setConfigVC}
\alias{setConfigVC,ConfigVC-method}
\alias{setConfigVC,method}

\alias{setConfigVC.Bayes}
\alias{setConfigVC.Bayes,ConfigVC-method}
\alias{setConfigVC.Bayes,method}

\alias{setConfigVC.Classical}
\alias{setConfigVC.Classical,ConfigVC-method}
\alias{setConfigVC.Classical,method}


\title{
ConfigVC-Class: Constructor resp. setter functions
}
\description{
These functions serve as constructors for \code{ConfigVC} objects, i.e. these should be used to assign values to object slots.
}
\usage{
\S4method{setConfigVC}{ConfigVC}(object,tr=2,stimTimes,nStimBF=1,VC,nVCknots=10,
		degreeVCspline=2,confounder=NULL,freqHPF=128,y,mask=NULL,grandMeanScaling=TRUE)
		
\S4method{setConfigVC.Classical}{ConfigVC}(object)

\S4method{setConfigVC.Bayes}{ConfigVC}(object,chainSize=1100,burnin=100,step=5,sigma2start=1,
		sigma2IG=list(a=0.001,b=0.001),VCcontrol=NULL,diagSample=NULL)		
}

\arguments{
  	\item{object}{Object of class \code{ConfigVC}.}
    \item{tr}{Repetition time of fMRI scans in seconds. Default: \code{tr=2}}
    \item{stimTimes}{Numeric vector with stimulus presentation times of one stimulus type which should have a varying coefficient.}
    \item{nStimBF}{Number of canonical HRF basis functions: 1 = just canonical HRF, 2 = plus time derivative, 3 = plus dispersion derivative. Default: \code{nStimBF=1}}
    \item{VC}{Numeric vector with variable values the VC depends on. }
    \item{nVCknots}{Number of knots for the varying coefficient (VC) spline. Default: \code{nVCknots=10} }
    \item{degreeVCspline}{ Degree of the VC spline. Default: \code{degreeVCspline=2} }
    \item{confounder}{Matrix with confounder values. }
    \item{freqHPF}{Frequency of high pass filter (HPF). Default: \code{freqHPF=128}}
    \item{y}{4d fMRI dataset in form of an object of class \code{"nifti"}. }
    \item{mask}{3D analysing mask in form of an object of class \code{"nifti"}. Per default a threshold mask is clculated. }
    \item{grandMeanScaling}{Indicator if a grand mean scaling should be applied to the fMRI data.}
	\item{chainSize}{Total number of MCMC iterations to be conducted. Default: \code{chainsize=1100}}
    \item{burnin}{Number of iterations after which convergency of Markov chains to their equilibrium distribution is assumed. Default: \code{burnin=100} }
    \item{step}{Thinning parameter:  After the burnin phase one out of \code{step}-th samples of the Markov chains is used for estimation. Default: \code{step=5}}
    \item{sigma2start}{Starting parameter for each voxels residual variance 'sigma2'.}
    \item{sigma2IG}{A list with two elements containing the inverse gamma shape and scale parameter of the 'sigma2' prior distribution.}
    \item{VCcontrol}{A list with with \code{nStimBF} elements. Each element is a list with configurations for one HRF basis function. See details. }
    \item{diagSample}{Vector with 1D indices of voxels for which all diagnostics should be safed.}
   
}
\details{
For each HRF basis function, which is requested to be contained in the model, the \code{VCcontrol}-list contains a list with the following default: list(gammaStart=matrix(0,nrow=p,ncol=1),a=0.001,b=0.001).
If non-default values should be used, the whole \code{VCcontrol} statement must be given, e.g. for a model with canonical HRF and its time derivative (\code{nStimBF=2}):

\code{VCcontrol}=list(list(gammaStart=matrix(gammaLSest1,nrow=p,ncol=1),a=0.05,b=0.001), list(gammaStart=matrix(gammaLSest2,nrow=p,ncol=1),a=0.1,b=0.1))

where parameter p is the number of VC coefficients of each HRF (p<-nVCknots+degreeVCspline-1) and gammaLSest1 and gammaLSest2 are for example the estimates from a model with constant stimulus effects.  
}
\value{
  \item{\code{object}}{Object of class \code{ConfigVC}.}
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

