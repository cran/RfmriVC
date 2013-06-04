\name{ConfigVC-class}
\Rdversion{1.1}
\docType{class}
\alias{ConfigVC-class}
%\alias{doFmriVC,ConfigVC-method}
%\alias{doFmriVCbayes,ConfigVC-method}
%\alias{doFmriVCclassical,ConfigVC-method}
%\alias{setConfigVC.Bayes,ConfigVC-method}
%\alias{setConfigVC.Classical,ConfigVC-method}
%\alias{setConfigVC,ConfigVC-method}
\alias{show,ConfigVC-method}


\title{Class \code{"ConfigVC"}}
\description{
This is a class for configurating the RfmriVC-algorithm. The slots of the \code{"ConfigVC"}-class contain all arguments and helper variables which are passed to the RfmriVC-algorithm.
Two setter functions are provided to structure the setting of slots resp. arguments with regard to content. If kept unspecified in calls of the setters, slots are assigned default values. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ConfigVC", ...)}.
To initializes an object use setter functions \code{setConfigVC} and afterwards, if a bayesian analysis is requested \code{setConfigVC.Bayes}. 
}
\section{Slots}{
  \describe{
    \item{\code{tr}:}{Repetition time of fMRI scans in seconds.}
    \item{\code{stimTimes}:}{Numeric vector with stimulus presentation times of one stimulus type which should have a varying coefficient.}
    \item{\code{nStimBF}:}{Number of canonical HRF basis functions: 1 = just canonical HRF, 2 = plus time derivative, 3 = plus dispersion derivative.}
    \item{\code{nVCknots}:}{Number of knots for the varying coefficient (VC) spline. }
    \item{\code{degreeVCspline}:}{ Degree of the VC spline. }
    \item{\code{VC}:}{Numeric vector with variable values the VC depends on. }
    \item{\code{confounder}:}{Matrix with confounder values. }
    \item{\code{freqHPF}:}{Frequency of high pass filter (HPF), e.g. 128s. }
    \item{\code{y}:}{4d fMRI dataset in form of an object of class \code{"nifti"}. }
    \item{\code{mask}:}{3D analysing mask in form of an object of class \code{"nifti"}. }
    \item{\code{mask2totalIndex}:}{Matrix with mapping from mask index to total index, i.e. information about the i-th selected voxel is contained in row i, which contains the 3D index of the corresponding voxel.}
    \item{\code{total2maskIndex}:}{Object of class \code{"nifti"}, which contains at i,j,k the index of the voxel within the set of selected voxels if selected by mask or zero otherwise.}
    \item{\code{grandMeanScaling}:}{Indicator if a grand mean scaling should be applied to the fMRI data.}
    \item{\code{U}:}{Design matrix part with nuisance covariables (HPF baseline and confounders). }
    \item{\code{X}:}{Design matrix part with  }
    \item{\code{penMat}:}{Penalty matrix for one VC spline. }
    \item{\code{method}:}{The method to be applied for estimation and inference. Either "classical" (penalized least squares) or "bayes". }
    \item{\code{methodParams}:}{If \code{method}=="bayes", a list with configuration settings for the bayesian method. }
  }
}
\section{Methods}{
	The following methods operate on objects of class \code{ConfigVC}. Details of these methods are described on separate manual pages.
  \describe{
    \item{doFmriVC}{\code{signature(object = "ConfigVC")}: ... }
    \item{doFmriVC.Bayes}{\code{signature(object = "ConfigVC")}: ... }
    \item{doFmriVC.Classical}{\code{signature(object = "ConfigVC")}: ... }
    \item{setConfigVC.Bayes}{\code{signature(object = "ConfigVC")}: ... }
    \item{setConfigVC.Classical}{\code{signature(object = "ConfigVC")}: ... }
    \item{setConfigVC}{\code{signature(object = "ConfigVC")}: ... }
    \item{show}{\code{signature(object = "ConfigVC")}: ... }
    \item{saveConfig}{\code{signature(object = "ConfigVC")}: ... }
	 }
}
\references{
Bothmann, L. (2012). Statistische Modellierung von EEG-abhaengigen Stimuluseffekten in der fMRT-Analyse. Diploma thesis. LMU Munich: Germany. 
}
\author{Ludwig Bothmann and Stefanie Kalus}

\seealso{
Setter functions for assigning slot values (constructors): \cr 
\code{\link{setConfigVC}}, \code{\link{setConfigVC.Bayes}}, \code{\link{setConfigVC.Classical}}\cr\cr
Running the estimation procedure: \cr
\code{\link{doFmriVC}}, \code{\link{doFmriVC.Classical}}, \code{\link{doFmriVC.Bayes}}\cr\cr  
Saving and loading an external representation of a ConfigVC-object: \cr 
\code{\link{saveConfig}}, \code{\link{loadConfig}}
}

\keyword{classes}
