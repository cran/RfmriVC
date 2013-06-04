\name{RfmriVC-package}
\alias{RfmriVC-package}
\docType{package}
\title{RfmriVC}
\description{
The RfmriVC package provides utilities for fitting a varying stimulus coefficient model within a voxelwise fMRI regression model. 
In each voxelwise model, the effects of all requested canonical HRF basis functions are allowed to vary depending on the values of a chosen variable.
The model fit can be done either with bayesian or classical (penalized least squares by \code{MGCV}) techniques.  
}
\details{
This package builds upon two \code{S4} classes which control the configuration, running and postprocessing of the varying coefficient fMRI model.
An object of class \code{ConfigVC} is used to set all configurations needed to run the algorithm (including output of preprocessing steps).
An Object of class \code{ResultsVC} is returned by a model run (return of \code{doFmriVC()}) and can be used to visualize results. 
}

\author{Ludwig Bothmann and Stefanie Kalus}

\references{
Bothmann, L. (2012). Statistische Modellierung von EEG-abhaengigen Stimuluseffekten in der fMRT-Analyse. Diploma thesis. LMU Munich: Germany. 
}

\keyword{ package }

\examples{
\dontrun{
library(RfmriVC)

#------- Read data ----------
# Read in fMRI data
y<-nifti.image.read("nifti_fmri.nii")

#Stimulus presentation times
stimTimes<-as.matrix(read.table("stimulusTimes.txt"))

#Confounders
confMat<-as.matrix(read.table("covariates.txt"))


#------- Classical ------------
# Generate configuration object for classical analysis (= default analysis)
configObj<-new("ConfigVC")
configObj<-setConfigVC(configObj,stimTimes=stimTimes,VC=stimTimes,y=y,
                        nStimBF=3,confounder=confMat)
# Do classical analysis
resObj<-doFmriVC(configObj)

#Look at results of first voxel (1,1,1)
par(mfrow=c(3,1),mar=c(3, 4, 0, 1) + 0.1,oma=c(0, 0, 5, 0))

#Plot stimulus varying coefficient depending on variable VC: beta(VC)
betaHatEst<-betaHat(resObj,voxelIndex3D=c(1,1,1))
matplot(resObj@VC.plot,betaHatEst,type="l",xlab="VC",ylab=expression(beta~(VC)))

#Plot HRF for VC=100
mat<-hrfHat(resObj,VC=100,t=seq(0,30,length.out=1000),voxelIndex3D=c(1,1,1))
matplot(resObj@t.plot,mat,type="l",xlab="t",ylab="hrf(t)")

#Plot HRF for t=5
matVC<-hrfHat(resObj,t=5,VC=seq(7,390,length.out=1000),voxelIndex3D=c(1,1,1))
matplot(resObj@VC.plot,matVC,type="l",xlab="VC",ylab="hrf(VC)")
title("Classical estimation",outer=TRUE,line=2)

#Image plot for HRF surface
windows()
t.grid<-seq(1,15,length.out=50)
VC.grid<-seq(7,670,length.out=50)
matVCt<-hrfHat(resObj,t=t.grid,VC=VC.grid,voxelIndex3D=c(1,1,1))
image(y=VC.grid,x=t.grid,t(matVCt[,,1]),main="Classical estimation")


#------- Bayes ------------
set.seed(1112)
# Generate configuration object for bayesian analysis 
# (= non-default analysis -> use extra setter)
configBayObj<-new("ConfigVC")
configBayObj<-setConfigVC(configBayObj,stimTimes=stimTimes,VC=stimTimes,y=y,
                          nStimBF=3,confounder=confMat)
configBayObj<-setConfigVC.Bayes(configBayObj)

#Do analysis
resBayObj<-doFmriVC(configBayObj)

#Look at results
windows()
par(mfrow=c(3,1),mar=c(3, 4, 0, 1) + 0.1,oma=c(0, 0, 5, 0))

#Plot stimulus varying coefficient depending on variable VC: beta(VC)
betaHatEst<-betaHat(resBayObj,voxelIndex3D=c(1,1,1))
beta.hat.plot(resBayObj@VC.plot,betaHatEst)

#Plot HRF for VC=100
mat<-hrfHat(resObj,VC=100,t=seq(0,30,length.out=1000),voxelIndex3D=c(1,1,1))
matplot(resObj@t.plot,mat,type="l",xlab="t",ylab="hrf(t)")

#Plot HRF for t05
matVC<-hrfHat(resObj,t=5,VC=seq(7,390,length.out=1000),voxelIndex3D=c(1,1,1))
matplot(resObj@VC.plot,matVC,type="l",xlab="VC",ylab="hrf(VC)")
title("Bayesian estimation",outer=TRUE,line=2)

#Image plot for HRF surface
windows()
t.grid<-seq(1,15,length.out=50)
VC.grid<-seq(7,670,length.out=50)
matVCt<-hrfHat(resObj,t=t.grid,VC=VC.grid,voxelIndex3D=c(1,1,1))
image(y=VC.grid,x=t.grid,t(matVCt[,,1]),main="Bayesian estimation")


}}
