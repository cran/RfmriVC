\name{calcDeviceHeight}
\alias{calcDeviceHeight}
\title{Calculates the optimal graphical device height for printing squared voxels}
\description{
This helper function calculates for a given width value, e.g., in \code{pdf(file,width, height, \dots)} the according height value, so that voxels appear as squares.
}
\usage{
calcDeviceHeight(widthVal, mfrowVal, maiVal, omiVal, niftiDim)
}

\arguments{
  \item{widthVal}{Intended \code{width} value for the graphical device.}
  \item{mfrowVal}{The par-parameter \code{mfrow} intended to be used for plotting.}
  \item{maiVal}{The par-parameter \code{mai} intended to be used for plotting.}
  \item{omiVal}{The par-parameter \code{omi} intended to be used for plotting.}
  \item{niftiDim}{A 3D vector with the 3D dimensions of the nifti image object.}
}
\value{
\item{heightVal}{The \code{height} value for a given \code{width} value, which leads to square voxels in a graphical device.}
}
\author{Stefanie Kalus <stefanie.kalus <at> stat.uni-muenchen.de>}
\examples{
#Plotting parameters
dimVec<-c(43,40,19)
mfrowVal<-c(4,5)
maiVec<-c(0,0, 0,0)+0.015
omiVec<-c(0,0,0.65,0)+0.15 #Outer title margin
    
#Calculate appropriate width and height, so that image ratios are preserved
widthVal<-7
heightVal<-calcDeviceHeight(widthVal,mfrowVal,maiVec,omiVec,dimVec)
heightVal
}
