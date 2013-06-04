\name{plotNifti}
\alias{plotNifti}
\title{Function to plot brain image slices}
\description{
A function to plot the layers of one 3D nifti image stored in a nifti object of library \code{Rniftilib}.
}
\usage{
plotNifti(nifti.img, zRange=NULL, zVec=NULL,  mfrowVal, plotTitle=TRUE,  
          pdfoutput=FALSE, pdfName=NULL, scale=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nifti.img}{
A nifti object of library \code{Rnfitilib} which stores the 3D (4D) nifti image.
}
  \item{zRange}{
An optional two-dimensional integer vector with a start and end layer to be plotted.
}
  \item{zVec}{
An optional integer vector with all layer numbers to be plotted.
}
  \item{mfrowVal}{
An optional 2D vector in form of \code{c(nRow,nCol)}. This is for \code{par(mfrow=mfrowVal)}.
}
  \item{plotTitle}{
A bool value stating whether a title should be plotted above image plots.
}
  \item{pdfoutput}{
A bool value stating whether a pdf should be generated.
}
  \item{pdfName}{
An optional string with a name for the pdf file.
}
  \item{scale}{
An optional scale value which can be used to scale the values of the plotted image.
}
\item{...}{Optional arguments which can be passed to \code{par(...)}}
}

\note{If a 4D image is given the nifti image only the first 3D nifti image is plotted.}

\author{Stefanie Kalus <stefanie.kalus <at> stat.uni-muenchen.de>}

\examples{
\dontrun{
require("Rniftilib",character.only=TRUE)
activMap<-nifti.image.read("activMap.nii")
plotNifti(activMap,pdfoutput=TRUE)
}
}


