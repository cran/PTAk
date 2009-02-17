\name{datasets}
\alias{timage12}
\alias{crimerate}
\alias{Zone_climTUN}
\title{data used for demo in SVDgen, PTA3}
\description{The \code{crimerate} dataset provides crime rates per 100,000 people in
 seven categories for each of the fifty states (USA) in 1977.  The \code{timage12} dataset
  is an image from fMRI analysis (one brain slice), it is a \emph{t}-statistic image over 12 subjects of the activation (verbal) parameter.
   The \code{Zone_climTUN} is an object of class Map  representing montly (\code{12}) measurements in  Tunisia  
of \code{10} climatic indicators. The grid of \code{2599} cells was stored previously as a shapefile and read using \code{read.shape}. }
\usage{
data(crimerate)
data(timage12)
data(Zone_climTUN)
}
\format{
\code{crimerate} is a matrix of \code{50 x 7} for the \code{crimerate} data.

\code{timage12} is a matrix \code{91 x 109} for \code{timage12} data.}

\source{\code{crimerate} comes from SAS. The \code{timage12} comes from FMRIB center, University of Oxford. 
 The \code{Zone_climTUN} comes from WorldCLIM database 2000 see references along with description of the indicators in Leibovici et al.(2007).
 }
 \references{
Leibovici D, Quillevere G, Desconnets JC (2007). A Method to Classify Ecoclimatic Arid
and Semi-Arid Zones in Circum-Saharan Africa Using Monthly Dynamics of Multiple
Indicators. IEEE Transactions on Geoscience and Remote Sensing, 45(12), 4000-4007.}

\keyword{datasets}
