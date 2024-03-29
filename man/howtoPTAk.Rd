\name{howtoPTAk}
\alias{howtoPTAk}
\title{howto for Principal Tensors Analysis of a k-modes Tensor}
\description{A mini guide to handle PTAk model decomposition}
 \usage{
howtoPTAk()}
 \details{
 The PTAk decomposition aims at building an approximation of a given multiway data, represented as a tensor, based on a variance criterion. This approximation is given by a set of rank one tensors, orthogonal to each other, in a nested algorithm process and so controlling the level of approximation by the amount of variability extracted and represented by the sum of squares of the singular values (associated to the rank one tensors). In that respect it offers a way of generalising PCA to tensors of order greater than 2.


The reference in JSS provides details about preparing a dataset and running a general PTAk and particularities for spatio-temporal data. Some aspects on FCAk can also be found in the NCRM publication.
 
 
 
 The license is GPL-3, support can be provided via http://c3s2i.free.fr, donations via Paypal to c3s2i@free.fr are welcome.}

\references{ Leibovici D and Sabatier R (1998) \emph{A Singular Value Decomposition of a k-ways
array for a Principal Component Analysis of multi-way data, the PTA-k}. Linear Algebra and its
Applications, 269:307-329 

Leibovici DG (2010) \emph{Spatio-temporal Multiway Decomposition using Principal Tensor Analysis on k-modes:the R package \pkg{PTAk}}. Journal of Statistical Software, 34(10), 1-34. \doi{10.18637/jss.v034.i10}

Leibovici DG and Birkin MH (2013) Simple, multiple and multiway correspondence analysis applied to spatial census-based population microsimulation studies using R. NCRM Working Paper. NCRM-n^o 07/13, Id-3178 \url{https://eprints.ncrm.ac.uk/id/eprint/3178}
}

\author{Didier G. Leibovici}
\seealso{\code{\link{PTA3}}, \code{\link{PTAk}} ,\code{\link{FCAk}}}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
