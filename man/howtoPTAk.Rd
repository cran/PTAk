\name{howtoPTAk}
\alias{howtoPTAk}
\title{howto for Principal Tensors Analysis of a k-modes Tensor}
\description{A mini guide to handle PTAk model decomposittion}
 \usage{
howtoPTAk()}

\arguments{}
 \details{The PTAk decomposition aims at building an approximation of a given multiway data, represented as a tensor, based on a variance criterion. This approximation is given by a set of rank one tensors, orthogonals to each other, in a nested algorithm process and so controlling the level of approximation by the amount of variability extracted and represented by the sum of squares of the singular values (associated to the rank one tensors). In that respect it offers a way of generalising PCA to tensors of order greater than 2. \cr
 \cr
 Full description of the licence is given in the file copyrights given in the distribution, or alternatively in the readme file  along with the history of changes. This function displays the copyright and licence information.}

\references{ Leibovici D and Sabatier R (1998) \emph{A Singular Value Decomposition of a k-ways
array for a Principal Component Analysis of multi-way data, the PTA-k}. Linear Algebra and its
Applications, 269:307-329 }

\author{Didier Leibovici \email{c3s2i@free.fr}}
\seealso{\code{\link{PTA3}}, \code{\link{PTAk}} ,\code{\link{FCAk}}}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
