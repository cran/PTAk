\name{TENSELE}
\alias{TENSELE}
\title{Elementary Tensor product}
\description{
Computes the Tensor Product of a list of vectors  (or matrices) according to
a given order.}
 }
 \usage{
 TENSELE(T,moins=NULL, asarray=TRUE,order=NULL,id=NULL)}
 \arguments{
   \item{T}{a list like a \code{PTAk} object and minimally  just contains \code{v} }
   \item{moins}{if not \code{NULL}, vector of indexes  (in the list T) to skip  }
   \item{asarray}{logical to specify the output form \code{TRUE} gives an array, \code{FALSE} gives a vector}
   \item{order}{
if not \code{NULL} vector of length \code{length(T)}, \code{NULL} is equivalent to
\code{length(T):1} as the function makes indexes in order run slowest to
fastest }
   \item{id}{
when \code{T} is a list of matrices, can be either a vector of
\code{length(T)} giving indexes of the vectors for each space (following
order) or a list of vectors of indexes.} }
\details{
The tensor product of the vectors (or matrices) in the list \code{T} is
computed, skipping or not the indexes in \code{moins}, the output tensor is
either in tensor form or in vector form. The way the tensor product is done
follows \code{order}.}
\value{
According to \code{asarray} the value is either an array, or a vector
representing the tensor product of the vectors (not in moins), the dimension
in \code{order[1]} running the slowest.}

\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
 \seealso{\code{\link{REBUILD}}}

\keyword{array }  \keyword{algebra}
