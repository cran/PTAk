\name{REBUILD}
\alias{REBUILD}
\title{Build an approximation of the tensor of any order}
\description{
Gives the approximation of a previously analysed tensor using its
given decomposition. }
\usage{
 REBUILD(solutions,nTens=1:2,testvar=1 ,redundancy=FALSE)}
\arguments{
   \item{solutions}{a \code{solutions.PTAk} object}
   \item{nTens}{a vector of identifying positions (numbers given in \code{summary}) for  the choice of Principal Tensors to use}
   \item{testvar}{control within \code{nTens} used Principal Tensor with  minimum percent of variability explained}
   \item{redundancy}{
logical to take into account (within \code{nTens}) PT \emph{tested}
redundant during analysis (seealso \code{RESUM}) if \code{TRUE}.}
      }
\details{
The function rebuilds the Principal Tensors, \emph{i.e.} rank-one
tensors of order the order of the tensor analysed, and add them up to
build an approximation of the tensor analysed (according to the
method used see \code{method}). This constitutes a best Least Squares
(ordinary or "weighted"  if metrics are used) approximation of
  \code{datanam} for a given \emph{orthogonal-rank r} (number of
 principal tensors used), if and only if the singular values used are
 the \emph{r} highest.
   }
   \value{A tensor with dimensions as \code{solutions[[k]][["datanam"]]}. }
 \note{
 This function can be called for \code{PARAFAC/CANDECOMP} and
 \code{PCAn}. A specific rebuilt is implemented for this last one. }

 \author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}

  \seealso{\code{\link{PTAk}}}

\keyword{multivariate}
\keyword{models}
