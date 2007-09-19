\name{APSOLUk}
\alias{APSOLUk}
\title{Associated k-modes Principal Tensors of a k-modes Principal Tensor}
\description{
Computes all the (k-1)-modes associated solutions to the given
Principal Tensor of the given tensor. Calls recursively PTAk. }
\usage{
 APSOLUk(X,solu,nbPT,nbPT2=1,
                       smoothing=FALSE,smoo=list(NA),
                        minpct=0.1,ptk=NULL,
                         verbose=getOption("verbose"),file=NULL,
                          modesnam=NULL)}

\arguments{

\item{X}{
a tensor (as an array) of order \emph{k}, if non-identity metrics are
used \code{X} is a list with \code{data}  as the array and
\code{met} a list of metrics}
\item{solu}{a \code{PTAk}  object }
\item{nbPT}{a number or a vector of dimension \emph{(k-2)}}
\item{nbPT2}{integer, if 0 no 2-modes solutions will be computed, 1 means all, >1 otherwise}
\item{smoothing}{see \code{\link{SVDgen}}}
\item{smoo}{see \code{\link{PTA3}}}
\item{minpct}{numerical 0-100 to control of computation of future solutions at this level and below}
\item{ptk}{a number identifying in solutions the Principal Tensor to use or the last (if \code{NULL})}
\item{verbose}{control printing}
 \item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
\item{modesnam}{character vector of the names of the modes, if \code{NULL} "\code{mo 1}" ..."\code{mo k}"}
}

\details{
For each component of the identified  Principal Tensor given in
\code{solutions}, a PTA-\emph{(k-1)}modes of the contracted product
of X and the component is done. This gives all the associated
Principal Tensors which updates  \code{solutions} supposed to contain
a Principal Tensors of X at the first place. For full description of
arguments see \code{\link{PTAk}}.}

\value{
an updated \code{PTAk} object }                                  

\references{
Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component Analysis of
multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329. }

\note{
Usually (\emph{i.e.} as in \code{PTA3} and \code{PTAk}) the
principal tensor used is the first Principal Tensor of
\code{X} (and is the last updated in \code{solutions}). If
it is another Principal Tensor, the obtained associated
solutions do not \emph{stricto sensu} refer to the
SVD-\emph{k}modes decomposition (because the orthogonality
is defined in the whole tensor space not necessarily on
each component space) but are still meaningful. This
function is usually called by \code{PTAk} but can be used
on its own to carry on a \code{PTAk} analysis if \code{X}
is the projected (of the original data) on the orthogonal
of all the \emph{k}modes Principal Tensor. In other words
the \code{ptk} rank-one tensor in \code{solutions} should
be the first best rank-one tensor approximating \code{X}
for this decomposition analysis to be called
PTA-\emph{k}modes.}

\author{Didier Leibovici \email{c3s2i@free.fr}}
\seealso{\code{\link{PTAk}} }
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
