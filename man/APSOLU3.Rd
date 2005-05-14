\name{APSOLU3}
\alias{APSOLU3}
\title{Associated 3-modes Principal Tensors of a 3-modes
Principal Tensor}
\description{ Computes all the 2-modes solutions associated to the given
Principal Tensor of the given tensor.  }
 \usage{
APSOLU3(X,solu,pt3=NULL,nbPT2=1,
                 smoothing=FALSE,smoo=list(NA),
                        verbose=getOption("verbose"),file=NULL )}

\arguments{
\item{X}{a tensor (as an array) of order \emph{3}, if non-identity metrics are
used \code{X} is a list with \code{data}  as the array and
\code{met} a list of metrics }
\item{solu}{a \code{\link{PTAk}} object}
\item{pt3}{a number identifying in \code{solu} the Principal Tensor to use or the
last (if \code{NULL})}
\item{nbPT2}{integer, if 1 all solutions will be computed otherwise at maximum
nbPT2  solutions}
\item{smoothing}{see \code{\link{SVDgen}}}
\item{smoo}{see \code{\link{PTA3}}}
\item{verbose}{control printing}
 \item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
   }
 \details{ For each component of the identified Principal Tensor
given in \code{solu}, an SVD of the contracted product of \code{X} and the component is done.
This gives all the associated Principal Tensors which updates \code{solu} supposed to contain
Principal Tensors of \code{X}. }
\value{an updated \code{\link{PTAk}} object}
\references{ Leibovici D and Sabatier R (1998) \emph{A Singular Value Decomposition of a k-ways
array for a Principal Component Analysis of multi-way data, the PTA-k}. Linear Algebra and its
Applications, 269:307-329 }

 \note{ Usually (i.e. as in \code{PTA3} and \code{PTAk}) the principal
tensor used is the first Principal Tensor of X (and is the last updated in \code{solu}). If
it is another Principal Tensor, the obtained associated solutions do not \emph{stricto sensu}
refer to the SVD-\emph{k}modes decomposition (because the orthogonality is defined in the whole
tensor space not necessarily on each component space) but are still meaningful. }
\author{Didier Leibovici \email{c3s2i@free.fr}}
\seealso{\code{\link{PTA3}}, \code{\link{APSOLUk}} }
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
