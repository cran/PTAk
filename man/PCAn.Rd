\name{PCAn}
\alias{PCAn}
\title{Principal Component Analysis on n modes }
\description{Performs the Tucker\emph{n} model using a space version of RPVSCC (\code{SINGVA}).}
\usage{
 PCAn(X,dim=c(2,2,2,3),test=1E-12,Maxiter=400,
              smoothing=FALSE,smoo=list(NA),
                verbose=getOption("verbose"),file=NULL,
                  modesnam=NULL,addedcomment="")
}
\arguments{
   \item{X}{
a tensor (as an array) of order \emph{k}, if non-identity metrics are
used \code{X} is a list with \code{data}  as the array and
\code{met} a list of metrics}
   \item{dim}{a vector of  numbers specifying the dimensions in each space }
\item{test}{control of convergence}
\item{Maxiter}{maximum number of iterations allowed for convergence}
\item{smoothing}{see \code{\link{SVDgen}}}
\item{smoo}{see \code{\link{PTA3}}}
\item{verbose}{control printing}
\item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
\item{modesnam}{character vector of the names of the modes, if \code{NULL}
  "\code{mo 1}" ..."\code{mo k}"}
\item{addedcomment}{character string printed after the title of the analysis}

   }
\details{
Looking for the best rank-one tensor approximation (LS) the three
methods described in the package are equivalent. If the number of
tensors looked for is greater then one the  methods differs:
PTA-\emph{k}modes will look for best approximation according to the
\emph{orthogonal rank} (\emph{i.e.} the rank-one tensors are
orthogonal), PCA-\emph{k}modes will look for best approximation
according to the \emph{space ranks} (\emph{i.e.} the rank of every
bilinear form, that is the number of components in each space),
PARAFAC/CANDECOMP will look for best approximation according to the
\emph{rank} (\emph{i.e.} the rank-one tensors are not necessarily
orthogonal). For the sake of comparisons the PARAFAC/CANDECOMP method
and the PCA-\emph{n}modes are also in the  package but complete
functionnality  of the use these methods  and more complete packages
may be fetched at the www site quoted below. }
\value{a \code{solutions.PCAn} (inherits \code{\link{solutions.PTAk}}) object}
\references{
Caroll J.D  and Chang J.J (1970) \emph{Analysis of individual
differences in multidimensional   scaling via n-way generalization of
"Eckart-Young" decomposition}. Psychometrika 35,283-319.

Harshman R.A (1970) \emph{Foundations of the PARAFAC procedure:
models and conditions for "an explanatory" multi-mode factor
analysis}. UCLA Working Papers in Phonetics, 16,1-84.

Kroonenberg P (1983) \emph{Three-mode Principal Component Analysis:
Theory and Applications}. DSWO press. Leiden.(related references in
\url{http://www.fsw.leidenuniv.nl/~kroonenb/})

Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component Analysis of
multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329.
              }
\note{The use of metrics (diagonal or not) and smoothing extend flexibility of analysis. }
\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
\keyword{models}
\keyword{array  algebra}
\keyword{multivariate}
