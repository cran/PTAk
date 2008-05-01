\name{CANDPARA}
\alias{CANDPARA}
\title{CANonical DECOMPosition analysis and PARAllel FACtor analysis}
\description{Performs the identical models known as PARAFAC or CANDECOMP model.}
 \usage{
CANDPARA(X,dim=3,test=1E-8,Maxiter=1000,
                     smoothing=FALSE,smoo=list(NA),
                      verbose=getOption("verbose"),file=NULL,
                       modesnam=NULL,addedcomment="")}
\arguments{
\item{X}{a tensor (as an array) of order \emph{k}, if non-identity metrics are
used \code{X} is a list with \code{data} as the array and
\code{met} a list of metrics.}
\item{dim}{a number specifying the number of rank-one tensors }
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
tensors looked for is greater then one the methods differs:
PTA-\emph{k}modes will look for best approximation according to the
\emph{orthogonal rank} (\emph{i.e.} the rank-one tensors are
orthogonal), PCA-\emph{k}modes will look for best approximation
according to the \emph{space ranks} (\emph{i.e.} the ranks of all
(simple) bilinear forms , that is the number of components in each
space), PARAFAC/CANDECOMP will look for best approximation according
to the \emph{rank} (\emph{i.e.} the rank-one tensors are not
necessarily orthogonal). For sake of comparisons the
PARAFAC/CANDECOMP method and the PCA-\emph{n}modes are also in the
package but complete functionnality of the use these methods and more
complete packages may be checked at the www site quoted below.
 }
\value{a \code{CANDPARA} (inherits from \code{\link{PTAk}}) object}
 \references{
Caroll J.D and Chang J.J (1970) \emph{Analysis of individual
differences in multidimensional scaling via n-way generalization of
'Eckart-Young' decomposition}. Psychometrika 35,283-319.

Harshman R.A (1970) \emph{Foundations of the PARAFAC procedure:
models and conditions for 'an explanatory' multi-mode factor
analysis}. UCLA Working Papers in Phonetics, 16,1-84.

Kroonenberg P (1983) \emph{Three-mode Principal Component Analysis:
Theory and Applications}. DSWO press. Leiden.(related references in
\url{http://three-mode.leidenuniv.nl})

Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component Analysis of
multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329.
 }
\note{
The use of metrics (diagonal or not) and smoothing extends
flexibility of analysis. This program runs slow! A PARAFAC orthogonal
can be done with PTAk looking only for k-modes Principal Tensors
\emph{i.e.} with the options \code{nbPT=c(rep(0,k-2),dim), nbPT2=0}.
It is identical to look in any \code{PTAk} decomposition only for the
\emph{k}modes solution but obviously with unecessary computations.}

\author{Didier Leibovici \email{c3s2i@free.fr}}

\keyword{array }  \keyword{algebra}
\keyword{multivariate}
