\name{INITIA}
\alias{INITIA}
\title{Initialisation used in SINGVA}
\description{Gives the first Tucker1 components of a given tensor.
}
\usage{
   INITIA(X,modesnam=NULL,method="Presvd",dim=1,...)}
\arguments{
    \item{X}{a tensor (as an array) of order \emph{k}}
    \item{modesnam}{a character vector of the names of the modes}
    \item{method}{
uses either the inbuilt SVD \code{method="svd"} or a power
algorithm giving only the first \code{method="Presvd"} or
any other function given applying to the column space of a
matrix and returning a list with \code{v} (in columns
vectors as in \code{svd}) and \code{d}. }
    \item{dim}{
default 1 in each space otherwise specify the number of dimensions
e.g. \code{c(2,3..,2)} (with \code{"Presvd"} dim is obviously 1)}
\item{...} {extra arguments of the method \code{method}: the first argument is fixed (see details). }}
 \details{
Computes the first (or \code{dim}) right singular vector  (or other
summaries) for every representation of the tensor as a matrix with
\code{dim(X)[i]} columns, \code{i=1...k}.}

\value{
a list (of length \emph{k}) of  lists  with arguments:
   \item{v}{the singular vectors in rows}
   \item{modesnam}{a character object naming the mode, \code{"m i"} otherwise}
   \item{n}{labels of mode \code{i} entries as given in \code{dimnames} of the data, can be \code{NULL}}
   \item{d}{the corresponding first singular values}
   }
\references{
Kroonenberg P.M (1983) \emph{Three-mode Principal Component Analysis:
Theory and Applications}.  DSWO Press, Leiden.

Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component Analysis of
multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329. }
\note{
The collection these eigenvectors, is known as the Tucker1 solution
or initialisation related to PCA-3modes or PCA-\emph{n}modes models.
If a function is given  it may include \code{dim} as argument.}
\author{Didier Leibovici \email{c3s2i@free.fr}}
\seealso{ \code{\link{SINGVA}}, \code{\link{PTAk}} }
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
