\name{CONTRACTION}
\alias{CONTRACTION}
\alias{CONTRACTION.list}
 \title{Contraction of two tensors }
\description{Computes the contraction product of two tensors as a generalisation of matrix product.
 }
\usage{
 CONTRACTION(X,z, Xwiz=NULL,zwiX=NULL,rezwiX=FALSE,usetensor=TRUE)
 CONTRACTION.list((X,zlist,moins=1,zwiX=NULL,usetensor=TRUE,withapply=FALSE)
}
 \arguments{
   \item{X}{a tensor(as an array) of any order}
   \item{z}{another tensor  (with at least one space in common)}
   \item{zlist}{
a list of lists like a \code{solution.PTAk} at least with
\code{v} for every list(here \code{v} can be any array) }
   \item{Xwiz}{
\code{Xwiz} is to specify the entries of \code{X} to
contract with entries of \code{z} specified by \code{zwiX},
if \code{Xwiz} \code{NULL} \code{dim(z)[zwiX]} matching
\code{dim(X)} will do without ambiguity (taking all
\code{z} dimensions if \code{zwiX} is \code{NULL}). In
\code{CONTRACTION.list} it is not set as one supposes the
contractions in the list to operate follow the dimensions
of X}
   \item{zwiX}{
idem as \code{Xwiz}. If both \code{Xwiz} and \code{zwiX}
are \code{NULL} \code{zwiX}is replaced by full
possibilities (\code{1:length(dimz)}) then \code{Xwiz} is
looked for. In \code{CONTRACTION.list} it is the vector for
dimensions in the \code{v} to contract with X. Only 1-way
dimension for each \code{v}.}
   \item{moins}{
the elements in \code{zlist} to skip (see also
\code{\link{TENSELE}})}
   \item{rezwiX}{
logical if \code{TRUE} (and zwiX is \code{NULL}) rematches
the dimensions in for \code{zwiX}: useful only if the
dimensions of z were not following the Xwiz order and are
not equals.}
   \item{usetensor}{if \code{TRUE} uses \code{tensor} (add-on package)}
    \item{withapply}{if \code{TRUE} (only for vectors in \code{zlist} uses \code{apply}}
   }
\details{
Like two matrices \emph{contract} according to the
appropriate dimensions (columns matching rows) when one
performs a matrix product, this operation does pretty much
the same thing for tensors(array) and specified contraction
dimensions given by \code{Xwiz} and \code{zwiX} which
should match. The function is actually written like:
transforms both tensors as matrices with the ``matching
tensor product" of their contraction dimensions in columns
(for higher order tensor) and rows (the other one),
performs the matrix product and rebuild the result as a
tensor(array). Without using \code{tensor}, if \code{Xwiz}
and/or \code{zwiX} are not specified the functions tries to
match all \code{z} dimensions onto the dimensions of X
where X is the higher order tensor (if it is not the case
in the arguments the function swaps them). }
\value{
A tensor of dimension \code{c(dim(X)[-Xwiz],dim(z)[-zwiX])}
if \code{X} has got a bigger order than \code{z}.}
\references{
Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component
Analysis of multi-way data, the PTA-k}. Linear Algebra and
its Applications, 269:307-329.

Schwartz L (1975) \emph{Les Tenseurs}. Herman, Paris. }

\note{
This operation generalises the \emph{matrix} product to the
\emph{contracted} product of any two tensors(arrays), and
should theoretically perform the tensor product if no
matching (no contraction) but has not been implemented. I
recently put the option of using \code{tensor} which does
exactly the same thing faster as well as it is from
\code{C}. When using \code{tensor} if \code{Xwiz} or
\code{zwiX} are \code{NULL} they are replaced by the full
possibilities.}

\author{Didier Leibovici \email{c3s2i@free.fr}}

 \seealso{\code{\link{PTAk}},  \code{\link{APSOLUk}} }

 \examples{
 library(tensor)
   z <-  array(1:12,c(2,3,2))
   X <- array(1:48,c(3,4,2,2))
   Xcz <- CONTRACTION(X,z,Xwiz=c(1,3,4),zwiX=c(2,3,1))
   dim(Xcz)   # 4
   Xcz1 <- CONTRACTION(X,z,Xwiz=c(3,4),zwiX=c(1,3))
   dim(Xcz1) # 3,4,3
   Xcz2 <- CONTRACTION(X,z,Xwiz=c(3,4),zwiX=c(3,1))
   Xcz1[,,1]
   Xcz2[,,1]
   #######
   sval0 <- list(list(v=c(1,2,3,4)),list(v=rep(1,3)),list(v=c(1,3)))
   tew <- array(1:24,c(4,3,2))
    CONTRACTION.list(tew,sval0,moins=1)
       #this is equivalent to the following which may be too expensive for big datasets
    CONTRACTION(tew,TENSELE(sval0,moins=1),Xwiz=c(2,3))
   ##
     CONTRACTION.list(tew,sval0,moins=c(1,2)) #must be equal to
     CONTRACTION(tew,sval0[[3]]$v,Xwiz=3)
   }
\keyword{array} \keyword{algebra}
