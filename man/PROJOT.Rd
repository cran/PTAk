\name{PROJOT}
\alias{PROJOT}
\title{Orthogonal Tensor projection }
\description{
Orthogonal-tensoriel projection of a tensor according to a rank-1
tensor, or a to bigger structure defined by kronecker product of
matrices. }
\usage{
 PROJOT(X,solu,numo=1,bortho=TRUE,Ortho=TRUE,metrics=NULL)
}
\arguments{
   \item{X}{a tensor(as an array) of any order }
   \item{solu}{an object like a \code{solutions.PTAk} object with at least \code{v}}
   \item{numo}{
a vector of numbers  or a list of vectors  (length the order of the
tensor) identifying for each space the structure to project onto, if
NULL for a specific space then no projection is done for this space}
\item{bortho}{list of logicals saying if the structures are othogonal
or not.}
  \item{Ortho}{
list of logicals telling the projectors on each space to be  on the
structure  or on its orthogonal.}
  \item{metrics}{NULL or list of metrics (either diagonal or not) for each entry of \code{X}}
 }
\details{
This function computes the \emph{tensorial orthogonal projection} of
\code{X} onto the \emph{tensorial structure} defined by \code{solu}
and \code{numo}. For each space (involved in the tensorial product
where from  \code{X} belongs), one defined the projector onto
\code{solu[[i]][["v"]][numo,]} (or on its orthogonal if
\code{Ortho[[i]]==TRUE}), then the result is the image of \code{X} by
the tensorial product of the projectors, i.e.
 \deqn{(P_{S1} \otimes P_{S2} \otimes \ldots \otimes P_{Sk})(X)}.
 }
\value{A tensor with dimensions as \code{X}}
\references{
Leibovici D (1993) \emph{Facteurs à Mesures Répétées et Analyses
Factorielles : applications à un suivi épidémiologique}. Université
de Montpellier II PhD Thesis in Mathématiques et Applications
(Biostatistiques).

Leibovici D (2000) \emph{Multiway Multidimensional Analysis for
Pharmaco-EEG Studies}.(submitted) \url{http://www.fmrib.ox.ac.uk/~didier/cv/recentpub.html}
                     }
\note{
For PTA-\emph{k}modes the projection used is only on rank-one tensors
(Principal Tensors), \emph{i.e.} \code{numo} is a number. The code
here can be used for any structure (on each spaces) and constitutes
the projector onto a tensorial structure, and can define the
PTAIV-\emph{k}modes (PTAk on Instrumental Variables Leibovici(1993).
(see other references  for tensorial product of linear operators in
Leibovici(2000) \emph{e.g.} Dauxois et al.(1994))}
\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
\seealso{\code{\link{PTAk}} }
\examples{

don <- array(1:360,c(5,4,6,3))
 don <- don + rnorm(360,10,2)

 ones <- list(list(v=rep(1,5)),list(v=rep(1,4)),list(v=rep(1,6)),list(v=rep(1,3)))
 donfc <- PROJOT(don,ones)

 apply(donfc,c(2,3,4),mean)
apply(donfc,c(1),mean)

  # implementation de PTAIVk with obvious settings
 PTAIVk <- function(X,STruct,...)
         {X <- PROJOT(X$data,STruct,numo=Struct[[1]]$numo,Ortho=Struct[[1]]$Ortho,metrics=X$met)
         PTAk(X,...)
         }

   }
\keyword{array  algebra}
\keyword{multivariate}
