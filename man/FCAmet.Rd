\name{FCAmet}
\alias{FCAmet}
\title{Tool used in Generalisation of Correspondence Analysis for k-way tables}
\description{
Computes the  ratio Observed/Expected under complete independence
with margins of the multiple contingency table (in frequencies) and
gives \code{chi2} statistic of lack of complete independence. }
\usage{ FCAmet(X,chi2=FALSE,E=NULL,No0margins=TRUE)}
\arguments{
   \item{X}{a multiple contingency table (array) of order \emph{k}}
   \item{chi2}{if \code{TRUE} prints the chi2 statistic information }
   \item{E}{if not \code{NULL} represent a model which
   would be used for an \code{\link{FCAk}} relatively to a model}
   \item{No0margins}{if TRUE, prevents zero margins in replacing cells involved by the min of the non-zero margins /nb of zero cells}
    }

\value{a list with
\item{data}{
an array \code{(X/count (-E))/Indepen} where \code{Indepen} is the array obtained
from he products of the margins}
\item{met}{
a list wherein each entry is the vector of the corresponding margins i.e.
\code{apply(X,i,sum)/count}}
\item{count}{is the total sum \code{sum(X).}}
 }
\note{The statistics and metrics do not depend on \code{E}. The statistic given measure only the lack of independence.}

\author{Didier G. Leibovici}
\seealso{ \code{\link{FCAk}} }
\keyword{ multivariate}
\keyword{models}
