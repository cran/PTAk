\encoding{utf-8}
\name{COS2-CTR}
\alias{COS2}
\alias{CTR}
\title{Interpretation summaries}
\description{
After a \code{\link{FCA2}}, a \code{\link{SVDgen}}, a \code{\link{FCAk}} or a \code{\link{PTAk}} computes the traditional guides for interpretations used in PCA and correspondence analysis: COS2 or the percentage of variability rebuilt by the component and CTR or the amount of contribution towards that component.
}
\usage{
COS2(solu, mod=1, solnbs=2:4)
CTR(solu, mod=1, solnbs=1:4)
}
\arguments{
   \item{solu}{an object inheriting from class \code{PTAk}, representing a generalised  singular value decomposition}
   \item{mod}{an integer representing the mode number entry, 1 is row, 2 columns, ...}
   \item{solnbs}{a vector of integers representing the tensor numbers in the listing summary}  
    }

\details{
	Classical measures helping to interpret the plots in PCA, FCA and in PTAk as well.
	The sum of the COS2 across all the components needed to rebuild fully the tensor analysed) would make 1000 and the sum pf the CTR across the entry mode would be 1000.
}
 \value{a matrix whose columns are the COS2 or CTR  as per thousands (‰) for the mode considered }
 

\references{
Escoufier Y (1985) \emph{L'Analyse des correspondances : ses
\enc{propriétés}{propretes} et ses extensions}. ISI 45th session Amsterdam.

Leibovici D(1993) \emph{Facteurs \enc{à}{a} Mesures \enc{Répétées}{Repetees} et Analyses Factorielles :
applications \enc{à}{a} un suivi \enc{Epidémiologique}{Epidemiologique}}. \enc{Université}{Universite} de Montpellier
II. PhD Thesis in \enc{Mathématiques}{Mathematiques} et Applications (Biostatistiques).

Leibovici DG (2010) \emph{Spatio-temporal Multiway Decomposition using Principal Tensor Analysis on k-modes:the R package \pkg{PTAk}}. Journal of Statistical Software, 34(10), 1-34. \url{http://www.jstatsoft.org/v34/i10/}
}
\author{Didier G. Leibovici \email{c3s2i@free.fr}}
\seealso{\code{\link{PTAk}}, \code{\link{FCA2}}, \code{\link{FCAk}}, \code{\link{summary.FCAk}}, \code{\link{plot.PTAk}} }

\examples{
 data(crimerate)
  cri.FCA2 <- FCA2(crimerate)
   summary(cri.FCA2)
    plot(cri.FCA2, mod = c(1,2), nb1 = 2, nb2 = 3) # unscaled
    plot(cri.FCA2, mod = c(1,2), nb1 = 2, nb2 = 3, coefi = 
    	list(c(0.130787,0.130787),c(0.104359,0.104359)) ) # symmetric-map biplot
   CTR(cri.FCA2, mod = 1, solnbs = 2:4)
   CTR(cri.FCA2, mod = 2, solnbs = 2:4)
 COS2(cri.FCA2, mod = 2, solnbs = 2:4)
 
 
}
\keyword{models}
\keyword{principal components analysis}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
