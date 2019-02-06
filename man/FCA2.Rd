\encoding{utf-8}
\name{FCA2}
\alias{FCA2}
\title{Correspondence Analysis for 2-way tables}
\description{
Performs a particular \code{SVDgen}  data as a ratio Observed/Expected
under complete independence with metrics as margins of the
contingency table (in frequencies).
}
\usage{
 FCA2(X, nbdim =NULL, minpct = 0.01, smoothing = FALSE,
         smoo = rep(list(function(u) ksmooth(1:length(u), u, kernel = "normal",
        bandwidth = 3, x.points = (1:length(u)))$y), length(dim(X))),
      verbose = getOption("verbose"), file = NULL, modesnam = NULL,
    addedcomment = "", chi2 = FALSE, E = NULL, ...)}
\arguments{
   \item{X}{a matrix table of positive values}
   \item{nbdim}{a number of dimension to retain, if \code{NULL} the default value of maximum possible number of dimensions is kept}
   \item{minpct}{numerical 0-100 to control of computation of future solutions at this level and below}

   \item{smoothing}{see \code{\link{SVDgen}}}
   \item{smoo}{see \code{\link{SVDgen}}}
   \item{verbose}{control printing}
   \item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
   \item{modesnam}{character vector of the names of the modes, if \code{NULL} "\code{mo 1}"
     ..."\code{mo k}"}
    \item{addedcomment}{character string printed if \code{printt} after the title of the analysis}
   \item{chi2}{print the chi2 information when computing margins in \code{\link{FCAmet}}}
    \item{E}{ if not \code{NULL} is a matrix with the same dimensions as X with the same margins}
\item{...}{any other arguments passed to SVDGen or other functions}
    }

\details{
Gives the SVD-\emph{2}modes decomposition of the \eqn{1+\chi^2/N} of
the contingency table of full count \eqn{N=\sum X_{ij}},
i.e. complete independence + lack of independence (including marginal
independences) as shown for example in Lancaster(1951)(see reference
in Leibovici(1993 or 2000)). Noting \eqn{P=X/N},  a \code{SVD} of the
\eqn{(3)}-uple is done, that is :
 \deqn{ ((D_I^{-1} \otimes D_J^{-1})..P, \quad D_I, \quad D_J)}
where the metrics are diagonals of the corresponding margins. For
full description of arguments see \code{\link{PTAk}}. If \code{E}
is not \code{NULL} an FCAk-modes relatively to a model is
done (see Escoufier(1985) and therin reference
Escofier(1984) for a 2-way derivation), e.g. for a three way contingency table
\eqn{k=3} the \emph{4}-tuple data and metrics is:
 \deqn{ ((D_I^{-1} \otimes D_J^{-1} \otimes D_K^{-1})(P-E), \quad D_I, \quad D_J, \quad D_K)}
  If \code{E} was the complete independence (product of the margins)
  then this would give an \code{AFCk} but without looking at the
  marginal  dependencies (i.e. for a three way table no two-ways lack of
  independence are looked for).
}
 \value{a \code{FCA2} (inherits \code{\link{FCAk}} and \code{\link{PTAk}}) object }
 
\author{Didier G. Leibovici}

\references{
Escoufier Y (1985) \emph{L'Analyse des correspondances : ses
\enc{propriétés}{propretes} et ses extensions}. ISI 45th session Amsterdam.

Leibovici D(1993) \emph{Facteurs \enc{à}{a} Mesures \enc{Répétées}{Repetees} et Analyses Factorielles :
applications \enc{à}{a} un suivi \enc{Epidémiologique}{Epidemiologique}}. \enc{Université}{Universite} de Montpellier
II. PhD Thesis in \enc{Mathématiques}{Mathematiques} et Applications (Biostatistiques).

Leibovici D (2000) \emph{Multiway Multidimensional Analysis for
Pharmaco-EEG Studies}.\url{http://www.fmrib.ox.ac.uk/analysis/techrep/tr00dl2/tr00dl2.pdf}

Leibovici DG (2010) \emph{Spatio-temporal Multiway Decomposition using Principal Tensor Analysis on k-modes:the R package \pkg{PTAk}}. Journal of Statistical Software, 34(10), 1-34. \url{http://www.jstatsoft.org/v34/i10/}
}

\seealso{\code{\link{PTAk}}, \code{\link{FCAmet}}, \code{\link{summary.FCAk}}  }

\examples{
 data(crimerate)
 cri.FCA2 <- FCA2(crimerate)
 summary(cri.FCA2)
  plot(cri.FCA2, mod = c(1,2), nb1 = 2, nb2 = 3) # unscaled
  plot(cri.FCA2, mod = c(1,2), nb1 = 2, nb2 = 3, coefi = 
  	list(c(0.130787,0.130787),c(0.104359,0.104359)) )# symmetric-map biplot
 CTR(cri.FCA2, mod = 1, solnbs = 2:4)
 CTR(cri.FCA2, mod = 2, solnbs = 2:4)
 COS2(cri.FCA2, mod = 2, solnbs = 2:4)
 
 
}
\keyword{models}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
