\name{FCAk}
\alias{FCAk}
\title{Generalisation of Correspondence Analysis for k-way tables}
\description{
Performs a particular \code{PTAk}  data as a ratio Observed/Expected
under complete independence with metrics as margins of the multiple
contingency table (in frequencies).
}
\usage{
 FCAk(X,nbPT=3,nbPT2=1,minpct=0.01,
               smoothing=FALSE,smoo=rep(list(
                       function(u)ksmooth(1:length(u),u,kernel="normal",
                       bandwidth=3,x.points=(1:length(u)))$y),length(dim(X))),
                     verbose=getOption("verbose"),file=NULL,
                       modesnam=NULL,addedcomment="",chi2=TRUE,E=NULL)}
\arguments{
   \item{X}{a multiple contingency table (array) of order \emph{k}}
   \item{nbPT}{a number or a vector of dimension \emph{(k-2)}}
   \item{nbPT2}{if 0 no 2-modes solutions will be computed, 1 =all, >1
     otherwise}
   \item{minpct}{numerical 0-100 to control of computation of future solutions at this level and below}

   \item{smoothing}{see \code{\link{SVDgen}}}
   \item{smoo}{see \code{\link{SVDgen}}}
   \item{verbose}{control printing}
   \item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
   \item{modesnam}{character vector of the names of the modes, if \code{NULL} "\code{mo 1}"
     ..."\code{mo k}"}
    \item{addedcomment}{character string printed if \code{printt} after the title of the analysis}
   \item{chi2}{print the chi2 information when computing margins in \code{\link{FCAmet}}}
    \item{E}{ if not \code{NULL} is an array with the same dimensions as X}
    }

\details{
Gives the SVD-\emph{k}modes decomposition of the \eqn{1+\chi^2/N} of
the multiple contingency table of full count \eqn{N=\sum X <- {ijk...}},
i.e. complete independence + lack of independence (including marginal
independences) as shown for example in Lancaster(1951)(see reference
in Leibovici(2000)). Noting \eqn{P=X/N},  a \code{PTAk} of the
\eqn{(k+1)}-uple is done, e.g. for a three way contingency table
\eqn{k=3} the \emph{4}-uple data and metrics is:
 \deqn{ ((D_I^{-1} \otimes D_J^{-1} \otimes D_K^{-1})P, \quad D_I, \quad D_J, \quad D_K)}
where the metrics are diagonals of the corresponding margins. For
full description of arguments see \code{\link{PTAk}}. If \code{E}
is not \code{NULL} an FCAk-modes relatively to a model is
done (see Escoufier(1985) and therin reference
Escofier(1984) for a 2-way derivation, e.g. for a three way contingency table
\eqn{k=3} the \emph{4}-tuple data and metrics is:
 \deqn{ ((D_I^{-1} \otimes D_J^{-1} \otimes D_K^{-1})(P-E), \quad D_I, \quad D_J, \quad D_K)}
  If \code{E} was the complete independence (product of the margins)
  then this would give an \code{AFCk} but without looking at the
  marginal  dependencies (i.e. for a three way table no two-ways lack of
  independence are looked for).
}
 \value{a \code{solutions.FCAk} (inherits \code{\link{solutions.PTAk}}) object }
\references{
Escoufier Y (1985) \emph{L'Analyse des correspondances : ses
propriétés et ses extensions}. ISI 45th session Amsterdam.

Leibovici D (1993) \emph{Facteurs à Mesures Répétées et Analyses
Factorielles : applications à un suivi épidémiologique}. Université
de Montpellier II PhD Thesis in Mathématiques et Applications
(Biostatistiques).

Leibovici D (2000) \emph{Multiway Multidimensional Analysis for
Pharmaco-EEG Studies}. (submitted)
\url{http://www.fmrib.ox.ac.uk/~didier/cv/recentpub.html}
}
\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
\seealso{\code{\link{PTAk}}, \code{\link{FCAmet}}, \code{\link{summary.solutions.FCAk}}  }

\examples{
# try the demo
demo.FCAk()
}
\keyword{models}
\keyword{array  algebra}
\keyword{multivariate}
