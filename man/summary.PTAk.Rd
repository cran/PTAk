\name{summary.PTAk}
\alias{summary.PTAk}
\alias{summary.FCAk}
\title{Summary of a PTA-k modes analysis}
\description{Print a summary listing of the decomposition obtained.
}
\usage{
\method{summary}{PTAk}(object,testvar=1,dontshow="*", \dots)
\method{summary}{FCAk}(object,testvar=0.5,dontshow="*", \dots)
}
\arguments{
   \item{object}{an object inheriting from class \code{PTAk}, representing a generalised  singular value decomposition}
   \item{testvar}{control within \code{nTens} used Principal Tensor with  minimum percent of variability explained}
   \item{dontshow}{
boolean criterion to remove Principal Tensors from the summary, or
default is a character "*" equivalent to the criterion:
\cr
\code{!substr(solution[[length(solution)]][["vsnam"]],1,1)=="*"} }
\item{...}{summary generic additional arguments not used here }
}
\details{
The function prints a listing of the decomposition with historical
order (instead of traditional singular value order). It is useful
before any plots or reconstruction, a screeplot (using
\code{plot.PTAk}) will be also useful. It is useful before any plots
r reconstruction, a screeplot (using \code{plot.PTAk}) will be also
useful. \code{summary.FCAk} is alike
\code{summary.PTAk} but \code{testvar} operates on the variability of
the lack of complete independence. }

\value{prints on the prompt with an invisible return of the summary table}
\references{
Leibovici D (2000) \emph{Multiway Multidimensional Analysis for
  Pharmaco-EEG Studies}.(submitted)
\url{http://c3s2i.free.fr/cv/recentpub.html}
}

\author{Didier G. Leibovici \email{GeotRYcs@gmail.com}}

\note{At the moment can be used for \code{PCAn},
  \code{CANDPRA}, better summaries will be in the next release.}
 \seealso{ \code{\link{plot.PTAk}} }

 \examples{
 data(crimerate)
 crimerate.mat <- sweep(crimerate,2,apply(crimerate,2,mean))
 crimerate.mat <- sweep(crimerate.mat,2,sqrt(apply(crimerate,2,var)),FUN="/")
 cri.svd <- SVDgen(crimerate.mat)
 summary(cri.svd,testvar=0)
  plot(cri.svd,scree=TRUE)
  par(new=TRUE)
  RiskJackplot(cri.svd,nbvs=1:7,mod=NULL,max=NULL,rescaled=TRUE,
        axes=FALSE,ann=FALSE)
  par(new=FALSE)

  # or equivalently

  plot(cri.svd,scree=TRUE,type="b",lty=3,RiskJack=1) #set mod=NULL or c(1,2)
  ###
   data(crimerate)
   criafc <- FCAmet(crimerate,chi2=TRUE)
   cri.afc <- SVDgen(criafc$data,criafc$met[[2]],criafc$met[[1]])
    summary(cri.afc)
   plot(cri.afc,scree=TRUE)
   plot(cri.afc,scree=TRUE,type="b",lty=3,RiskJack=1,method="FCA")


   }

\keyword{array }  \keyword{algebra}
\keyword{multivariate}
