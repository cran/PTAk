\name{summary.solutions.PTAk}
\alias{summary.solutions.PTAk}
\alias{summary.solutions.FCAk}
\title{Summary of a PTA-k modes analysis}
\description{Print a summary listing of the decomposition obtained.
}
\usage{
summary.solutions.PTAk(object,testvar=1,dontshow="*")
summary.solutions.FCAk(object,testvar=1,dontshow="*")
summary(object,testvar=1,dontshow="*")
}
\arguments{
   \item{object}{a \code{solutions.PTAk} object}
   \item{testvar}{control within \code{nTens} used Principal Tensor with  minimum percent of variability explained}
   \item{dontshow}{
boolean criterion to remove Principal Tensors from the summary, or
default is a character "*" equivalent to the criterion:
\cr
\code{!substr(object[[length(object)]][["vsnam"]],1,1)=="*"} } }
\details{
The function prints a listing of the decomposition with historical
order (instead of traditional singular value order). It is useful
before any plots or reconstruction, a screeplot (using
\code{plot.PTAk}) will be also useful. It is useful before any plots
r reconstruction, a screeplot (using \code{plot.solutions.PTAk}) will be also
useful. \code{summary.solutions.FCAk} is alike
\code{summary.solutions.PTAk} but \code{testvar} operates on the variability of
the lack of complete independence. }

\value{prints on the prompt}
\references{
Leibovici D (2000) \emph{Multiway Multidimensional Analysis for
  Pharmaco-EEG Studies}.(submitted)
\url{http://www.fmrib.ox.ac.uk/~didier/cv/recentpub.html}
}

\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}

\note{At the moment can be used for \code{solutions.PCAn},
  \code{solutions.CANDPRA}, better summaries will be in the next release.}
 \seealso{ \code{\link{plot.solutions.PTAk}} }

 \examples{
 data(crimerate)
 crimerate.mat <- sweep(crimerate,2,apply(crimerate,2,mean))
 crimerate.mat <- sweep(crimerate.mat,2,sqrt(apply(crimerate,2,var)),FUN="/")
 cri.svd <- SVDgen(crimerate.mat)
 summary(cri.svd,testvar=0)
  plot(cri.svd,scree=TRUE)
  par(new=TRUE)
  RiskJack.plot(cri.svd,nbvs=1:7,mod=NULL,max=NULL,rescaled=TRUE,
        axes=FALSE,ann=FALSE)
  par(new=FALSE)

  # or equivalently

  plot(cri.svd,scree=TRUE,type="b",lty=3,RiskJack=1) #set mod=NULL or c(1,2)
  ###
demo.FCA <- function(){
   data(crimerate)
   criafc <- FCAmet(crimerate,chi2=TRUE)
   cri.afc <- SVDgen(criafc$data,criafc$met[[2]],criafc$met[[1]])
    summary(cri.afc)
   plot(cri.afc,scree=TRUE)
   plot(cri.afc,scree=TRUE,type="b",lty=3,RiskJack=1,method="FCA")
   }

   }

\keyword{array  algebra}
\keyword{multivariate}
