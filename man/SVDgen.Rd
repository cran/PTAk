\name{SVDgen}
\alias{SVDgen}
\title{SVD with metrics and smoothing approximation}
\description{
Computes the generalised Singular Value Decomposition, \emph{i.e.} with
non-identity metrics. A smooth approximation can be asked to constraint the
components (\code{u} and \code{v}) to be smooth.
 }

\usage{
  SVDgen(Y,D2=1,D1=1,smoothing=FALSE,nomb=dim(Y)[2],
                  smoo=list(function(u)ksmooth(1:length(u),u,kernel="normal",
                           bandwidth=3,x.points=(1:length(u)))$y)) }
 \arguments{
   \item{Y}{a matrix \eqn{n \times p}}
   \item{D2}{metric in \eqn{R^p} either a vector (\eqn{p \times 1}) or a matrix (\eqn{p \times p})}
   \item{D1}{metric in \eqn{R^n} either a vector (\eqn{n \times 1}) or a matrix (\eqn{n \times n})}
   \item{smoothing}{logical if \code{TRUE}  the smoothing methods in \code{smoo} are used.}
   \item{smoo}{
list of lists of smoothing functions on a vector of the approriate dimension; if on one dimension it is
\code{NA} no smoothing will be done for this one; if the length of a list is one the function is used
for all components. If only one list in the list it will be used for both dimensions. } }

 \details{
The function computes the decomposition \eqn{X=UL^{1/2}V'} where \eqn{U'D_1U=Id_p} and
\eqn{V'D_2V=Id_p} and the diagonal matrix \eqn{L} containing no zeros squared singular values. If
\code{smoothing} a \emph{constraint}  on Least Squares solution  is used, then the
above decomposition becomes an approximation (unless \code{X} belongs to the space defined by the constraints). A \emph{Power Method} algorithm to compute each
principal tensor is used wherein Alternated Least Squares are always followed by a \emph{smoothed
version} of the updated vectors. If a Spline smoothing was used the algorithm would be equivalent
to use the traditional \emph{penalised least squares} at each iteration and could be called
\emph{Penalised Power Method} or Splined Alternated Least Squares Algorithm (SALSA is already  an
acronym used by Besse and Ferraty (1995) in where  a similar idea is developped: but smoothing
operates only on variables, and is\emph{ model based} as the Alternating operates on te whole
approximation \emph{i.e.} given the choice of the dimension reduction). }
\value{a \code{solutions.PTAk} object}

\note{\code{SVDgen} makes use of a non-identity version \code{svd} (inbuilt) or
\code{svdksmooth} which outputs like the inbuilt \code{svd}. The smoothing
option is also implemented in PTA-kmodes, FCA-kmodes, PCAn and
CANDECOMP/PARAFAC. The use of metrics (diagonal or not) allows flexibility of
analysis like \emph{e.g.} correspondence analysis, discriminant analysis,
robust analysis. Smoothing option extends the analysis towards functional data
analysis, and or outliers protection. It is theoretically valid for Principal
Tensors (here order 2) belonging to a tensor product of separable Hilbert
spaces (\emph{e.g.} Sobolev spaces) see Leibovici and El Maach (1997). The
function offers the choice to change of smoothing (method or parameters) as the
number of components grows as in Ramsay and Silverman (1997).}
 \references{
Leibovici D and  El Maache H (1997) \emph{Une décomposition en Valeurs Singulières d'un élément
d'un produit  Tensoriel de k espaces de Hilbert Séparables.} Compte Rendus de l'Académie des
Sciences tome  325, série I, Statistiques (Statistics) & Probabilités (Probability Theory):
779-782.

Besse P and Ferraty F (1995) \emph{Curvilinear fixed effect model}. Computational Statistics, 10:339-351.

Leibovici  D (2001) \emph{Metric choices for fMRI Multiway Data
Analysis}. (to be submitted)

Leibovici  D (2001) \emph{A Penalised algorithm for SVD and Multiway
   functional  methods}. (to be submitted)

Ramsay J.O. and Silverman B.W.  (1997) \emph{Functional Data
    Analysis}. Springer Series in Statistics.
}



\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
\seealso{\code{\link{PTAk}},\code{\link{PCAn}}, \code{\link{CANDPARA}} }

\examples{
library(modreg)
 library(tensor)

 # on smoothing

 data(longley)
 long <- as.matrix(longley[,1:7])

 long.svd <- SVDgen(long,smoothing=FALSE)
  summary.solutions.PTAk(long.svd,testvar=0)
   # X11(width=4,height=4)
  plot.solutions.PTAk(long.svd,scree=TRUE,RiskJack=0,type="b",lty=3)

 long.svdo <- SVDgen(long,smoothing=TRUE,
  smoo=list(function(u)ksmooth(1:length(u),
      u,kernel="normal",bandwidth=3,x.points=(1:length(u)))$y,NA))
  summary.solutions.PTAk(long.svdo,testvar=0)
  #  X11(width=4,height=4)
  plot.solutions.PTAk(long.svdo,scree=TRUE,RiskJack=0,type="b",lty=3)
     ####
comtoplot <- function(com=1,solua=long.svd,solub=long.svdo,openX11s=FALSE,...)
         {
  if(openX11s)X11(width=4,height=4)
 yla <- paste(solua[[2]]$vsnam[com],round((100*(solua[[2]]$d[com])^2)/
     solua[[2]]$ssX[1],2),"",
     round((100*(solub[[2]]$d[com])^2)/solua[[2]]$ssX[1],2))

limi <- range(c(solua[[1]]$v[com,],solub[[1]]$v[com,]))
  plot(solua,nb1=com, mod=1,type="b",lty=3,lengthlabels=4,cex=0.4,
   ylimit=limi,ylab=yla,col=2,...)
 par(new=TRUE)

  plot.solutions.PTAk(solub,nb1=com,mod=1,labels=F,type="b",lty=1,
  lengthlabels=4,cex=0.6,ylimit=limi,xylab=FALSE,...)
  par(new=FALSE)
}   ####
 comtoplot(com=1)
#  on using non-diagonal metrics

 data(crimerate)
  crimerate.mat <- sweep(crimerate,2,apply(crimerate,2,mean))
  crimerate.mat <- sweep(crimerate.mat,2,sqrt(apply(crimerate.mat,2,var)),FUN="/")
   metW <- Powmat(CauRuimet(crimerate.mat),(-1))
   # inverse of the within "group" (to play a bit more you could set m0 relating
   # the neighbourhood of states (see CauRuimet)

  cri.svd <- SVDgen(crimerate.mat,D2=1,D1=1)
  summary(cri.svd,testvar=0)
   plot(cri.svd,scree=TRUE,RiskJack=0,type="b",lty=3)
  cri.svdo <- SVDgen(crimerate.mat,D2=metW,D1=1)
   summary(cri.svdo,testvar=0)
   plot(cri.svdo,scree=TRUE,RiskJack=0,type="b",lty=3)
  # X11(width=8,height=4)
  par(mfrow=c(1,2))
   plot(cri.svd,nb1=1,nb2=2,mod=1,lengthlabels=3)
  plot(cri.svd,nb1=1,nb2=2,mod=2,lengthlabels=4,main="canonical")
  # X11(width=8,height=4)
  par(mfrow=c(1,2))
 plot(cri.svdo,nb1=1,nb2=2,mod=1,lengthlabels=3)
 plot(cri.svdo,nb1=1,nb2=2,mod=2,lengthlabels=4,
       main=expression(paste("metric ",Wg^{-1})))
# try the demo
demo.svdsmoo1()

    }

\keyword{smooth}
\keyword{multivariate}
