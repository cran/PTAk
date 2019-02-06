\name{PTA3}
\alias{PTA3}
\title{Principal Tensor Analysis on 3 modes}
\description{
Performs a truncated SVD-\emph{3}modes analysis with or
without specific metrics, penalised or not.
 }
 \usage{PTA3(X,nbPT=2,nbPT2=1,
                smoothing=FALSE,
                smoo=list(function(u)ksmooth(1:length(u),u,kernel="normal",
                          bandwidth=4,x.points=(1:length(u)))$y,
                          function(u)smooth.spline(u,df=3)$y,
                          NA),
                minpct=0.1,verbose=getOption("verbose"),file=NULL,
                modesnam=NULL,addedcomment="", ...)}

\arguments{
   \item{X}{
a tensor (as an array) of order \emph{3}, if non-identity metrics are
used \code{X} is a list with \code{data}  as the array and
\code{met} a list of metrics}
   \item{nbPT}{a number specifying the number of \emph{3}modes Principal Tensors requested}
   \item{nbPT2}{if 0 no \emph{2}-modes solutions will be computed, 1 =all, >1 otherwise}
   \item{smoothing}{logical to consider smoothing or not}
   \item{smoo}{a list of length 3 with lists of functions operating on
     vectors component for the appropriate dimension (see details) }
   \item{minpct}{
numerical 0-100 to control of computation of future solutions at this
level and below}
\item{verbose}{control printing}
 \item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
   \item{modesnam}{
character vector of the names of the modes, if \code{NULL} "\code{mo 1}"
..."\code{mo k}"}
 \item{addedcomment}{character string printed after the title of the analysis}
\item{...}{any other arguments passed to SVDGen or other functions}

   }

 \details{
According to the decomposition described in Leibovici(1993) and
Leibovici and Sabatier(1998) the function gives a generalisation of
the SVD (2 modes) to \emph{3} modes. It is the same algorithm used
for \code{PTAk} but simpler as the recursivity implied by the
\emph{k} modes analysis is reduced only to one level \emph{i.e} for
every \emph{3}-modes  Principal Tensors, 3 SVD are performed for
every contracted product with one the three components of the
\emph{3}-modes  Principal Tensors (see \code{\link{APSOLU3}},
\code{\link{PTAk}}). \cr

Recent work from Tamara G Kolda showed on an example that \emph{orthogonal rank}
decompositions are not necesseraly nested. This makes PTA-\emph{3}modes a model with
nested decompositions not giving the exact \emph{orthogonal rank}.
 So PTA-\emph{3}modes will look for best approximation according to orthogonal tensors in a nested approximmation process. PTA3 decompositions is "a" generalisation of SVD but not the ...
 \cr
With the \code{smoothing} option \code{smoo} contain a list of (lists) of functions to
apply on vectors of component (within the algorithm, see
\code{\link{SVDgen}}). For a given dimension (1,2,or 3) a list of
functions is given. If this list consists only of one function (no list
needed) this
function will be used at any level all the time : if one want to smooth
only for the first Principal Tensor, put \code{list(function,NA)}. Now
you start to understand this list will have a maximum length of
\code{nbPT} and the corresponding function will be used for the
corresponding \emph{3}mode Principal Tensor. To smooth differently the
associated solutions one have to put another level of nested lists
otherwise the function given at the \emph{3}mode level will be used for
all. These rules are te same for \code{\link{PTAk}}.
}

\value{a \code{\link{PTAk}} object}
\references{
Leibovici D(1993) \emph{Facteurs <e0> Mesures R<e9>p<e9>t<e9>es et Analyses Factorielles :
applications <e0> un suivi <e9>pid<e9>miologique}. Universit<e9> de Montpellier
II. PhD Thesis in Math<e9>matiques et Applications (Biostatistiques).

Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component Analysis of
multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329.

Kolda T.G (2003)\emph{ A Counterexample to the Possibility of an Extension of the Eckart-Young Low-Rank Approximation Theorem for the Orthogonal Rank Tensor Decomposition}. SIAM J. Matrix Analysis, 24(2):763-767, Jan. 2003.
 }

  \note{
The use of metrics (diagonal or not) allows flexibility of analysis like in \emph{2} modes
\emph{e.g.} correspondence analysis, discriminant analysis, robust analysis.  Smoothing option
extends the analysis towards functional data analysis, and or outliers "protection" is
theoretically valid for tensors belonging to a tensor product of separable Hilbert spaces
(\emph{e.g.} Sobolev spaces) (see references in \code{PTAk}, \code{SVDgen}).}

 \author{Didier G. Leibovici}
\seealso{\code{\link{SVDgen}}, \code{\link{FCAk}}, \code{\link{PTAk}}, \code{\link{summary.PTAk}}}
 \examples{
 # example using Zone_climTUN dataset
 #  
# library(maptools)
# library(RColorBrewer)
# Yl=brewer.pal(11,"PuOr")
# data(Zone_climTUN)
## in fact a modified version of plot.Map was used
# plot(Zone_climTUN,ol=NA,auxvar=Zone_climTUN$att.data$PREC_OCTO)
##indicators 84 +3 to repeat
# Zone_clim<-Zone_climTUN$att.data[,c(2:13,15:26,28:39,42:53,57:80,83:95,55:56)]
# Zot <-Zone_clim[,85:87] ;temp <-colnames(Zot) 
# Zot <- as.matrix(Zot)%x%t(as.matrix(rep(1,12)))
# colnames(Zot) <-c(paste(rep(temp [1],12),1:12),paste(rep(temp [2],12),1:12),
#	 paste(rep(temp [3],12),1:12))
# Zone_clim <-cbind(Zone_clim[,1:84],Zot)
 
# Zone3w <- array(as.vector(as.matrix(Zone_clim)),c(2599,12,10))
## preprocessing
#Zone3w<-Multcent(dat=Zone3w,bi=NULL,by=3,centre=mean, 
#	centrebyBA=c(TRUE,FALSE),scalebyBA=c(TRUE,FALSE))
# Zone3w.PTA3<-PTA3(Zone3w,nbPT=3,nbPT2=3)
## summary and plot
# summary(Zone3w.PTA3)
#plot(Zone3w.PTA3,mod=c(2,3),nb1=1,nb2=11,lengthlabels=5,coefi=list(c(1,1,1),c(1,-1,-1)))
#plot(Zone_climTUN,ol=NA,auxvar=Zone3w.PTA3[[1]]$v[1,],nclass=30)
#plot(Zone_climTUN,ol=NA,auxvar=Zone3w.PTA3[[1]]$v[11,],nclass=30)
 
 ##############
 cat(" A little fun using iris3 and matching randomly 15 for each iris sample!","\n")
 cat("   then performing a PTA-3modes.  If many draws are done, plots")
 cat("   show the stability of the first and third Principal Tensors.","\n")
 cat("iris3 is centered and reduced beforehand for each original variables.","\n")
#  demo function 
# source(paste(R.home(),"/library/PTAk/demo/PTA3.R",sep=""))
# demo.PTA3(bootn=10,show=5,openX11s=FALSE)
 
}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
