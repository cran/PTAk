\name{PTAk-internal}

\alias{Ginv}
\alias{PPMA}
\alias{Powmat}
\alias{RaoProd}
\alias{REBUILDPCAn}
\alias{RESUM}
\alias{svdsmooth}
\alias{toplist}
\alias{svd.p}
\title{Internal PTAk functions}
\description{
Internal PTAk functions }
\usage{
Ginv(A)
PPMA(X,test=1E-10,pena=list(function(u)ksmooth(1:length(u),u,kernel="normal",
                     bandwidth=3,x.points=(1:length(u)))$y
                     ,NA) ,ini=mean,vsmin=1E-20,Maxiter=2000, ...)
Powmat(A,pw,eltw=FALSE)
RaoProd(A,B)
REBUILDPCAn(solu)
RESUM(solb,sola=NULL,numass=NULL,verbose=getOption("verbose"),file=NULL
                        ,summary=FALSE,testvar=0.1,with=TRUE)
svdsmooth(X,nomb=min(dim(X)),
                   smooth=list(function(u)ksmooth(1:length(u),u,kernel="normal",
                               bandwidth=3,x.points=(1:length(u)))$y),vsmin=1E-16, ...)
toplist(li)
svd.p(X,...)
}
\arguments{
These functions are not supposed to be called directly.
        \item{X}{a matrix}
         \item{test}{a zero limit number}
         \item{pena}{ list of functions to be used as smoother}
         \item{ini}{initialisation method over the dual dimension}
         \item{vsmin}{zero limit for singular value}
          \item{Maxiter}{limit number of iteration}
           \item{A}{a matrix}
         \item{pw}{ power value number}
          \item{eltw}{ boolean to perform power elementwise or matrix power }
            \item{B}{a matrix}
              \item{solb}{an object inheriting from class \code{PTAk} }
                \item{sola}{an object inheriting from class \code{PTAk} }
                 \item{solu}{an object inheriting from class \code{PTAk} }
                  \item{numass}{position number of the associated solution, NULL is equivalent to the last in \code{sola}}
                  \item{verbose}{boolean playing a verbose role}
                     \item{file}{string pointing a destination of  file output}
                  \item{summary}{boolean to show the summary or not }
                      \item{testvar}{threshold control for minimum percent of variability explained}
                       \item{with}{boolean expression to give a supplementary selection criterion}
                 \item{nomb}{integer giving the number of components to fit}
                 \item{smooth}{idem as pena}
                 \item{li}{any list}
\item{...}{any other arguments passed to functions}
                 }
 \author{Didier G. Leibovici}
  \seealso{\code{\link{PTAk}}}
\keyword{misc}
