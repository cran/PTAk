\name{PTAk-internal}

\alias{Ginv}
\alias{PPMA}
\alias{Powmat}
\alias{RaoProd}
\alias{REBUILDPCAn}
\alias{RESUM}
\alias{svdsmooth}
\alias{toplist}

\title{Internal PTAk functions}
\description{
Internal PTAk functions }
\usage{
Ginv(A)
PPMA(X,test=1E-10,pena=list(function(u)ksmooth(1:length(u),u,kernel="normal",
                     bandwidth=3,x.points=(1:length(u)))$y
                     ,NA) ,ini=mean,vsmin=1E-20,Maxiter=2000)
Powmat(A,pw,eltw=FALSE)
RaoProd(A,B)
REBUILDPCAn(solu)
RESUM(solb,sola=NULL,numass=NULL,verbose=getOption("verbose"),file=NULL
                        ,summary=FALSE,testvar=0.1,not=NULL)
svdsmooth(X,nomb=min(dim(X)),
                   smooth=list(function(u)ksmooth(1:length(u),u,kernel="normal",
                               bandwidth=3,x.points=(1:length(u)))$y),vsmin=1E-16)
toplist(li)
}

\details{
These functions are not supposed to be called directly. }

 \author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
  \seealso{\code{\link{PTAk}}}
\keyword{misc}
