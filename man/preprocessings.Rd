\name{preprocessings}
\alias{preprocessings}

\title{Few useful functions for preprocessing arrays }
\description{Choices of centering or detrending  and scaling are important preprocessings for multiway analysis.}
\usage{
Multcent(dat=X,bi=c(1,2),by=3,
          centre=mean,
           centrebyBA=c(TRUE,FALSE),scalebyBA=c(TRUE,FALSE))
IterMV(n=10,dat=X,Mm=c(1,3),Vm=c(2,3),
               fFUN=mean,usetren=FALSE,
                  tren=function(x)smooth.spline(as.vector(x),df=5)$y,
                   rsd=TRUE)
Detren(dat,Mm=c(1,3),rsd=TRUE,
         tren=function(x)smooth.spline(as.vector(x),df=5)$y )
Susan1D(y,x=NULL,sigmak=NULL,sigmat=NULL,
          ker=list(function(u)return(exp(-0.5*u**2))))
}
\arguments{
  \item{---}{\emph{function} \code{Multcent}}
  \item{dat}{array}
  \item{bi}{vector defining the "centering, bicentering or multi-centering" one wants
    to operate crossed with \code{by}}
  \item{by}{number or vector defining the entries used "with" in
    the other operations}
  \item{centre}{function used as \code{FUN} in applying
    "multi-centering"}
  \item{centrebyBA}{a bolean vector for "centering" with \code{centre} Before and After
    according to \code{by}}
  \item{scalebyBA}{idem as centrebyBA, for scaling operation}
  \cr
  \item{---}{\emph{function} \code{IterMV}}
  \item{n}{number of iterations between "centering" and scaling}
  \item{dat}{array}
  \item{Mm}{margins to performs \code{Detren} or \code{fFUN} on}
  \item{Vm}{margins to scale}
  \item{fFUN}{function to use as \code{FUN} if \code{usetren} is
    \code{FALSE}}
  \item{usetren}{logical, to use \code{Detren}}
  \item{tren}{function to use in \code{Detren}}
  \item{rsd}{logical passed into \code{Detren} (only) to detrend or not}
  \cr
  \item{----}{\emph{function} \code{Detren}}
  \item{dat}{array}
  \item{Mm}{as above}
  \item{rsd}{as above}
  \item{tren}{the function to use}
  \cr
  \item{----}{\emph{function} \code{Susan1D}}
  \item{y}{vector (length \code{n})}
  \item{x}{vector of same length, if \code{NULL} it is \code{1:n} }
  \item{sigmak}{parameter related to kernel bandwidth with \code{y}
    values (default is \code{1/2*range}}
  \item{sigmat}{parameter related to kernel bandwidth with  \code{x}
    values (default value is \code{8*n^{-1/5}}, with a minimum number of
  neigbours set as one apart)}
  \item{ker}{a list of two kernels \code{list("t"=function "k"=function
      ) } for each weightings (if only one given it is used for
      both)}
  }
\details{\code{Multcent} performs in order "centering" by \code{by};
  "multicentering" for every \code{bi} with \code{by}; then scale
  (standard deviation) to one by \code{by}.

   \code{IterMV} performs an iterative "detrending" and scaling
   according to te margins defined (see Leibovici(2000) and references
   in it).

    \code{Detren} detrends (or smooths if \code{rsd} is \code{FALSE})
    the data accoding to th margins given.

    \code{Susan1D} performs a non-linear kernel smoothing  of \code{y}
    against \code{x} (both reordered in the function according to orders
     of \code{x}) with an usual kernel (\code{t}) as for kernel
     regression and a kernel (\code{t}) for the values of \code{y} (the
     product of the kernels constitutes the non-linear weightings. This
     function is adapted from SUSAN algorithm (see references).



}
  \author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}

\references{
 Smith S.M. and J.M. Brady (1997) \emph{SUSAN - a new approach to low level image processing}.
     International Journal of Computer Vision, 23(1):45-78, May 1997.
   }

\keyword{multivariate}
\keyword{smooth}
