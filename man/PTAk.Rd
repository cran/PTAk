\name{PTAk}
\alias{PTAk}
\title{Principal Tensor Analysis on k modes}
\description{
Performs a truncated SVD-\emph{k}modes analysis with or without
specific metrics, penalised or not. }
\usage{
 PTAk(X,nbPT=2,nbPT2=1,minpct=0.1,
                smoothing=FALSE,
                   smoo=list(NA),
                    verbose=getOption("verbose"),file=NULL,
                     modesnam=NULL,addedcomment="", ...)}

\arguments{
   \item{X}{
a tensor (as an array) of order \emph{k}, if non-identity metrics are
used \code{X} is a list with \code{data}  as the array and
\code{met} a list of metrics}
   \item{nbPT}{
integer vector of length \emph{(k-2)}  specifying the  maximum number
of Principal Tensors requested for the \emph{(3,...,k-1, k)} modes
levels (see details), if it is not a vector every levels would have
the same given nbPT value }
   \item{nbPT2}{if 0 no 2-modes solutions
     will be computed, 1 =all, >1 otherwise}
   \item{minpct}{
numerical 0-100 to control of computation of future solutions at this
level and below}
   \item{smoothing}{see \code{\link{PTA3}}, \code{\link{SVDgen}}}
   \item{smoo}{see \code{\link{PTA3}}}
\item{verbose}{control printing}
 \item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
   \item{modesnam}{character vector of the names of the modes, if \code{NULL} \code{mo 1} ...\code{mo k}}
   \item{addedcomment}{character string printed if \code{printt} after the title of the analysis}
\item{...}{any other arguments passed to other functions}
   }

 \details{
According to the decomposition described in Leibovici(1993) and
Leibovici and Sabatier(1998) the function gives a generalisation of
the SVD (2 modes) to \emph{k} modes. The algorithm is recursive,
calling \code{APSOLUk} which calls \code{PTAk} for \emph{(k-1)}.
\code{nbPT}, \code{nbPT2} and \code{minpct} control the number of
Principal Tensors desired. For example \code{nbPT=c(2,4,3)} means a
tensor of order \emph{5} is analysed, the maximum number of 5-modes
PT is set to \emph{3}, for \emph{each of them} one sets a maximum of
\emph{4} associated \emph{4-modes} (for each of the five components),
for \emph{each of these later} a maximum of \emph{2} associated
\emph{3-modes} PT is asked (for each of the four components). Then
\code{nbPT2} complete for \emph{2-modes} associated or not. Overall
\code{minpct} controls to carry on the algorithm at any level and
lower, \emph{i.e.} stops if \eqn{100(vs^2/ssx)<minpct} (where
\eqn{vs} is the singular value, and \emph{ssx} is the total sum of
squares of the tensor \eqn{X} or the "metric transformed" \eqn{X}).
Putting a \eqn{0} at a given level in \code{nbPT} obviously
automatically puts \eqn{0}  in \code{nbPT} at lower levels. Putting
high values in \code{nbPT} allows control only on \code{minpct}
helping to reach the full decomposition. All these controls allow to
truncate the full decomposition in a level-controlled fashion. Notice
the full decomposition always contains any possible choice of
truncation, \emph{i.e.} the solutions are not dependant on the
truncation scheme (Generalised Eckart-Young Theorem).
\cr
Recent work from Tamara G Kolda showed on an example that \emph{orthogonal rank}
decompositions are not necesseraly nested. This makes PTA-\emph{k}modes a model with
nested decompositions not giving the exact \emph{orthogonal rank}.
 So PTA-\emph{k}modes will look for best approximation according to orthogonal tensors in a nested approximmation process.}

  \value{a \code{\link{PTAk}} object which consist of a list of lists. Each mode has a list in which is listed:
\item{$v}{matrix of components for the given mode}
\item{$iter}{vector of iterations numbers where maximum was reach}
\item{$test}{vector of test values at maximum}
\item{$modesnam}{name of the mode}
\item{$v}{matrix of components for the given mode}
 	The last mode list has also some additional information on the analysis done:
\item{$d}{vector of singular values}
\item{$pct}{percentage of sum of squares for each quared singular value}
\item{$ssX}{vector of local sum of squares \emph{i.e.} of the current tensor with the rescursive algorithm}
\item{$vsnam}{vector of names given to the singular value according to a recursive data dependent scheme}
\item{$datanam}{data reference}
\item{$method}{call applied: could be PTAk or CANDPARA or PCAn or even SVDgen, with parameters choices}
\item{$addedcomment}{ the addedcomment (repeated) given in the call}
You will notice that methods other than PTAk may not have all list elements but the essential ones such as: $v, $d, $ssX, and may also have additional ones like $coremat for PCAn (the core array).
}

\author{Didier Leibovici \email{c3s2i@free.fr}}
\references{
Leibovici D(1993) \emph{Facteurs <e0> Mesures R<e9>p<e9>t<e9>es et Analyses Factorielles :
applications <e0> un suivi <e9>pid<e9>miologique}. Universit<e9> de Montpellier
II. PhD Thesis in Math<e9>matiques et Applications (Biostatistiques).

Leibovici D and El Maache H (1997) \emph{Une d<e9>composition en Valeurs Singuli<e8>res d'un <e9>l<e9>ment
d'un produit Tensoriel de k espaces de Hilbert S<e9>parables}. Compte Rendus de l'Acad<e9>mie des
Sciences tome 325, s<e9>rie I, Statistiques (Statistics) & Probabilit<e9>s (Probability Theory):
779-782.

Leibovici D and Sabatier R (1998) \emph{A Singular Value Decomposition of a k-ways array for a
Principal Component Analysis of multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329.

Leibovici D (2008) \emph{Spatio-temporal Multiway Decomposition using Principal Tensor Analysis on k-modes:the R package \pkg{PTAk}}
. to be submitted soon at Journal of Statisticcal Software.

Leibovici  D (2008) \emph{A Simple Penalised algorithm for SVD and Multiway
   functional  methods}. (to be submitted)

Kolda T.G (2003)\emph{ A Counterexample to the Possibility of an Extension of the Eckart-Young Low-Rank Approximation Theorem for the Orthogonal Rank Tensor Decomposition}. SIAM J. Matrix Analysis, 24(2):763-767, Jan. 2003.
}
\note{
The use of metrics (diagonal or not) allows flexibility of analysis like in \emph{2} modes
\emph{e.g.} correspondence analysis, discriminant analysis, robust analysis.  Smoothing option
extending the analysis towards functional data analysis is theoretically
valid for Principal Tensors  belonging to a tensor product of separable
Hilbert spaces (\emph{e.g.} Sobolev spaces) see Leibovici and El Maach (1997).}
 \seealso{\code{\link{REBUILD}}, \code{\link{FCAk}}, \code{\link{PTA3}} \code{\link{summary.PTAk}}
      }
 \examples{

  # don <- array((1:3)%x%rnorm(6*4)%x%(1:10),c(10,4,6,3))

don <- array(1:360,c(5,4,6,3))
 don <- don + rnorm(360,1,2)

 dimnames(don) <- list(paste("s",1:5,sep=""),paste("T",1:4,sep=""),
          paste("t",1:6,sep=""),c("young","normal","old"))
   # hypothetic data on learning curve at different age and period of year

 ones <-list(list(v=rep(1,5)),list(v=rep(1,4)),list(v=rep(1,6)),list(v=rep(1,3)))

 don <- PROJOT(don,ones)
 don.sol <- PTAk(don,nbPT=1,nbPT2=2,minpct=0.01,
               verbose=TRUE,
                modesnam=c("Subjects","Trimester","Time","Age"),
                 addedcomment="centered on each mode")

don.sol[[1]] # mode Subjects results and components
don.sol[[2]] # mode Trimester results and components
don.sol[[3]] # mode Time results and components
don.sol[[4]] # mode Age results and  components with additional information on the call

 summary(don.sol,testvar=2)
  plot(don.sol,mod=c(1,2,3,4),nb1=1,nb2=NULL,
     xlab="Subjects/Trimester/Time/Age",main="Best rank-one approx" )
  plot(don.sol,mod=c(1,2,3,4),nb1=4,nb2=NULL,
      xlab="Subjects/Trimester/Time/Age",main="Associated to Subject vs1111")

#  demo function
 # demo.PTAk()

  }
 \keyword{array }  \keyword{algebra}
 \keyword{multivariate}
