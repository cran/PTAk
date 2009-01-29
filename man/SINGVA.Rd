\name{SINGVA}
\alias{SINGVA}
\title{Optimisation algorithm RPVSCC}
\description{Computes the best rank-one approximation using the RPVSCC algorithm.}
\usage{
SINGVA(X,test=1E-12,PTnam="vs111",Maxiter=2000,
                  verbose=getOption("verbose"),file=NULL,
                    smoothing=FALSE,smoo=list(NA),
                     modesnam=NULL,
                      Ini="Presvd",sym=NULL)}
\arguments{
   \item{X}{
a tensor (as an array) of order \emph{k}, if non-identity metrics are
used \code{X} is a list with \code{data}  as the array and
\code{met} a list of metrics}
   \item{test}{numerical value to stop optimisation}
   \item{PTnam}{character giving the name of the \emph{k}-modes Principal Tensor}
   \item{Maxiter}{if \code{iter > Maxiter} prompts to carry on or not, then do it
     every other 200 iterations}
   \item{verbose}{control printing}
\item{file}{output printed at the prompt if \code{NULL}, or printed in the given  \file{file}}
   \item{smoothing}{logical to use smooth functiosns or not (see
     \code{\link{SVDgen}}) }
   \item{smoo}{list of functions returning smoothed vectors (see
     \code{\link{PTA3}}) }
   \item{modesnam}{character vector of the names of the modes, if \code{NULL} "\code{mo 1}" ..."\code{mo k}"}
   \item{Ini}{method used for initialisation of the algorithm (see \code{INITIA})}
   \item{sym}{description of the symmetry of the tensor \emph{e.g.} c(1,1,3,4,1) means the
     second mode and the fifth are identical to the first }
}

\details{
The algorithm termed \emph{RPVSCC}  in Leibovici(1993) is implemented
to compute the first Principal Tensor (rank-one tensor with its
singular value) of the given tensor \code{X}. According to the
decomposition described in Leibovici(1993) and Leibovici and
Sabatier(1998), the function gives a generalisation  to \emph{k}
modes of the \emph{best rank-one approximation} issued from SVD whith
2 modes. It is identical to  the PCA-\emph{k}modes if only 1
dimension is asked in each space, and to PARAFAC/CANDECOMP if the
rank of the approximation is fixed to 1. Then the methods differs,
PTA-\emph{k}modes will look for best approximation according to the
\emph{orthogonal rank} (\emph{i.e.} the rank-one tensors (of the
decomposition) are orthogonal), PCA-\emph{k}modes will look for best
approximation according to the \emph{space ranks} (\emph{i.e.}  ranks
of every bilinear form deducted from the original tensor, that is the
number of components in each space), PARAFAC/CANDECOMP will look for
best approximation according to the \emph{rank} (\emph{i.e.} the
rank-one tensors are not necessarily orthogonal).
\cr
Recent work from Tamara G Kolda showed on an example that \emph{orthogonal rank}
decompositions are not necesseraly nested. This makes PTA-\emph{k}modes a model with
nested decompositions not giving the exact \emph{orthogonal rank}.
 So PTA-\emph{k}modes will look for best approximation according to orthogonal tensors in a nested approximmation process.
 }
\value{a \code{\link{PTAk}} object (without \code{datanam method})}
\references{
Kroonenberg P (1983) \emph{Three-mode Principal Component Analysis:
Theory and Applications}. DSWO press. Leiden.(related references in
\url{http://three-mode.leidenuniv.nl})
                         
Leibovici D(1993) \emph{Facteurs <e0> Mesures R<e9>p<e9>t<e9>es et Analyses Factorielles :
applications <e0> un suivi <e9>pid<e9>miologique}. Universit<e9> de Montpellier
II. PhD Thesis in Math<e9>matiques et Applications (Biostatistiques).


Leibovici D and Sabatier R (1998) \emph{A Singular Value
Decomposition of a k-ways array for a Principal Component Analysis of
multi-way data, the PTA-k}. Linear Algebra and its Applications,
269:307-329.

De Lathauwer L, De Moor B and Vandewalle J (2000) \emph{On the best
rank-1 and rank-(R1,R2,...,Rn) approximation of higher-order
tensors}. SIAM J. Matrix Anal. Appl. 21,4:1324-1342.


Kolda T.G (2003)\emph{ A Counterexample to the Possibility of an Extension of the Eckart-Young Low-Rank Approximation Theorem for the Orthogonal Rank Tensor Decomposition}. SIAM J. Matrix Analysis, 24(2):763-767, Jan. 2003.
}
\note{
The algorithm was derived in generalising the \emph{transition
formulae} of SVD (Leibovici 1993), can also be understood as a
generalisation of the \emph{power method} (De Lathauwer et al.
2000). In this paper they also use a similar algorithm  to build
bases in each space, reminiscent of three-modes and \emph{n}-modes
PCA (Kroonenberg(1980)), \emph{i.e.} defining what they called a
rank-(R1,R2,...,Rn) approximation (called here \emph{space ranks},
see \code{PCAn}). \emph{RPVSCC} stands for  \emph{R}echerche de la \emph{P}remi<e8>re
\emph{V}aleur \emph{S}inguli<e8>re par \emph{C}ontraction
\emph{C}ompl<ea>te.
 }
\author{Didier Leibovici \email{c3s2i@free.fr}}
\seealso{\code{\link{INITIA}}, \code{\link{PTAk}}, \code{\link{PCAn}},
\code{\link{CANDPARA}} }

\keyword{smooth}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
