\name{is.solutions.PTAk}
\alias{solutions.PTAk}
\title{Output object from a PTA-kmodes analysis and other
methods implemented in this package}
\description{
The object class  \code{solutions.PTAk} is list of list.
\code{solutions[[i]]} is a list with arguments the following
arguments: }

\usage{
is.solutions.PTAk(solutions)
}
\arguments{
   \item{v}{a matrix \code{nb of solutions} \eqn{\times } \code{dimension of the mode} \eqn{i} }
   \item{modesnam}{a character object naming the mode, \code{"m i"} otherwise}
   \item{n}{labels of mode \code{i} entries as given in \code{dimnames} of the data \cr }
   \item{-}{}
   \item{d}{singular values}
   \item{pct}{local percent of sum of square for the above singular values \eqn{100 d^2/ssX}}
   \item{ssX}{local sum of squares, i.e. for the tensor currently analysed (nested projections of \code{X},
    \code{ssX[1]} is the total sum of squares}
   \item{smoocheck}{bolean matrix with \emph{k} rows and \code{length(vsnam)} columns to check
   if a smooth componnet was computed}
   \item{vsnam}{names of the singular values, see details}
   \item{datanam}{name of the tensor analysed}
   \item{method}{the \code{\link{match.call}} of the function which gave the resulting solutions.PTAk object.}
   }

 \details{
The object \code{solutions.PTAk} is output from \code{PTAk} is a list
of length the order of the tensor used in the analysis. Each element
of the list is a list with the above arguments (the last arguments
starting from \code{d} are given only for the last member of the
list ). For a given row number \eqn{t} the tensor product of the mode
components \code{v[t,]} time the corresponding \code{d[t]}
constitute a Principal Tensor of \code{datanam}. The corresponding
name of the singular value follows the rule: \emph{e.g.} \code{vs1111},
\code{vs2222} are \emph{k}modes solutions (\emph{4}-modes analysis);
in a \emph{3}-modes analysis \code{12 vs111 5 12} means an associated
solution to the \emph{k}modes solutions (\code{vs111}) according to
the mode of dimension \eqn{12}, the other dimensions \eqn{5} and
\eqn{12} are given only when the \code{vs} results from a \emph{2}modes
analysis; associated of associated have a form like \code{9-15
vs1111 10 888} (from a \emph{4}-modes analysis).}
\value{
The class  object is a list of a lists with the specified
arguments and right dimensions. The class object is used in
\code{RESUM}, \code{summary}, \code{plot}, \code{REBUILD}, \code{APSOLU3} and \code{APSOLUk}.
   }
\note{The function \code{\link{SVDgen}} output an object of class \code{solutions.PTAk}.}
 \author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
\seealso{\code{RESUM}}, \code{\link{summary.solutions.PTAk}},
\code{\link{plot.solutions.PTAk}},\code{\link{REBUILD}}. }

\keyword{array  algebra}
\keyword{list}
 \keyword{classes}
