\name{plot.PTAk}
\alias{plot.PTAk}
\alias{RiskJackplot}
\title{Plot a PTAk object}
\description{
Screeplot of singular values or superposed plot of modes for one or two components (1 dimensional scatterplot with spread labels or scatterplot on two dimensions).}
\usage{
\method{plot}{PTAk}(x, labels = TRUE, mod = 1, nb1 = 1, nb2 = NULL,
    coefi = list(NULL, NULL), xylab = TRUE, ppch = (1:length(solution)),
    lengthlabels = 2, scree = FALSE, ordered = TRUE,
    nbvs = 40, RiskJack = NULL, method = "",ZoomInOut=NULL, Zlabels=NULL, Zcol=NULL,
    poslab=c(2,1,3,3), signedCTR = FALSE, relCTR = TRUE,\dots)
RiskJackplot(x, nbvs = 1:20, mod = NULL, max = NULL, rescaled=TRUE, \dots)
                              }
\arguments{
\item{x}{an object inheriting from class \code{PTAk}, representing a generalised  singular value decomposition}
\item{labels}{logical if \code{TRUE} plots the labels given in \code{solution[[mod]]["n"]}}
\item{mod}{vectors of the modes numbers to be plotted}
\item{nb1}{number identifying the Principal Tensor to display on the vertical
axe, can be checked using \code{\link{summary.PTAk}} }
\item{nb2}{as nb1 to be displayed on the horizontal axe, if \code{NULL} the horizontal
axe will be used  as \code{Index} (see \code{\link{plot.default}}) }
\item{coefi}{coefficients to multiply components for all modes (not just the one in \code{mod}) for rescaling or changing signs
purposes; each element of the list correspond to \code{nb1} and \code{nb2}
and are vectors of dimentions the tensor order}
\item{xylab}{logical to display axes labels}
\item{ppch}{a vector of length at least \code{length(mod)} used for \code{pch=}}
\item{lengthlabels}{
a number or a vector of numbers of characters in labels to be used for display}
\item{scree}{
logical to display a screeplot of squared singular values as percent of total variation}
\item{ordered}{logical used when displaying the screeplot  with sorted
values (TRUE) or the order is given by output listing from
\code{summary.PTAk}}
\item{nbvs}{a maximum number of singular values to display on the screeplot or a vector of ranks}
\item{max}{is the number of singular values to be considered as giving the perfect fit, NULL is the max possible in x}
\item{rescaled}{boolean to rescale the y axis to 0-100}
\item{RiskJack}{
if not \code{NULL} is an integer,  scree is \code{TRUE} and ordered is \code{TRUE}, plots
on top of the scree plot a Risk plot with \code{nbvs = 1:RiskJack}. It is
possible to use directly the function \code{RiskJackplot}: the default
maximum dimension (argument \code{max}) is  \code{length(solution[[k]][["d"]])}. }
\item{method}{default is \code{""}, a value \code{"FCA"} is to be used only
if \code{solution} is after an FCA with \code{SVDgen}}
\item{ZoomInOut}{ list used as [[1]] for xlim and [[2]] for ylim in xy-plots instead of max and min range}
\item{Zlabels}{ used as labels instead of \code{x[[mod]]$n}, it is a list with the same length as all modes. For example on 3 modes changing the labels of the second mode only will have to set \code{Zlabels=list(NULL,rep("a",length(x[[2]]$n) ), NULL )} }
\item{Zcol}{ list of vectors of colours for Zlabels}
\item{poslab}{integer or vector for 'pos' parameter, position of labels}
\item{signedCTR}{logical to plot signed-CTR instead of coordinates, see \code{\link{CTR}}}
\item{relCTR}{logical if \code{signedCTR = TRUE} use relative CTRs, expected contribution \code{1} if uniform or equal contributions.}
\item{...}{plot arguments can be passed (except \code{xlim, ylim, ylab,pch,xaxt} for component plot,
and \code{xlab, ylab} for screeplot). For example to have normed plot one can use \code{asp=1}}
    }
 \details{
Plot components of one or two Principal Tensors, modes are superposed if more
than one is asked, or gives a screeplot. As it is using \code{plot.default} at
some point some added features can be used in the ... part, especially
\code{xlab= } may be useful when \code{nb2=NULL}. Plots are superposed as they
correspond to the same Principal Tensor and so this gives insight to
interpretation of it, but careful is recommended as only overall
interpretation, once the Principal Tensor has been rebuilt mentally
(\emph{i.e.} product of signs ...) to work out oppositions or associations. The
risk plot on top of a screeplot is an approximation of the Jacknife estimate of
the MSE in the choice of number of dimensions (see Besse et al.(1997)). }

\author{Didier G. Leibovici}
\references{
Besse, P Cardot, H and Ferraty, F (1997) \emph{Simultaneous
non-parametric regressions of unbalanced longitudinal
data}. Computational Statistics and Data Analysis, 24:255-270.

Leibovici D (2000) \emph{Multiway Multidimensional Analysis for Pharmaco-EEG
Studies}. \url{https://www.researchgate.net/publication/216807619_Multiway_Multidimensional_Analysis_for_Pharmaco-EEG_Studies}
             }
\note{This function is used all for \code{FCAk}, and
  \code{CANDPARA}, \code{PCAn} objjects notheless for this
  last object other interesting plots known as jointplots have not been implemented.}

\seealso{ \code{\link{PTAk}}, \code{\link{PTA3}}, \code{\link{CTR}},
\code{\link{FCAk}},\code{\link{SVDgen}}}
\examples{
#  see the demo function   source(paste(R.home(),"/ library/PTAk/demo/PTA3.R",sep=""));
# or    source(paste(R.home(),"/ library/PTAk/demo/PTAk.R",sep=""));
 
 # demo.PTA3()
 }

\keyword{hplot}
\keyword{multivariate}
