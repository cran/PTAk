\name{CauRuimet}
\alias{CauRuimet}
\title{Robust estimation of within group covariance}
\description{
Gives a robust estimate of an unknown within group covariance, aiming
either to look for dense groups or to sparse groups (outliers)
according to \emph{local variance and weighting function} choice.
 }
\usage{
 CauRuimet(X,ker=1,m0=1,withingroup=TRUE,
              loc=substitute(apply(X,2,mean,trim=.1)),matrixmethod=TRUE)

        }
 \arguments{
   \item{X}{matrix}
   \item{ker}{either numerical or a function:
if numerical the weighting function is \eqn{e^{(-ker \;t)}}, otherwise
\cr ker=function(t)\{return(expression)\} is a positive decreasing function. }
   \item{m0}{
is a graph of neighbourhood or another proximity matrix, the hadamard
product of the proximities will be operated}
   \item{withingroup}{
logical,if \code{TRUE} the aim is to give a robust estimate for dense
groups, if \code{FALSE} the aim is to give a robust estimate for outliers}
  \item{loc}{
a vector of locations or a function using mean, median, to give an
estimate of it}
\item{matrixmethod}{if \code{TRUE} (only with \code{withingroup}) uses some matrix computation rather
  than double looping as suggests the formula below }
 }

 \details{
When withingroup is \code{TRUE}, local(defined by the weighting) variance
formula is used, aiming at finding dense groups:
\deqn{W_g=\frac{\sum_{i=1}^{n-1}\sum_{j=i+1}^n
ker(d^2_{S^-}(X_i,X_j))(X_i-X_j)'(X_i-X_j)}{\sum_{i=1}^{n-1}\sum_{j=i+1}^n
ker(d^2_{S^-}(X_i,X_j))}} where \eqn{d^2_{S^-}( . , .)} is the
squared euclidian distance with \eqn{S^-} the inverse of a robust
sample covariance (i.e. using \code{loc} instead of the mean) ; if
\code{FALSE} weighted global variance is used:
 \deqn{W_o=\frac{\sum_{i=1}^nker(d^2_{S^-}(X_i,\tilde{X}))(X_i-\tilde{X})'(X_i-\tilde{X})}
 {\sum_{i=1}^n  ker(d^2_{S^-}(X_i,\tilde{X}))}}
where \eqn{\tilde{X}} is the vector \code{loc}.
\cr
If \code{m0} is a graph of neighbourhood and ker is the function
returning 1 (no proximity due to distance is used) the function will
return (when \code{withingroup=TRUE}) the \emph{local
variance-covariance} matrix as define in Lebart(1969). }
 \value{a matrix}
 \references{
Caussinus, H and Ruiz, A (1990) \emph{Interesting Projections of
Multidimensional Data by Means of Generalized Principal Components
Analysis}. COMPSTAT90, Physica-Verlag, Heidelberg,121-126.

Faraj, A (1994) \emph{Interpretation tools for Generalized
Discriminant Analysis}.In: New Approches in Classification and Data
Analysis, Springer-Verlag, 286-291, Heidelberg.

Lebart, L (1969) \emph{Analyse statistique de la
contiguit�e}.Publication de l'Institut de Statistiques Universitaire
de Paris, XVIII,81-112.

            }

\note{
As mentioned by Caussinus and Ruiz a good strategy to reveal dense
groups with generalised PCA would be to reveal outliers first using
the metric \eqn{W_o^{-1}} and remove them before using the metric
\eqn{W_g^{-1}}.}
\author{Didier Leibovici \email{didier@fmrib.ox.ac.uk}}
\seealso{\code{\link{SVDgen}}}
\examples{

data(iris)
iris2 <- as.matrix(iris[,1:4])
dimnames(iris2)[[1]] <- iris[,5]

D2 <- CauRuimet(iris2,ker=6,withingroup=TRUE)
D2 <- Powmat(D2,(-1))
iris2 <- sweep(iris2,2, apply(iris2,2,mean))
res <- SVDgen(iris2,D2=D2,D1=1)
plot(res,nb1=1,nb2=2,cex=0.5)
summary(res,testvar=0)

 # the same in a demo function
demo.SVD.CauRui(ker=4,withingroup=TRUE,openX11s=FALSE)
demo.SVD.CauRui(ker=0.15,withingroup=FALSE,openX11s=FALSE)

}
\keyword{robust}
\keyword{multivariate}