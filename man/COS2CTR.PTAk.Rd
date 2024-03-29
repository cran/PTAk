\encoding{utf-8}
\name{COS2-CTR}
\alias{COS2}
\alias{CTR}
\title{Interpretation summaries}
\description{
After a \code{\link{FCA2}}, a \code{\link{SVDgen}}, a \code{\link{FCAk}} or a \code{\link{PTAk}} computes the traditional guides for interpretations used in PCA and correspondence analysis: COS2 or the percentage of variability rebuilt by the component and CTR or the amount of contribution towards that component.
}
\usage{
COS2(solu, mod=1, solnbs=2:4)
CTR(solu, mod=1, solnbs=1:4, signed = TRUE, mil = TRUE)
}
\arguments{
   \item{solu}{an object inheriting from class \code{PTAk}, representing a generalised  singular value decomposition}
   \item{mod}{an integer representing the mode number entry, 1 is row, 2 columns, ...}
   \item{solnbs}{a vector of integers representing the tensor numbers in the listing summary}
   \item{signed}{logical to use signed-CTR from affect the sign of corresponding value in \code{solu[[mod]]$v[,c ]}, c defined by solnbs.}
   \item{mil}{logical}
    }

\details{
	Classical measures helping to interpret the plots in PCA, FCA and in PTAk as well.
	The sum of the COS2 across all the components needed to rebuild fully the tensor analysed) would make 1000 and the sum pf the CTR across the entry mode would be 1000.
}
 \value{a matrix whose columns are the COS2 or CTR  as per thousands (‰) for the mode considered }
 

\references{
Escoufier Y (1985) \emph{L'Analyse des correspondances : ses
\enc{propriétés}{propretes} et ses extensions}. ISI 45th session Amsterdam.

Leibovici D(1993) \emph{Facteurs \enc{à}{a} Mesures \enc{Répétées}{Repetees} et Analyses Factorielles :
applications \enc{à}{a} un suivi \enc{Epidémiologique}{Epidemiologique}}. \enc{Université}{Universite} de Montpellier
II. PhD Thesis in \enc{Mathématiques}{Mathematiques} et Applications (Biostatistiques).

Leibovici DG (2010) \emph{Spatio-temporal Multiway Decomposition using Principal Tensor Analysis on k-modes:the R package \pkg{PTAk}}. Journal of Statistical Software, 34(10), 1-34. \doi{10.18637/jss.v034.i10}
}
\author{Didier G. Leibovici}
\seealso{\code{\link{PTAk}}, \code{\link{FCA2}}, \code{\link{FCAk}}, \code{\link{summary.FCAk}}, \code{\link{plot.PTAk}} }

\examples{
 data(crimerate)
  cri.FCA2 <- FCA2(crimerate)
   summary(cri.FCA2)
    plot(cri.FCA2, mod = c(1,2), nb1 = 2, nb2 = 3) # unscaled
    plot(cri.FCA2, mod = c(1,2), nb1 = 2, nb2 = 3, coefi = 
    	list(c(0.130787,0.130787),c(0.104359,0.104359)) ) # symmetric-map biplot
   CTR(cri.FCA2, mod = 1, solnbs = 2:4)
   CTR(cri.FCA2, mod = 2, solnbs = 2:4)
 COS2(cri.FCA2, mod = 2, solnbs = 2:4)
 ##### useful fonctions
 ##selecting and sorting out dimensions positive and negative sides
 "ctrcos2" <-function(Ta, mod=1, dim=2, NegPos=TRUE, select=c("avg",12,"none"),nbdig=2, 
 			cos2min=333){
	dim=c(dim, dim+1)
	ctr=CTR(Ta,mod=mod,solnbs=dim);cos2=COS2(Ta,mod=mod,solnbs=dim)
	val=round(Ta[[mod]]$v[dim,][1,],digits=nbdig)
	oo=order(ctr[,1],decreasing=TRUE)
	if(NegPos)oo=order(val,decreasing=TRUE)
	out=cbind(ctr[oo,1],cos2[oo,1],val[oo])
	colnames(out)=c("ctr","cos2",paste0("dim",dim[1]))
	if(select[1]=="none") sout=0
	if(select[1]=="avg") sout= 1000/length(val)
	if(is.numeric(select[1])) sout=select[1]
return(out[ out[,1]>=sout | out[,2]>=cos2min, ])	 
}#ctrcos2

## plot ctr cos2
"plotctrcos2"<-function(sol,mod12=c(1,2),dim=2, ratio2avg=TRUE, col=c(1,2),pch=c("-","°"), 
				posi=c(2,3),reposi=TRUE, cos2min=333,select="avg",...){
	### ctrcos2 ini Ta,mod=1,dim=2, NegPos=TRUE, select=c("avg",12,"none"), nbdig=2, cos2min=333
	pre<-function(mod=1,soldim=2, ...){
	diim=length(sol[[mod]]$v[1,])
	ctrov= ctrcos2(sol,dim= soldim,mod=mod,...)	
	x=ctrov[,1]*sign(ctrov[,3])
	if(ratio2avg)x=round(x/(1000/diim),2)
	y=ctrov[,2]
	lab=rownames(ctrov)
	len=dim(ctrov)[1]
	return(list("x"=x,"y"=y,"len"=len,"lab"=lab))
	}
	if(length(col)<length(mod12))col=rep(col,length(mod12))
	if(length(pch)<length(mod12))pch=rep(pch,length(mod12))
	x=NULL;y=NULL;coul=NULL;pchl=NULL;lab=NULL;poslab=NULL
	for(m in mod12){
		prep=pre(mod=m,soldim=dim,...)
		x=c(x,prep$x);y=c(y,prep$y);
        coul=c(coul,rep(col[m],prep$len));pchl=c(pchl,rep(pch[m],prep$len))
		repos=rep(posi[m],prep$len); if(reposi)repos=sample(1:4,prep$len, replace=TRUE)
		lab=c(lab,prep$lab);poslab=c(poslab,repos)
	}
	summsol=summary(sol)
	
	if(match("FCA2" ,class(sol),nomatch=0)>0) xlabe=paste0( "Global pct ",
	round(summsol[dim,4],2), " FCA pct",round(summsol[dim,5],2)) else 
	xlabe=paste0( " local pxt",round(summsol[dim,4],2), " Global pct", round(summsol[dim,5],2))
	dimi=paste0("dim",dim)
	if(ratio2avg)ctrlab="CTR (signed ctr /(uniform ctr))" else ctrlab="CTR (signed)"
	if(!is.null(cos2min))cos2lab=paste("COS2 (> ",cos2min,")")else cos2lab= "COS2"
	plot(x,y,xlab=ctrlab, main=paste(dimi, xlabe ),ylab= cos2lab,col=coul,
	             pch=pchl,ylim=c(min(y),1050),xlim=c(min(x-0.5),max(x+0.5)))
	abline(v=0,col=4,lty=2)
	abline(v=1,col=3,lty=2)
	abline(v=-1,col=3,lty=2)
	text(x,y,lab,pos=poslab, col=coul)
return(cbind(x,y,lab,coul,pchl,poslab))
}#plotctrcos2

ctrcos2(cri.FCA2,mod=1)
ctrcos2(cri.FCA2,mod=2)
plotctrcos2(cri.FCA2)
 
}
\keyword{models}
\keyword{principal components analysis}
\keyword{array }  \keyword{algebra}
\keyword{multivariate}
